module subs
contains
subroutine check(status)
    USE NETCDF
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

end subroutine check

SUBROUTINE BilinealInterp(xP,yP,NETCDFValues,meltInterp, dxarg, dyarg, xInitarg, yInitarg)

        USE types
        USE CoordinateSystems
        USE SolverUtils
        USE ElementDescription
        USE DefUtils

        REAL(KIND=dp), INTENT(IN):: xP,yP
        REAL(Kind=dp), INTENT(IN) :: NETCDFValues(:,:)
        REAL(KIND=dp), INTENT(OUT) :: meltInterp
        REAL(KIND=dp) :: iFract, jFract, melt11, melt12, melt21, melt22
        INTEGER :: iIndex, jIndex
        REAL(KIND=dp) :: DX , DY, xInit, yInit
        REAL(KIND=dp), INTENT(IN), OPTIONAL ::  DXarg , DYarg, xInitarg, yInitarg


        !Default values for MISOMIP
        if( .not. present(DXarg)) then
            DX=2000.0
        else
            DX = DXarg
        end if

        if( .not. present(DYarg)) then
            DY=2000.0
        else
            DY = DYarg
        end if

        if( .not. present(xInitarg)) then
            xInit = 319000.0
        else 
            xInit = xInitarg
        end if

        if( .not. present(yInitarg)) then
            yInit = -1000.0
        else
            yInit = yInitarg
        end if

        !Where the node is found in the NetCDF grid?

        iIndex = int(floor((xP - xInit ) / DX)+1)
        iFract = (xP - xInit) - (iIndex-1) * DX

        jIndex = int(floor((yP - (yInit)) / DY)+1)
        jFract = (yP - (yInit)) - (jIndex-1) * DY

        if (iIndex .gt. 0) then
                melt11 = NETCDFValues(iIndex,jIndex)
                melt12 = NETCDFValues(iIndex,jIndex+1)
                melt21 = NETCDFValues(iIndex+1,jIndex)
                melt22 = NETCDFValues(iIndex+1,jIndex+1)
                meltInterp = 1/(DX*DY) * ( &
                        melt11*(DX-iFract)*(DY-jFract) + &
                        melt21*(iFract)*(DY-jFract) + &
                        melt12*(DX-iFract)*(jFract) + &
                        melt22*(iFract)*(jFract) &
                        )
        else    
                meltInterp = 0.0_dp
        end if


End subroutine BilinealInterp
end module subs

SUBROUTINE MISOMIP_Melt( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils
  USE subs
  USE NETCDF

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Solver_t),POINTER :: PSolver
  TYPE(Variable_t),POINTER :: MeltVar=>NULL(), GMVar=>NULL()
 
  INTEGER , POINTER :: MeltPerm(:), GMPerm(:)
  REAL(KIND=dp) , POINTER :: Melt(:),GM(:)
  REAL(KIND=dp) , POINTER ::DATAPointer(:,:)
  INTEGER :: NMax, ncid, node, ncidDraft
  REAL(KIND=dp) ::  xP , yP, meltInt

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL,SAVE :: ExtrudedMesh=.False.
  LOGICAL :: Found, Got

  CHARACTER(len = 200) :: FILE_NAME 
  CHARACTER(len = 200) :: FILE_NAME_DRAFT
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='InitMELTMISOMIP'

  CHARACTER(len = 200), parameter :: meltname='fwfisf', xname='x', yname='y'
  REAL(KIND=dp), allocatable, target :: meltvarNC_Av(:,:), meltvarNC(:,:,:), xVarNC(:), yVarNC(:)
  INTEGER :: varXid, varYid, varid, dimid1, dimid2, dimid3, lenX, lenY, lenTime, status1, res

  REAL(KIND=dp) :: x_NC_Init, x_NC_Fin, y_NC_Init, y_NC_Fin, x_NC_Res, y_NC_Res, facunits


!------------------------------------------------------------------------------
  Mesh => Model % Mesh

  NMax = Solver % Mesh % NumberOfNodes

!!! get required variables Zb,Zs,H

  MeltVar => VariableGet( Model % Mesh % Variables, 'Melt')
  IF (.NOT.ASSOCIATED(MeltVar)) THEN
     Message='Melt not found'
     CALL FATAL(SolverName,Message)
  END IF
  GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
  IF (.NOT.ASSOCIATED(GMVar)) THEN
     Message='GroundedMask not found'
     CALL FATAL(SolverName,Message)
  END IF

  FILE_NAME_DRAFT = GetString(Solver % Values,'Draft file',Got)

  IF (.NOT. Got) then
     Message='Draft File not found'
     CALL FATAL(SolverName,Message)
  END IF

  FILE_NAME = GetString(Solver % Values,'Melt rates file',Got)

  IF (.NOT. Got) then
     Message='Melt Rates File not found'
     CALL FATAL(SolverName,Message)
  END IF

  facunits = GetCReal( Model % Constants, 'Factor melt units', Found )
  IF (.NOT.Found) &
             CALL FATAL(SolverName,'<Factor melt units> not found')

  MeltPerm => MeltVar % Perm
  Melt => MeltVar % Values

  GMPerm => GMVar % Perm
  GM => GMVar % Values

  ! GET NTCDF DIMENSIONS FOR ALLOCATION
  CALL check(nf90_open(FILE_NAME,NF90_NOWRITE,ncid))
  CALL check(nf90_open(FILE_NAME_DRAFT,NF90_NOWRITE,ncidDraft))
  status1=nf90_inq_dimid(ncid,"x",dimid1)
  status1=nf90_inquire_dimension(ncid,dimid1,len=lenX)

  status1=nf90_inq_dimid(ncid,"y",dimid2)
  status1=nf90_inquire_dimension(ncid,dimid2,len=lenY)

  status1=nf90_inq_dimid(ncid,"time_counter",dimid3)
  status1=nf90_inquire_dimension(ncid,dimid3,len=lenTime)

  allocate(meltvarNC(lenX,lenY,lenTime))
  allocate(meltvarNC_Av(lenX,lenY))
  allocate(xVarNC(lenX))
  allocate(yVarNC(lenY))


  !GET Variables

  status1=nf90_inq_varid(ncid,meltname,varid)

  status1=nf90_get_var(ncid,varid,meltVarNC)

  status1=nf90_inq_varid(ncidDraft,xname,varXid)

  status1=nf90_get_var(ncidDraft,varXid,xVarNC)

  status1=nf90_inq_varid(ncidDraft,yname,varYid)

  status1=nf90_get_var(ncidDraft,varYid,yVarNC)

  x_NC_Init = MINVAL(xVarNC)
  x_NC_Fin = MAXVAL(xVarNC)
  y_NC_Init = MINVAL(yVarNC)
  y_NC_Fin = MAXVAL(yVarNC)

  x_NC_Res = xVarNC(2)-xVarNC(1)
  y_NC_Res = yVarNC(2)-yVarNC(1)

  MeltVarNC_Av = SUM(meltvarNC,dim=3)/lenTime

  DO node=1, nMax 
        xP =  Mesh % Nodes % x(node)
        yP =  Mesh % Nodes % y(node)

        if (xP .gt. x_NC_Fin .or. xP .lt. x_NC_Init) then
                Melt(MeltPerm(node)) = 0.0_dp   
                cycle
        end if

        if (yP .gt. y_NC_Fin .or. yP .lt. y_NC_Init) then
                Melt(MeltPerm(node)) = 0.0_dp
                cycle
        end if

        if (GM(GMPerm(node)) .lt. 0.5) then
                CALL BiLinealInterp(xP,yP,meltvarNC_Av, meltInT, x_NC_Res, y_NC_Res, x_NC_Init, y_NC_Init)
                Melt(MeltPerm(node)) = meltInt * facunits ! from kg/m^2/s to meter of ice / yr
        else
                Melt(MeltPerm(node)) = 0.0_dp
        end if
  end do

  CALL check(nf90_close(ncid))
  CALL check(nf90_close(ncidDraft))

END SUBROUTINE MISOMIP_Melt

