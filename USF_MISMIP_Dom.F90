
FUNCTION bedrock(x,y) RESULT(bed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE
   REAL(KIND=dp) :: x, y, Bx, By, bed

!!!Parameters MISMIP+!!!!!

   REAL(KIND=dp) :: Bmax= 720.0
   REAL(KIND=dp) :: B0 = -150.0
   REAL(KIND=dp) :: B2 = -728.8
   REAL(KIND=dp) :: B4 = 343.91
   REAL(KIND=dp) :: B6 = -50.57
   REAL(KIND=dp) :: xchar = 300.0e3
   REAL(KIND=dp) :: fc = 4.0e3
   REAL(KIND=dp) :: dc = 500.0
   REAL(KIND=dp) :: wc = 24.0e3
   REAL(KIND=dp) :: xcalv = 640.0e3
   REAL(KIND=dp) :: Ly = 80.0e3
   REAL(KIND=dp) :: Lx = 640.0e3

   Bx = B0 &
        + B2*(x/xchar)*(x/xchar) &
        + B4*(x/xchar)*(x/xchar)*(x/xchar)*(x/xchar) &
        + B6*(x/xchar)*(x/xchar)*(x/xchar)*(x/xchar)*(x/xchar)*(x/xchar)

   By = dc/ ( 1 + exp(-2 *( y - (Ly/2) - wc) / fc)) &
        +dc/ ( 1 + exp( 2 * (y -(Ly/2) + wc) / fc))

   IF ((BX+BY) .gt. -1.0 * Bmax) THEN
        bed = Bx + By
   ELSE
        bed = -1.0 * BMax
   END IF

END FUNCTION bedrock

FUNCTION bedIni ( Model, nodenumber, x) RESULT(Zbed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x(2), Bx, By, Zbed, bedrock       

   Zbed = bedrock(x(1),x(2))
   
END FUNCTION bedIni

FUNCTION ZtIni ( Model, nodenumber, x) RESULT(Ztop)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x(2), Bx, By, Zbed, Ztop, IceFloatingTop, IceFloatingBot, bedrock
   LOGICAL :: FirstTime=.True.

   REAL(KIND=dp) :: rhoi = 928.0   
   REAL(KIND=dp) :: rhow = 1028.0

   REAL(KIND=dp) :: initThick  = 100.0

   IceFloatingTop = initThick-100*rhoi/rhow
   IceFloatingBot = IceFloatingTop-initThick
   
   Zbed = bedrock(x(1),x(2))
  
   IF ( Zbed .lt. IceFloatingBot) THEN
	Ztop = IceFloatingTop
   ELSE
	Ztop = initThick + Zbed
   END IF

END FUNCTION ZtIni

FUNCTION ZbIni ( Model, nodenumber, x) RESULT(Zb)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x(2), Bx, By, Zbed, Ztop, Zb, IceFloatingTop, IceFloatingBot, bedrock
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)
   LOGICAL :: FirstTime=.True.

   REAL(KIND=dp) :: rhoi = 928.0
   REAL(KIND=dp) :: rhow = 1028.0

   REAL(KIND=dp) :: initThick  = 100.0

   IceFloatingTop = initThick-100*rhoi/rhow
   IceFloatingBot = IceFloatingTop-initThick

   Zbed = bedrock(x(1),x(2))

   IF ( Zbed .lt. IceFloatingBot) THEN
        Ztop = IceFloatingTop
   ELSE
        Ztop = initThick + Zbed
   END IF


   Zb = Ztop - initThick

END FUNCTION ZbIni
