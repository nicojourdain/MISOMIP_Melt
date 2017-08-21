## MISOMIP Melt

Elmer Solvers to feed Elmer with NEMO's melt rates in the context of MISOMIP.

* MISOMIP\_Melt.F90 : Elmer solver to interpolate NEMO's melt rates onto Elmer/Ice's grid (regular NEMO grid only).
* MISOMIP\_Melt\_Consv.F90 : Same as MISOMIP\_Melt.F90 but ensures conservative melt exchange (recommended).
* MISOMIP\_Melt\_Consv\_Evolv.F90 : Same as MISOMIP\_Melt\_Consv.F90 but with extrapolation at GL (beta version).
* ThicknessSolver\_Alpha.F90 : Nacho's version of ThicknessSolver.F90 (see his PhD thesis).
