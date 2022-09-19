!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   Fig6 :   AUTO file to reproduce Figure 6 of the article
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!--------- ---- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)

  DOUBLE PRECISION R, V, S, K, Q
  DOUBLE PRECISION DELTA, ETABAR, J, TAU, EPS, T, PI
  
       ! Define the state variables
       R = U(1)
       V = U(2)
       S = U(3)
       K = U(4)
       Q = U(5)

       
       ! Define the system parameters
       DELTA  = PAR(1)
       ETABAR = PAR(2)
       J      = PAR(3)
       TAU    = PAR(4)
       EPS    = PAR(5)
       T      = PAR(11)

       ! Define PI
       PI = 4.0d0 * ATAN(1.0d0)

       ! Define the right-hand sides
       F(1) =  T*(DELTA/PI + 2*R*V)
       F(2) =  T*(V**2 - PI**2*R**2 + J*S + K)
       F(3) =  T*(-S+R)/TAU
       F(4) =  T*EPS*Q
       F(5) = -T*EPS*(K-ETABAR)
       
END SUBROUTINE FUNC
!---------------------------------------------------------------------- 

SUBROUTINE STPNT(NDIM,U,PAR,T)

!---------------------------------------------------------------------- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  PAR(1)  =   1.00d0   !  PARAMETER DELTA
  PAR(2)  =  -6.50d0   !  PARAMETER ETABAR
  PAR(3)  =  15.00d0   !  PARAMETER J
  PAR(4)  =   0.02d0   !  PARAMETER TAU
  PAR(5)  =   0.05d0   !  PARAMETER EPS
  PAR(6)  =   3.00d0   !  Q(0)=A
  PAR(11) = 125.6637d0 !  PARAMETER T
  
END SUBROUTINE STPNT

SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP
      !INTEGER NDX,NCOL,NTST

      ! Set PAR(7) equal to the maximum of U(1) - minimum of U(1)
      PAR(7)=GETP('MAX',1,U)-GETP('MIN',1,U)
      

END SUBROUTINE PVLS

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

      ! DEFINE BOUNDARY CONDITIONS

      FB(1) = U1(1) - U0(1)
      FB(2) = U1(2) - U0(2)
      FB(3) = U1(3) - U0(3)
      FB(4) = U0(4) - PAR(2)   ! K(0)=ETABAR
      FB(5) = U1(4) - PAR(2)   ! K(1)=ETABAR
      FB(6) = U0(5) - PAR(6)   ! Q(0)=A
      
END SUBROUTINE BCND

SUBROUTINE ICND
END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT
!---------------------------------------------------------------------- 
