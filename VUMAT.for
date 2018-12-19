!     # of PROPS = 19
!     # of STATE = 24
!      
!     THE STATE VARIABLES ARE STORED AS:
!
!     STATE(*, 1) = void ratio
!     STATE(*, 2) = plastic strain 11
!     STATE(*, 3) = plastic strain 22
!     STATE(*, 4) = plastic strain 33
!     STATE(*, 5) = plastic strain 12
!     STATE(*, 6) = plastic strain 23
!     STATE(*, 7) = plastic strain 31
!     STATE(*, 8) = p_s
!     STATE(*, 9) = beta
!     STATE(*, 10) = po
!     STATE(*, 11) = qo
!     STATE(*, 12) = p
!     STATE(*, 13) = q
!     STATE(*, 14) = pimg
!     STATE(*, 15) = qimg
      

      SUBROUTINE VUMAT(
      !READ ONLY -
     1 NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3 PROPS, DENSITY, STRAININC, RELSPININC,
     4 TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5 STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6 TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
      !WRITE ONLY -
     7 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)

      !INCLUDE 'VABA_PARAM.INC'

      REAL*8 PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     1 CHARLENGTH(NBLOCK), STRAININC(NBLOCK, NDIR+NSHR),
     2 RELSPININC(NBLOCK, NSHR), TEMPOLD(NBLOCK),
     3 STRETCHOLD(NBLOCK, NDIR+NSHR),DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     4 FIELDOLD(NBLOCK, NFIELDV), STRESSOLD(NBLOCK, NDIR+NSHR),
     5 STATEOLD(NBLOCK, NSTATEV), ENERINTERNOLD(NBLOCK),
     6 ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     7 STRETCHNEW(NBLOCK, NDIR+NSHR),DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     8 FIELDNEW(NBLOCK, NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR),
     9 STATENEW(NBLOCK, NSTATEV), ENERINTERNNEW(NBLOCK),
     1 ENERINELASNEW(NBLOCK)

      CHARACTER*80 CMNAME 
      
      ! DEFORMATION GRADIENT
      REAL*8 DefGrad_old(1:3,1:3)
      REAL*8 DefGrad_new(1:3,1:3)
      REAL*8 DefGrad_mid(1:3,1:3)
      
      !------------------
      ! Model Parameters
      !------------------

      ! CRITICAL STATE STRESS RATIO
      REAL*8 m

      ! NORMARLY CONSOLIDATED LINE AND RECOMPRESSION LINE
      REAL*8 lambda, kappa

      ! HARDENING MOUDULUS
      REAL*8 ho, mu

      ! ELASTIC MODULI: G, K, AND nu
      REAL*8 G, K, nu
      
      ! TAIL LENGTH
      REAL*8 lt

      ! SIZE AND ROTATION OF YIELD SURFACE
      REAL*8 ps_ini, beta_ini
      
      ! INITIAL VOID RATIO
      REAL*8 voidratio_ini

      !--------------------------------
      ! PROPERTIES: STRESS-INTEGRATION
      !--------------------------------

      ! MACHINE EPSILON
      REAL*8 EPS

      ! TOLERANCE FOR STRESS INTEGRATION ALGORITHM, YIELD FURFACE, BOUNDING SURFACE
      REAL*8 STOL, FTOL, LTOL

      ! OLD STATES (BEFORE INCREMENT)
      REAL*8 voidratio_old    ! VOID RATIO
      REAL*8 stress_old(1:3,1:3)  ! STRESS
      REAL*8 Ep_old(1:3,1:3)      ! PLASTIC STRAIN      
      REAL*8 ps_old, beta_old ! INTERNAL VARIABLES
      REAL*8 po_old, qo_old   ! PROJECTION CENTER
      REAL*8 p_old, q_old     ! STRESS INVARIANTS
      REAL*8 Epv_old, Epq_old ! PLASTIC STRAIN INVARIANTS
      
      ! NEW STATE (OR UPDATED STATE)
      REAL*8 voidratio_new        ! VOID RATIO
      REAL*8 stress_new(1:3,1:3)  ! STRESS
      REAL*8 Ep_new(1:3,1:3)      ! PLASTIC STRAIN
      REAL*8 ps_new, beta_new ! INTERNAL VARIABLES
      REAL*8 po_new, qo_new   ! PROJECTION CENTER      
      REAL*8 p_new, q_new     ! STRESS INVARIANTS
      REAL*8 Epv_new, Epq_new ! PLASTIC STRAIN INVARIANTS

      ! STRAIN INCREMENTS
      REAL*8 dE(1:3,1:3), dEv, dEq    ! VOLUMETRIC AND DEVIATORIC PART OF STRAIN INCREMENT
      REAL*8 dEe(1:3,1:3), dEev, dEeq ! VOLUMETRIC AND DEVIATORIC PART OF ELASTICSTRAIN INCREMENT
      REAL*8 dEp(1:3,1:3), dEpv, dEpq ! VOLUMETRIC AND DEVIATORIC PART OF PLASTIC STRAIN INCREMENT
      
      ! IMAGE POINT
      REAL*8 pimg, qimg
      
      ! GRADIENT OF BOUNDING SURFACE (OR PLASTIC POTENTIAL)
      REAL*8 dfdp, dfdq, dfdps, dfdbeta
      
      ! HARDENING MODULUS
      REAL*8 hb, h
      
      !---------------------------------------------
      ! IDENTITY, ZERO MATRICES, AND PARAMETERS
      !---------------------------------------------

      REAL*8 IDENTITY(1:3,1:3), ZEROS(1:3,1:3)
      REAL*8 ONE_OVER_THREE
      REAL*8 TWO_OVER_THREE
      REAL*8 SQRT_TWO
      REAL*8 SQRT_THREE
      REAL*8 SQRT_SIX
      REAL*8 SQRT_ONE_OVER_THREE
      REAL*8 SQRT_TWO_OVER_THREE
      REAL*8 SQRT_THREE_OVER_TWO

      PARAMETER (PI =                   3.14159265358979D0,
     1           ONE_OVER_THREE =       3.33333333333333D-1,
     2           TWO_OVER_THREE =       6.66666666666666D-1,
     3           SQRT_TWO =             1.41421356237310D0,
     4           SQRT_THREE =           1.73205080756888D0,
     5           SQRT_SIX =             2.44948974278318D0,
     6           SQRT_ONE_OVER_THREE =  5.77350269189626D-1,
     7           SQRT_TWO_OVER_THREE =  8.16496580927726D-1,
     8           SQRT_THREE_OVER_TWO =  1.22474487139159D0)

      ! SETTING IDENTITY MATRIX AND ZERO MATRIX TO USE LATER
      IDENTITY(1,1)=1.0D0
      IDENTITY(2,2)=1.0D0
      IDENTITY(3,3)=1.0D0
      IDENTITY(1,2)=0.0D0
      IDENTITY(2,1)=0.0D0
      IDENTITY(1,3)=0.0D0
      IDENTITY(3,1)=0.0D0
      IDENTITY(2,3)=0.0D0
      IDENTITY(3,2)=0.0D0
      
      ZEROS(1,1)=0.0D0
      ZEROS(2,2)=0.0D0
      ZEROS(3,3)=0.0D0
      ZEROS(1,2)=0.0D0
      ZEROS(2,1)=0.0D0
      ZEROS(1,3)=0.0D0
      ZEROS(3,1)=0.0D0
      ZEROS(2,3)=0.0D0
      ZEROS(3,2)=0.0D0

      ! PROPERTIES - SET MODEL PROPERTIES
      nu = PROPS(1)
      m = PROPS(2)
      lambda = PROPS(3)
      kappa = PROPS(4)
      ho = PROPS(5)
      mu = PROPS(6)
      Lt = PROPS(7)
      ps_ini = PROPS(8)
      beta_ini = PROPS(9)
      voidratio_ini = PROPS(10)
      
      !PRINT*, 'G:', G
      !PRINT*, 'nu:', nu
      !PRINT*, 'm:', m
      !PRINT*, 'lambda:', lambda
      !PRINT*, 'kappa:', kappa
      !PRINT*, 'ho:', ho
      !PRINT*, 'mu:', mu
      !PRINT*, 'Lt:', Lt      
      
      ! SET TOLERANCES
      STOL = 1.D-3
      FTOL = 1.D-8
      LTOL = 1.D-3

      ! SET MACHINE EPSILON
      ! TODO: PRE-CALCULATION OF EPS, IT TAKES LONG TIME!!!!
      EPS = 1.0D0
      DO WHILE((1.0D0+EPS) .NE. 1.0D0)
          EPS = EPS/2.0D0
      END DO

      !-------------------
      ! STRESS UPDATE
      !-------------------

      ! FOR ALL QUADRATURE (OR MATERIAL) POINTS
      DO I = 1, NBLOCK

          ! SET OLD STRESS
          stress_old(1,1) = -STRESSOLD(I,1)
          stress_old(2,2) = -STRESSOLD(I,2)
          stress_old(3,3) = -STRESSOLD(I,3)
          stress_old(1,2) = -STRESSOLD(I,4)
          stress_old(2,1) = -STRESSOLD(I,4)
          IF (NSHR == 3) THEN
              stress_old(2,3) = -STRESSOLD(I,5)
              stress_old(3,2) = -STRESSOLD(I,5)
              stress_old(1,3) = -STRESSOLD(I,6)
              stress_old(3,1) = -STRESSOLD(I,6)
          ELSE
              stress_old(2,3) = 0.D0
              stress_old(3,2) = 0.D0
              stress_old(1,3) = 0.D0
              stress_old(3,1) = 0.D0
          END IF
        
          ! SET STRAIN INCREMENT
          dE(1,1) = -STRAININC(I,1)
          dE(2,2) = -STRAININC(I,2)
          dE(3,3) = -STRAININC(I,3)
          dE(1,2) = -STRAININC(I,4)
          dE(2,1) = -STRAININC(I,4)
          IF (NSHR == 3) THEN
              dE(2,3) = -STRAININC(I,5)
              dE(3,2) = -STRAININC(I,5)
              dE(1,3) = -STRAININC(I,6)
              dE(3,1) = -STRAININC(I,6)
          ELSE
              dE(2,3) = ZERO
              dE(3,2) = ZERO
              dE(1,3) = ZERO
              dE(3,1) = ZERO
          END IF

          ! SET OLD DEFORMATION GRADIENT
          DefGrad_old(1,1) = DEFGRADOLD(I,1)
          DefGrad_old(2,2) = DEFGRADOLD(I,2)
          DefGrad_old(3,3) = DEFGRADOLD(I,3)
          DefGrad_old(1,2) = DEFGRADOLD(I,4)
          IF(NSHR == 1) THEN
              DefGrad_old(2,1) = DEFGRADOLD(I,5)
              DefGrad_old(2,3) = 0.D0
              DefGrad_old(3,1) = 0.D0
              DefGrad_old(3,2) = 0.D0
              DefGrad_old(1,3) = 0.D0
          ELSE 
              DefGrad_old(2,3) = DEFGRADOLD(I,5)
              DefGrad_old(3,1) = DEFGRADOLD(I,6)
              DefGrad_old(2,1) = DEFGRADOLD(I,7)
              DefGrad_old(3,2) = DEFGRADOLD(I,8)
              DefGrad_old(1,3) = DEFGRADOLD(I,9)
          END IF

          ! SET NEW DEFORMATION GRADIENT
          DefGrad_new(1,1) = DEFGRADNEW(I,1)
          DefGrad_new(2,2) = DEFGRADNEW(I,2)
          DefGrad_new(3,3) = DEFGRADNEW(I,3)
          DefGrad_new(1,2) = DEFGRADNEW(I,4)
          IF(NSHR == 1) THEN
              DefGrad_new(2,1) = DEFGRADNEW(I,5)
              DefGrad_new(2,3) = 0.D0
              DefGrad_new(3,1) = 0.D0
              DefGrad_new(3,2) = 0.D0
              DefGrad_new(1,3) = 0.D0
          ELSE 
              DefGrad_new(2,3) = DEFGRADNEW(I,5)
              DefGrad_new(3,1) = DEFGRADNEW(I,6)
              DefGrad_new(2,1) = DEFGRADNEW(I,7)
              DefGrad_new(3,2) = DEFGRADNEW(I,8)
              DefGrad_new(1,3) = DEFGRADNEW(I,9)
          END IF

          ! SET VOID RATIO AND DENSITY
          voidratio_old = STATEOLD(I,1)

          ! SET PLASTIC STRAIN
          Ep_old(1,1) = STATEOLD(I,2)
          Ep_old(2,2) = STATEOLD(I,3)
          Ep_old(3,3) = STATEOLD(I,4)
          Ep_old(1,2) = STATEOLD(I,5)
          Ep_old(2,1) = STATEOLD(I,5)
          Ep_old(2,3) = STATEOLD(I,6)
          Ep_old(3,2) = STATEOLD(I,6)
          Ep_old(1,3) = STATEOLD(I,7)
          Ep_old(3,1) = STATEOLD(I,7)
        
          ! SET INTERNAL VARIABLES
          ps_old = STATEOLD(I,8)
          beta_old = STATEOLD(I,9)
        
          ! SET PROJECTION CENTER
          po_old = STATEOLD(I,10)
          qo_old = STATEOLD(I,11)

          ! GET STRESS INVARIANTS
          p_old = GetMean(stress_old) 
          q_old = GetMisesStress(stress_old)
                
          ! SET INITIAL CONDITION
          IF (STATEOLD(I,1) .EQ. 0.D0) THEN

              voidratio_old = voidratio_ini            
              ps_old = ps_ini
              beta_old = beta_ini            
              po_old = p_old - 1.D0
              qo_old = q_old
                      
          END IF
          
          p_new = p_old
          q_new = q_old
          po_new = po_old
          qo_new = qo_old
          ps_new = ps_old
          beta_new = beta_old
          Epv_new = Epv_old
          Epq_new = Epq_old

          ! UPDATE VOID RATIO
          voidratio_new=(voidratio_ini+1.0D0)*GetI3(DefGrad_new)-1.D0
          
          ! GET STRAIN INCREMENT INVARIANTS
          dEv = GetI1(dE)
          dEq = GetShearStrain(dE)
          
          !PRINT*, '---'
          !PRINT*, 'PRINT CONDITION BEFORE STRESS INTEGRATION'
          !PRINT*, 'p, q:', p_old, q_old
          !PRINT*, 'po, qo:', po_old, qo_old
          !PRINT*, 'ps, beta:', ps_old, beta_old
          !PRINT*, 'void ratio:', voidratio_new
          !PRINT*, 'dE11, dE22, dE33:', dE(1,1), dE(2,2), dE(3,3)
          !PRINT*, 'dEv, dEq:', dEv, dEq
          !PRINT*, '---'
          
          ! UPDATE ELASTO-PLASTIC
          CALL StateUpdate_ElastoPlastic_TX()

          ! CALCULATE STRESS COMPONENTS FROM NEW p and q
          stress_new(1,1) = p_new + 2.D0/3.D0*q_new
          stress_new(2,2) = p_new - 1.D0/3.D0*q_new
          stress_new(3,3) = p_new - 1.D0/3.D0*q_new
      
          ! CALULATE PLASTICN STARIN COMPONENTS
          
          ! SET STRESS NEW
          STRESSNEW(I,1) = -stress_new(1,1)
          STRESSNEW(I,2) = -stress_new(2,2)
          STRESSNEW(I,3) = -stress_new(3,3)
          STRESSNEW(I,4) = 0.D0
          IF (NSHR == 3) THEN
              STRESSNEW(I,5) = 0.D0
              STRESSNEW(I,6) = 0.D0
          END IF

          ! SET NEW STATE VARIABLES
          STATENEW(I,1) = voidratio_new

          STATENEW(I,2) = Ep_new(1,1)
          STATENEW(I,3) = Ep_new(2,2)
          STATENEW(I,4) = Ep_new(3,3)
          STATENEW(I,5) = Ep_new(1,2)
          STATENEW(I,6) = Ep_new(2,3)
          STATENEW(I,7) = Ep_new(3,1)
        
          ! SET INTERNAL VARIABLES
          STATENEW(I,8) = ps_new
          STATENEW(I,9) = beta_new
        
          ! PROJECTION CENTER
          STATENEW(I,10) = po_new
          STATENEW(I,11) = qo_new
        
          ! PROJECTION CENTER
          STATENEW(I,12) = p_new
          STATENEW(I,13) = q_new
          
          ! PROJECTION CENTER
          STATENEW(I,14) = pimg
          STATENEW(I,15) = qimg
        
      END DO

      RETURN

      ! CONTAINES FOLLOWING FUNCTIONS AND SUBROUTINES
      CONTAINS

      ! UPDATE STATE FOR ELASTOPLASTICITY
      SUBROUTINE StateUpdate_ElastoPlastic_TX()

      IMPLICIT NONE
      
      ! LOCAL VARIABLES
      REAL*8 p_tri, q_tri, po_tri, qo_tri, ps_tri, beta_tri
      REAL*8 Epv_tri, Epq_tri

      ! PSEUDO TIME
      REAL*8 T, dT, qf, qs

      ! BOOL
      LOGICAL success, fail

      ! ERROR ESTIMATORS
      REAL*8 RE, RE_s, RE_ps, RE_beta

      ! STATE 1
      REAL*8 dp1, dq1, dpo1, dqo1, dps1, dbeta1
      REAL*8 dTdEpv1, dTdEpq1

      ! STATE 2
      REAL*8 dp2, dq2, dpo2, dqo2, dps2, dbeta2
      REAL*8 dTdEpv2, dTdEpq2

      ! LOCAL VARIABLES
      REAL*8 dTdEv, dTdEq
      
      ! BS FUNCTION VALUE
      REAL*8 f
        
      ! SET TRIAL STRESS, STATE VARIABLES, DENSITY AND DEFORMATION GRADIENT
      p_tri = p_old
      q_tri = q_old
      po_tri = po_old
      qo_tri = qo_old
      ps_tri = ps_old
      beta_tri = beta_old
      Epv_tri = Epv_old
      Epq_tri = Epq_old
      
      ! INITIALIZE PSEUDO TIME AND ITS INCREMENT
      T = 0.D0
      dT = 1.D0

      ! REPEAT INTEGRATION UNTIL PSEUDO TIME IS EQUAL TO ONE
      DO WHILE(T < 1.D0)

          success = .FALSE.
          fail = .FALSE.
          RE = 0.D0

          ! REPEAT SUB-INCREMENT UNTIL IT IS SUCCESSFUL
          DO WHILE(success == .FALSE.)

              ! COMPUTE STRAIN INCREMENT WITH PSEUDO TIME INCREMENT
              dTdEv = dT * dEv
              dTdEq = dT * dEq
               
              ! STATE INCREMENT 1
              CALL GetStateIncrement(
     1                dTdEv, dTdEq, 
     2                p_tri, q_tri, po_tri, qo_tri, ps_tri, beta_tri,
     3                dTdEpv1, dTdEpq1, 
     4                dp1, dq1, dpo1, dqo1, dps1, dbeta1)
              
              ! UPDATE TRIAL STRESS AND STATE VARIABLES
              p_tri = p_tri + dp1
              q_tri = q_tri + dq1
              Epv_tri = Epv_tri + dTdEpv1
              Epq_tri = Epq_tri + dTdEpq1
              po_tri = po_tri + dpo1
              qo_tri = qo_tri + dqo1
              ps_tri = ps_tri + dps1
              beta_tri = beta_tri + dbeta1
              
              ! STATE 2
              CALL GetStateIncrement(
     1                dTdEv, dTdEq, 
     2                p_tri, q_tri, po_tri, qo_tri, ps_tri, beta_tri,
     3                dTdEpv2, dTdEpq2, 
     4                dp2, dq2, dpo2, dqo2, dps2, dbeta2)
              
              ! ADJUST TRIAL STRESS AND STATE VARIABLES
              p_tri = p_tri - 0.5D0 * dp1 + 0.5D0 * dp2
              q_tri = q_tri - 0.5D0 * dq1 + 0.5D0 * dq2
              Epv_tri = Epv_tri - 0.5D0 * dTdEpv1 + 0.5D0 * dTdEpv2
              Epq_tri = Epq_tri - 0.5D0 * dTdEpq1 + 0.5D0 * dTdEpq2
              po_tri = po_tri - 0.5D0 * dpo1 + 0.5D0 * dpo2
              qo_tri = qo_tri - 0.5D0 * dqo1 + 0.5D0 * dqo2
              ps_tri = ps_tri - 0.5D0 * dps1 + 0.5D0 * dps2
              beta_tri = beta_tri - 0.5D0 * dbeta1 + 0.5D0 * dbeta2

              ! ESTIMATE RELATIVE ERRORS
              RE_s = SQRT((dp1 - dp2)**2.D0 + (dq1 - dq2)**2.D0)/
     1                MAX(SQRT(p_new*p_new + q_new*q_new), EPS)
              RE_ps = ABS(dps1 - dps2)/ MAX(ps_tri, EPS)
              RE_beta = ABS(dbeta1 - dbeta2)
              RE = MAX(RE_s, MAX(RE_ps, RE_beta))
              
              ! CHECK ERROR
              IF(RE <= STOL) THEN
                  ! SUB-STEP IS SUCCESSFUL
                  success = .TRUE.
                  fail = .FALSE.
              ELSE
                  ! SUB-STEP IS FIALED
                  success = .FALSE.
                  fail = .TRUE.

                  qf = MAX(0.9D0*SQRT(STOL/RE), 0.1D0)
                  dT = dT * qf
                  dT = MAX(dT, EPS)
                  
                  ! MOVE BACK TO PREVIOUS VALUE
                  p_tri = p_new
                  q_tri = q_new
                  po_tri = po_new
                  qo_tri = qo_new
                  ps_tri = ps_new
                  beta_tri = beta_new
                  Epv_tri = Epv_new
                  Epq_tri = Epq_new
            
              END IF

          END DO

          ! CHECK YIELD FUNCTION VALUE
          f = GetBSFunction(p_tri, q_tri, ps_tri, beta_tri)
          
          ! IF YIELD FUNCTION IS VIOLATED
          IF (f > FTOL) THEN

              ! CORRECT DRIFT FROM THE YIELD SURFACE AFTER EXPLICIT STEPS
              CALL MoveBackToBSSurface(p_tri, q_tri, po_tri, qo_tri, 
     1                                 ps_tri, beta_tri)            

          END IF

          ! UPDATE CURRENT VARIABLES USING TRIAL VARIABLES
          p_new = p_tri
          q_new = q_tri
          po_new = po_tri
          qo_new = qo_tri
          ps_new = ps_tri
          beta_new = beta_tri
          Epv_new = Epv_tri
          Epq_new = Epq_tri

          ! UPDATE PSEUDO TIME AND RESIZE PSEUDO TIME INCREMENT
          qs = MIN(0.9D0 * SQRT(STOL/RE), 1.1D0)

          ! IF PREVIOUS STEP IS FAILED, RESTRICT TIME INCREMENT
          IF(fail == .TRUE.) THEN
            qs = MIN(qs, 1.D0)
          END IF
          
          ! UPDATE TIME INCREMENT AND TIME
          T = T + dT
          dT = dT * qs
          dT = MAX(dT, EPS)
          dT = MIN(dT, 1.D0 - T)
          
      END DO

      END SUBROUTINE

      ! RETURN STATE INCREMENT
      SUBROUTINE GetStateIncrement(dTdEv, dTdEq,
     1                             p, q, po, qo, ps, beta,
     2                             dTdEpv, dTdEpq,
     3                             dp, dq, dpo, dqo, dps, dbeta)

      IMPLICIT NONE

      ! INPUT VARIABLES
      REAL*8 dTdEv, dTdEq
      REAL*8 p, q, po, qo, ps, beta       

      ! RETURN VARIABLES
      REAL*8 dTdEpv, dTdEpq
      REAL*8 dp, dq, dpo, dqo, dps, dbeta

      ! LOCAL VARIABLES IN THIS SUBROUTINE
      REAL*8 deno, nume, pm
         
      ! TAIL LENGTH
      REAL*8 l
        
      ! UPDATE MODEL VARIABLES
      CALL UpdateModelVariables(p, q, po, qo, ps, beta)

      ! COMPUTE PLASTIC MULTIPLIER
      nume = dfdp*K*dTdEv + 3.D0*dfdq*G*dTdEq
      deno = h + 3.D0*dfdq*G*dfdq + dfdp*K*dfdp
      pm = MAX(nume/deno, 0.D0)

      ! COMPUTE PLASTIC STRAIN INCREMENT
      dTdEpv = pm*dfdp
      dTdEpq = pm*dfdq
        
      ! COMPUTE STRESS INCREMENT
      dp = K*(dTdEv - dTdEpv)
      dq = 3.D0*G*(dTdEq - dTdEpq)
        
      ! COMPUTE INTERNAL VARIABLES' INCREMENT
      dps = pm*(1.D0 + voidratio_ini)/(lambda - kappa)*ps*dfdp
      dbeta = pm*(1.D0 + voidratio_ini)/(lambda - kappa)*q/p*dfdp

      ! COMPUTE PROJECTION CENTER'S INCREMENT
      dpo = p + dp - po
      dqo = q + dq - qo
      l = SQRT(dpo*dpo + dqo*dqo)
      dpo = MAX(l - lt, 0.D0)/l * dpo
      dqo = MAX(l - lt, 0.D0)/l * dqo
      
      END SUBROUTINE
      
      ! MOVE BACK TO BOUNDING SURFACE
      SUBROUTINE MoveBackToBSSurface(p, q, po, qo, ps, beta)
      
      IMPLICIT NONE  
            
      ! INPUT & RETURN VARIABLES
      REAL*8 p, q, po, qo, ps, beta

      ! LOCAL VARIABLES IN THIS SUBROUTINE
      REAL*8 p_tri, q_tri, po_tri, qo_tri, ps_tri, beta_tri
      !REAL*8 p_new, q_new, po_new, qo_new, ps_new, beta_new
      REAL*8 f, f1, nume, pm, dEpv, dEpq, dpo, dqo, l
      INTEGER count

      ! SET TRIALS VARIABLES
      p_tri = p
      q_tri = q
      po_tri = po
      qo_tri = qo
      ps_tri = ps
      beta_tri = beta
       
      ! SET RETURN VARIABKES
      !p_new = p
      !q_new = q
      !po_new = po
      !qo_new = qo
      !ps_new = ps
      !beta_new = beta

      ! COMPUTE CURRENT YIELD FUNCTION VALUE
      f = GetBSFunction(p_tri, q_tri, ps_tri, beta_tri)
        
      ! INITIALIZE COUNTER
      count = 0

      ! ELASTO PLASTIC RANGE
      DO WHILE(f > FTOL)

          ! UPDATE COUNTER
          count = count + 1

          ! TOO MANY ATTEMPTS ERROR
          IF(count > 30) THEN
              PRINT*,"TOO MANY ATTEMPTS"
              PAUSE
          END IF

          ! UPDATE MODEL VARIABLES
          CALL UpdateModelVariables(p_tri, q_tri, 
     1                              po_tri, qo_tri, ps_tri, beta_tri)

          ! CALCULATE PLASTIC MULTIPLIER
          nume = h + 3.D0*dfdq*G*dfdq + dfdp*K*dfdp
          pm = MAX(f / nume, 0.D0)
          
          !IF (pm .EQ. 0.D0) THEN
          !    PRINT*, 'pm:', pm
          !    PRINT*, 'nume:', nume
          !    PRINT*, 'h:', h
          !    PRINT*, 'dfdp:', dfdp
          !    PRINT*, 'dfdq:', dfdq
          !    PRINT*, 'f:', f
          !END IF
          
          ! COMPUTE PLASTIC STRAIN INCREMENT
          dEpv = pm*dfdp
          dEpq = pm*dfdq
        
          ! COMPUTE STRESS INCREMENT
          p_tri = p_tri - K*dEpv
          q_tri = q_tri - 3.D0*G*dEpq
        
          ! COMPUTE INTERNAL VARIABLES' INCREMENT
          ps_tri = ps_tri + pm*(1.D0 + voidratio_ini)
     1               /(lambda - kappa)*ps*dfdp
          beta_tri = beta_tri
     1        + pm*(1.D0 + voidratio_ini)/(lambda - kappa)
     2         *q_new/p_new*dfdp

          ! COMPUTE PROJECTION CENTER'S INCREMENT
          dpo = p_tri - po_tri
          dqo = q_tri - qo_tri
          l = SQRT(dpo*dpo + dqo*dqo)
          dpo = MAX(l - lt, 0.D0)/l * dpo
          dqo = MAX(l - lt, 0.D0)/l * dqo
          
          po_tri = po_tri + dpo
          qo_tri = qo_tri + dqo
            
          ! UPDATE RETURN VALUES
          !p_new = p_tri
          !q_new = q_tri
          !po_new = po_tri
          !qo_new = qo_tri
          !ps_new = ps_tri
          !beta_new = beta_tri
          
          f = GetBSFunction(p_tri, q_tri, ps_tri, beta_tri)
          
      END DO
        
      p = p_tri
      q = q_tri
      po = po_tri
      qo = qo_tri
      ps = ps_tri
      beta = beta_tri
        
      END SUBROUTINE

      ! UPDATE MODEL VARIABLES (IMAGE POINT, HARDENING MODULUS, 
      SUBROUTINE UpdateModelVariables(p, q, po, qo, ps, beta)

      IMPLICIT NONE

      ! INPUT POINTER
      REAL*8 p, q
      REAL*8 ps, beta
      REAL*8 po, qo
      REAL*8 check

      ! LOCAL VARIABLE
      REAL*8 dp, dq, dc, db
      REAL*8 a, b, c
        
      ! SET ELASTIC MODULI
      CALL SetElasticModuli(p)
        
      ! GET GAPS BETWEEN CURRENT STRESS AND PROJECTION CENTER 
      dp = p - po
      dq = q - qo
        
      ! CALCULATE COEFFICIENTS OF THE EQUATION FOR DISTANCE TO IMAGE POINT
      a = (dq - beta*dp)*(dq - beta*dp) + M*M*dp*dp
      b = 2.D0*((q - beta*p)*(dq - beta*dp) - M*M*(ps - p)*dp)
      c = (q - beta*p)*(q - beta*p) - M*M*(2.D0*ps-p)*p
      
      ! CALCULATE DISTANCE TO IMAGE POINT
      dc = MAX(SQRT(dp*dp + dq*dq),EPS)
      
      IF (b*b - 4.D0*a*c < 0.D0) THEN
          db = 0.D0
      ELSE
          db = (-b + SQRT(b*b - 4.D0*a*c))/2.D0/a
      END IF
      
      ! SET IMAGE POINT (GLOBAL VARIABLES)
      pimg = p + db*dp
      qimg = q + db*dq
          
      ! CALCULATE GRADIENT OF PLASTIC POTENTIAL (GLOBAL VARIABLES)
      dfdp = 2.D0*m*m*(pimg - ps) - 2.D0*beta*(qimg-pimg*beta)
      dfdq = 2.D0*(qimg - pimg*beta)
      dfdps = -2.D0*m*m*pimg
      dfdbeta = -2.D0*pimg*(qimg - pimg*beta)
        
      ! CALCULATE HARDENING MODULUS (GLOBAL VARIABLES)      
      hb = -1.D0*(1.D0 + voidratio_ini)/(lambda - kappa)*dfdp*
     1        (ps*dfdps + q/p*dfdbeta)
      
      h = hb + ho*(db/dc)**mu
      
      END SUBROUTINE

      ! CALCULATE ELASTIC G, K
      SUBROUTINE SetElasticModuli(p)

      IMPLICIT NONE
      REAL*8 p, e
      
      K = (1.D0 + voidratio_ini)*p / kappa
      G = 3.D0*(1.D0 - 2.D0*nu)*K/2.D0/(1.D0 + nu) 
      !K = 2.0D0*(1.0D0 + nu) * G/(3.D0*(1.D0-2.D0*nu))

      END SUBROUTINE

      ! BOUNDING SURFACE FUNCTION VALUE
      FUNCTION GetBSFunction(p, q, ps, beta)

      IMPLICIT NONE
      
      REAL*8 p, q, ps, beta
      REAL*8 GetBSFunction

      GetBSFunction = (q - p * beta)*(q - p * beta)
     1                - M * M * p * (2.D0 * ps- p)
      
      END FUNCTION

      ! I1 (OR TRACE)
      FUNCTION GetI1(mat)

      IMPLICIT NONE

      REAL*8 GetI1
      REAL*8 mat(1:3,1:3)

      GetI1 = mat(1,1) + mat(2,2) + mat(3,3)

      END FUNCTION

      ! I2 (OR DETERMINENT)
      FUNCTION GetI3(mat)

      IMPLICIT NONE

      REAL*8 GetI3
      REAL*8 mat(1:3,1:3)

      GetI3 = mat(1,1)*mat(2,2)*mat(3,3)
     1        - mat(1,1)*mat(2,3)*mat(3,2)
     2        - mat(1,2)*mat(2,1)*mat(3,3)
     3        + mat(1,2)*mat(2,3)*mat(3,1)
     4        + mat(1,3)*mat(2,1)*mat(3,2)
     5        - mat(1,3)*mat(2,2)*mat(3,1)

      END FUNCTION

      ! MEAN (OR TRANCE/3)
      FUNCTION GetMean(mat)

      IMPLICIT NONE

      REAL*8 GetMean      ! RETURN
      REAL*8 mat(1:3,1:3) ! PARAMETER

      GetMean = GetI1(mat) * ONE_OVER_THREE

      END FUNCTION

      ! DEVIATORIC STRESS
      FUNCTION GetDev(mat)

      IMPLICIT NONE

      REAL*8 GetDev(1:3,1:3)  ! RETURN
      REAL*8 mat(1:3,1:3)     ! PARAMETER

      GetDev = mat - GetMean(mat) * IDENTITY

      END FUNCTION

      ! NORMALIZED DEVIATORIC 
      FUNCTION GetNormDevMatrix(matrix)

      IMPLICIT NONE

      REAL*8 GetNormDevMatrix(1:3,1:3), matrix(1:3,1:3)

      GetNormDevMatrix = GetDev(matrix) / MAX(GetMean(matrix),EPS)

      END FUNCTION

      ! FROBINEOUS NORM
      FUNCTION GetFNorm(mat)

      IMPLICIT NONE

      REAL*8 GetFNorm, mat(1:3,1:3)

      GetFNorm =
     1   SQRT(mat(1,1)**2.0D0 + mat(2,2)**2.0D0+ mat(3,3)**2.0D0
     2      + mat(1,2)**2.0D0 + mat(2,1)**2.0D0
     3      + mat(1,3)**2.0D0 + mat(3,1)**2.0D0
     4      + mat(2,3)**2.0D0 + mat(3,2)**2.0D0)

      END FUNCTION

      
      FUNCTION GetJ2(mat)

      IMPLICIT NONE

      REAL*8 GetJ2            ! RETURN
      REAL*8 mat(1:3,1:3)     ! PARAMETERS
      REAL*8 dev_mat(1:3,1:3) ! LOCAL

      dev_mat = GetDev(mat)
      GetJ2 = 0.5D0 * GetFNorm(dev_mat)**2.0D0

      END FUNCTION


      FUNCTION GetMisesStress(matrix)

      IMPLICIT NONE

      REAL*8 GetMisesStress, matrix(1:3,1:3)

      GetMisesStress = SQRT(3.D0 * GetJ2(matrix))

      END FUNCTION
      
      
      FUNCTION GetShearStrain(matrix)

      IMPLICIT NONE

      REAL*8 GetShearStrain, matrix(1:3,1:3)
      
      GetShearStrain = SQRT(4.D0/3.D0 * GetJ2(matrix))

      END FUNCTION

      ! DOUBLE INNER PRODUCT
      FUNCTION GetDIProduct(mat1, mat2)

      IMPLICIT NONE

      REAL*8 GetDIProduct                     ! RETURN
      REAL*8 mat1(1:3,1:3), mat2(1:3,1:3)     ! PARAMETER
      REAL*8 tmp(1:3,1:3)                     ! LOCAL

      GetDIProduct = GetI1(MATMUL(mat1, TRANSPOSE(mat2)))

      END FUNCTION

      ! INVERSE MATRIX
      FUNCTION GetInvMat(mat)

      IMPLICIT NONE

      REAL*8 GetInvMat(1:3,1:3)   ! RETURN
      REAL*8 mat(1:3,1:3)         ! PARAMETER
      REAL*8 tmp(1:3,1:3)         ! LOCAL

      ! IF SINGULAR, PRINT ERROR
      IF(ABS(GetI3(Mat)) .LE. EPS) THEN
          PRINT*, "SINGULAR MATRIX"
          PRINT*, mat(1,1), mat(1,2), mat(1,3)
          PRINT*, mat(2,1), mat(2,2), mat(2,3)
          PRINT*, mat(3,1), mat(3,2), mat(3,3)
          PAUSE
      END IF

      tmp(1,1) = +(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))
      tmp(1,2) = -(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))
      tmp(1,3) = +(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))
      tmp(2,1) = -(mat(1,2)*mat(3,3)-mat(1,3)*mat(3,2))
      tmp(2,2) = +(mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1))
      tmp(2,3) = -(mat(1,1)*mat(3,2)-mat(1,2)*mat(3,1))
      tmp(3,1) = +(mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2))
      tmp(3,2) = -(mat(1,1)*mat(2,3)-mat(1,3)*mat(2,1))
      tmp(3,3) = +(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))

      GetInvMat = TRANSPOSE(tmp) / GetI3(mat)

      END FUNCTION

      END SUBROUTINE