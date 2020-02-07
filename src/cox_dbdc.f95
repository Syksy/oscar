        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |      DBDC - THE PROXIMAL DOUBLE BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION      | | 
        !| |                                     for                                          | |
        !| |                     COX's PROPORTIONAL HAZARDS MODEL                             | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                       by Kaisa Joki (last modified 2019)                         | |
        !| |                                                                                  | |
        !| |      Features :                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |           * Possibility to use simple stepsize determination after               | |
        !| |             each 'main iteration'.                                               | |
        !| |                                                                                  | |         
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |     The software is free for academic teaching and research purposes but I       | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    Utilizes the new version of PLQDF1 by Ladislav Luksan as a quadratic solver.      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !|                                                                                      |
        !|                                                                                      |
        !|   Codes include:                                                                     |
        !|                                                                                      |
        !|   cox_dbdc.f95       - Main call for DBDC (this file)                                |
        !|   constants.f95      - Double precision (also some parameters)                       |
        !|   bundle1.f95        - Bundle of DC component f_1                                    |
        !|   bundle2.f95        - Bundle of DC component f_2                                    |
        !|   functions.f95      - User-specified DC components f_1 and f_2 together with        |
        !|                        subgradients of DC components. Contains also user-specified   |
        !|                        initial values for parameters                                 |
        !|   dbdc.f95           - DBDC method                                                   |
        !|                                                                                      |
        !|   plqdf1.f           - Quadratic solver by Ladislav Luksan                           |
        !|                                                                                      |
        !|   Makefile           - Makefile                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software CALL SUBROUTINE in cox_dbdc.f95                                |
        !|                                                                                      |       
        !|   To ALTER the default parameters of DBDC change                                     |
        !|                  cox_dbdc.95    and    functions.f95                                 |
        !|   as needed.                                                                         |
        !|                                                                                      |
        !|                                                                                      |
        !|   References:                                                                        |
        !|                                                                                      |
        !|   [1] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. Mäkelä:               |
        !|       "A proximal bundle method for nonsmooth DC optimization utilizing              |
        !|       nonconvex cutting planes. J. Glob. Optim. 68 (2017), pp. 501-535,              | 
        !|       https://doi.org/10.1007/s10898-016-0488-3                                      | 
        !|                                                                                      |
        !|                                                                                      |
        !|   [2] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa, Marko M. Mäkelä ja Sona Taheri:   |
        !|       Double bundle method for finding Clarke stationary points in nonsmooth         |
        !|       DC programming. SIAM J. Optim., 28 (2018), s. 1892–1919.                       |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       


      MODULE cox_dbdc
      
         USE constants, ONLY : dp, ip  ! C type double precision and C type integer (i.e. accuracies)
         USE functions                 ! INFORMATION from the USER
         USE bundle1                   ! Bundle 1
         USE bundle2                   ! Bundle 2
         USE dbdc                      ! DBDC method
         
        IMPLICIT NONE   
        
        PUBLIC :: coxdbdc_loop_kits
        
        CONTAINS 
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      coxdbdc_loop_kits   : The double bundle method for Cox's model              | |                                                                                   | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
             
         
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |        THE DOUBLE BUNDLE ALGORITHM FOR COX'S PROPORTIONAL HAZARDS MODEL        |  |
        !  |                                                                                |  |
        !  |                     SOLUTION FOR EVERY NUMBER OF KITS                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE coxdbdc_loop_kits( nft, nrecord, nkits,  &
                               & in_vX, in_vY, in_vC, in_vK,  &
                               & f_for_k, beta_for_k, CPUTime, iprint, start)
			!				   
            ! * 'problem = 3'        : the used objective function is Cox's proportional hazards model with L0-norm for kits
            !                          (parametric approach since the user needs to fix the number of nonzero elements) 
            !
            !   Solves in a loop all L0-norm problems with the fixed number of k=1,...,nkits of kits
            !
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'nft'            : The dimension of the problem = the number of features in a predictor, INTEGER
            !         * 'nrecord'        : The number of records (data points), INTEGER
            !         * 'nkits'          : the (maximum) number of kits, INTEGER
            !
            !         NOTICE: DATA IS GIVEN IN VECTOR FORMAT ! Due to this values for observations/kits are given in these vector one after another
            !
            !         * 'in_vX'          : the matrix X of input variables in vector form --> will be stored to in_mX (each row is one observation)
            !                               REAL, DIMENSION(nrecord*nft)  
            !         * 'in_vY'          : the matrix Y of y_time and y_event in vector form --> will be stored to in_mY (each row cosists of time and event)
            !                               INTEGER, DIMENSION(nrecord*2) 
            !         * 'in_vC'          : the cost vector C for kits --> will be stored to in_mC
            !                               REAL, DIMENSION(knits)
            !         * 'in_vK'          : the matrix K of kit structures in vector form --> will be stored to in_mK (each row is one kit)
            !                               INTEGER, DIMENSION(knits*nft)
            !
            !         * 'iprint'         : specifies the print, INTEGER 
            !         * 'start'          : specifies how starting points are selected when the L0-norm problem is solved for the fixed number of kits 
            !                              and in which order the L0-norm problems are solved, INTEGER
            !
            ! OUTPUT: * 'f_for_k'       : The vector containing objective function values at for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(nkits) 
            !         * 'beta_for_k'    : The vector containing solution vectors beta obtained for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(nft*nkits) 
            !         * 'CPUtime'       : the CPU time (in seconds) REAL
            !
            ! NOTICE: * The dimension of vectors 'f_for_k' has to be 'nkits' ('nkits' is the maximum number of kits)
            !         * The dimension of vectors 'beta_solution' has to be 'nft*nkits' ('nft' is the dimension of the problem and 'nkits' is the number of kits)
            !
            ! OPTIONAL INPUT (CAN BE INCLUDED AS INPUT PARAMETERS IF NEEDED. AT THE MOMENT DEFAULT VALUES ARE USED FOR THEM.):
            !         * 'mit'            : The maximum number of 'main iterations', INTEGER                  
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration', INTEGER
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Clarke stationary' algorithm, INTEGER
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm, LOGICAL           
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm, LOGICAL
            !         * 'scale_in_use'   : If .TRUE. then the data is scaled
            !
            ! NOTICE: * 'iprint' has to be 0, 1, 2 or 3. If it is NOT then DEFAULT value 1 is used.             
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER *************************************           
               REAL(KIND=dp), DIMENSION(nft*nkits), INTENT(OUT) :: beta_for_k  ! the solution vectors beta obtained for the problem 3 with different k
               REAL(KIND=dp), DIMENSION(nkits), INTENT(OUT)     :: f_for_k     ! the objective function values for the problem 3 with different k
               
               REAL(KIND=dp), DIMENSION(nrecord*nft), INTENT(IN) :: in_vX    ! predictor matrix in vector format
               INTEGER(KIND=ip), DIMENSION(nrecord*2), INTENT(IN) :: in_vY            ! observed times and labels matrix in vector format
               INTEGER(KIND=ip), DIMENSION(nkits*nft), INTENT(IN) :: in_vK            ! kit matrix in vector format
               REAL(KIND=dp), DIMENSION(nkits), INTENT(IN) :: in_vC          ! kit costs in vector format              
               
               REAL(KIND=dp), INTENT(OUT) :: CPUtime                         ! the CPU time (in seconds)
                                                       
               INTEGER(KIND=ip), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=ip), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=ip), INTENT(IN) :: nkits                   ! the number of kits
               
               INTEGER(KIND=ip), INTENT(INOUT) :: iprint               ! specifies the print
                                                              !   iprint = 0 : print is suppressed
                                                              !   iprint = 1 : prints for each L0-norm problem the final value of Cox's proportional hazard model
                                                              !   iprint = 2 : prints for each L0-norm problem the final value of Cox's proportional hazard model obtained from each starting point
                                                              !   iprint = 3 : prints for each starting point and L0-norm problem all intermediate results together the final results

               INTEGER(KIND=ip), INTENT(IN) :: start                ! specifies how starting points are selected when L0-norm problem is solved for the fixed number of kits
                                                              !   start = 1  : only one starting point for L0-norm problem with k nonzero components 
                                                              !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                                              !                (starting point is the solution obtained for L0-norm problem with k-1 nonzero components)
                                                              !
                                                              !   start = -1 : only one starting point for L0-norm problem with k nonzero components 
                                                              !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1
                                                              !                (starting point is the solution obtained for L0-norm problem with k+1 nonzero components)  
                                                              !
                                                              !   start = 2  : for L0-norm problem with k nonzero components uses 'nkit-k+1' starting points
                                                              !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                                              !                (they are generated utilizing the solution of L0-norm problem with k-1 nonzero components)
                                                              !
                                                              !   start = -2 : for L0-norm problem with k nonzero components uses 'k+1' starting points
                                                              !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1 
                                                              !                (they are generated utilizing the solution of L0-norm problem with k+1 nonzero components)
                                                              !
                                                              !   start = 3  : for L0-norm problem with k nonzero components uses 'nkit' starting points 
                                                              !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                                              !                (they are generated utilizing the solution of L0-norm problem with k-1 nonzero components)
                                                              !
                                                              !   start = -3 : for L0-norm problem with k nonzero components uses 'nkit' starting points 
                                                              !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1
                                                              !                (they are generated utilizing the solution of L0-norm problem with k+1 nonzero components)
                                                              !
                                                              !   start = 4  : for L0-norm problem with k nonzero components uses 'start_max' starting points 
                                                              !                L0-problems are solved in order: k = 1, 2, ... , nkits (randomly generated)
                                                              ! 
                                                              !   start = -4 : for L0-norm problem with k nonzero components uses 'start_max' starting points 
                                                              !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1 (randomly generated)                                                              
                               
                                           
           !***************************** LOCAL VARIABLES ************************************      
              
               REAL(KIND=dp), DIMENSION(nft) :: beta_solution    ! the solution vector beta obtained for the problem
               REAL(KIND=dp) :: f_solution                       ! the objective function value at the solution 'beta_solution'

               REAL(KIND=dp), DIMENSION(nft,nkits) :: points     ! the beta_solutions for problem 3 for fixed k ('nkits' different starting points) 
               REAL(KIND=dp), DIMENSION(nkits) ::     f_points   ! the objective function values for problem 3 for fixed k ('nkits' different starting points)
               REAL(KIND=dp), DIMENSION(nft) ::       x_koe      ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=dp), DIMENSION(nft) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller
               
               REAL(KIND=dp), DIMENSION(nrecord,nft) :: in_mX    ! predictor matrix (row is an observation)
               INTEGER(KIND=ip), DIMENSION(nrecord,2) :: in_mY   ! observed times and labels matrix (row is an observation)  
               INTEGER(KIND=ip), DIMENSION(nkits,nft) :: in_mK   ! kit matrix (row is a kit)
               REAL(KIND=dp), DIMENSION(nkits) :: in_mC          ! kit costs                   

               INTEGER(KIND=ip), DIMENSION(nkits) :: kits_beta     ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=ip), DIMENSION(nkits) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'

               REAL(KIND=dp), DIMENSION(nft) :: x_0                ! the starting point
               
               REAL(KIND=dp), DIMENSION(nrecord) :: mTimes         ! Times for the observations   
               INTEGER(KIND=ip), DIMENSION(nrecord) :: mTimesInd   ! Labels of times for the observations (0=alive, 1=death)  
              
               REAL(KIND=dp), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=ip) :: nstart                 ! the number of starting point
               INTEGER(KIND=ip) :: start_max              ! the number of starting point when 'start = 5'

               INTEGER(KIND=ip) :: termination            ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=ip), DIMENSION(8) :: counter      ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=ip) :: user_n                ! the dimension of the problem
                       
               INTEGER(KIND=ip) :: mit                   ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=ip) :: mrounds               ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=ip) :: mrounds_clarke        ! the maximum number of rounds during one 'Clarke stationary' algorithm
                                                         ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
               
               INTEGER(KIND=ip) :: iprint_DBDC  ! variable that specifies print option in DBDC method:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result 
                                                !   iprint = -1: basic print of final result (without the solution vector)
                                                !   iprint = 2 : extended print of final result 
                                                !   iprint = -2: extended print of final result (without the solution vector)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results
                                                !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                                !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                                !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                                
                                                ! If 'iprint' <= -5 .OR. 'iprint' >= 5 then DEFAULT value 'iprint'=1 is used    

       
               LOGICAL :: agg_used              ! .TRUE. if aggregation is used in the algorithm. Otherwise .FALSE.                                                           
               LOGICAL :: stepsize_used         ! .TRUE. if simple stepsize determination is used in the algorithm. Otherwise .FALSE.
               
               LOGICAL :: scale_in_use          ! If .TRUE. data is scaled
               
               LOGICAL :: kit_in_use            ! If .TRUE. kit is used in the solution
               LOGICAL :: run_stop              ! If .TRUE. run is stopped for selected k
               LOGICAL :: mukana                ! If .TRUE. specific kit is in the solution
               
               LOGICAL :: ed_sol_in_pen         ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
               LOGICAL :: new_start             ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
       
               REAL(KIND=dp) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=dp) :: cost            ! The cost of solution beta
               REAL(KIND=dp) :: pieni           ! The cost of solution beta
               
               REAL(KIND=dp) :: s_time, f_time  ! The start and finish times
               REAL(KIND=dp) :: cpu             ! The cpu time
               
               REAL(KIND=dp) :: f1_current      ! The value of f1
               REAL(KIND=dp) :: f2_current      ! The value of f2
               
               REAL(KIND=dp) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=dp) :: random_num                        ! Random number
               REAL(KIND=dp), DIMENSION(nkits) :: mRand           ! Random number matrix
               INTEGER(KIND=ip), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=ip) :: problem1              ! The DC component f1 
               INTEGER(KIND=ip) :: problem2              ! The DC component f2 
        
               INTEGER(KIND=ip) :: laskuri
               INTEGER(KIND=ip) :: num_rho
               INTEGER(KIND=ip) :: nremoved
               INTEGER(KIND=ip) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=ip) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii
                
               CALL cpu_time(s_time)   ! Start CPU timing     
              
               ! The default print is used if user specifieed value is not acceptable
               IF (iprint < 0 .OR. iprint > 4) THEN
                  iprint = 1
               END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5
               
               mrounds = 5000            ! maximum number of rounds during one 'main iterations'
               mit = 5000                ! maximum number of 'main iteration'
               mrounds_clarke = 500      ! maximum number of rounds during one 'Clarke stationary' algorithm
          
               iprint_DBDC = 0           ! basic print of intermediate results and extended print of final results
          
               agg_used = .TRUE.         ! Aggregation is used
               stepsize_used = .FALSE.   ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.     ! The scaling of data is used
            
               ed_sol_in_pen = .FALSE.
               !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem

               !Values for parameter rho used in penalization problem 
               !mrho = (/0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp /) 
               mrho = (/0.2_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp /) 

                ! Problem 
                problem1 = 3
                problem2 = 3                
               
                user_n = nft
                
                tol_zero = (10.0_dp)**(-6)
       
               IF (ed_sol_in_pen) THEN 
                  WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
                              & hyodynnetaan edellisen kierroksen ratkaisua uutena lahtopisteena'
               ELSE
                  WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
                              & ei hyodynneta edellisen kierroksen ratkaisua uutena lahtopisteena'
               END IF
               
               IF (start > 0) THEN 
                  WRITE(*,*) 'Tehtavat lapi pienimmasta suurimpaan'
               ELSE   
                  WRITE(*,*) 'Tehtavat lapi suurimmasta pienimpaan'
               END IF
               
               WRITE(*,*) 'Lahtopisteet generoidaan tavalla:', start 
               WRITE(*,*) 'Penalisointitehtavan ratkaisemiseen kaytetyt rho-arvot:', mrho
               
               WRITE(*,*) 
                           
               WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
       
               
              !--------------------------------------------------------------------------------------
          
              ! The starting points for the test problems presented in MODULE functions.f95 and in the articles [1] and [2]  
        
                  x_0 = 0.0_dp
                                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------
 
               ! Populate the matrix in_mX
               ind = 0
               DO j = 1, nft
                  DO i = 1, nrecord
                     ind = ind +1
                     in_mX(i,j) = in_vX(ind)
                  END DO
               END DO
               
               ! Populate the matrix in_mY
               ind = 0
               DO j = 1, 2
                  DO i = 1, nrecord
                     ind = ind +1
                     in_mY(i,j) = in_vY(ind)
                  END DO
               END DO              
 
               ! Populate the matrix in_mK
               ind = 0
               DO j = 1, nft
                  DO i = 1, nkits
                     ind = ind +1
                     in_mK(i,j) = in_vK(ind)
                  END DO
               END DO           

               ! Populate the vector in_mC
               DO j = 1, nkits
                  in_mC(j) = in_vC(j)
               END DO
            
               WRITE(*,*) 'Matrix X populated:', in_mX
               WRITE(*,*) 'Matrix Y populated:', in_mY
               WRITE(*,*) 'Matrix K populated:', in_mK
               WRITE(*,*) 'Vector C populated:', in_mC
           
           
              !---------------------------------------------------------------------------
              !               ORDERING DATA IN INCREASING ORDER BASED ON time 
              !---------------------------------------------------------------------------
           
               DO i = 1, nrecord
                   mTimes(i) = in_mY(i,1)*1.0_dp
                   mTimesInd(i) = i
                END DO        
            
                CALL heapsort_ind(mTimes,mTimesInd)
                    
                    
               ! Allocation of data matrices in function.f95
               CALL allocate_data(nft,nrecord,nkits,user_n)   
               
              !---------------------------------------------------------------------------
              !                     STORING DATA MATRICES 
              !---------------------------------------------------------------------------
              ! Notice: Matrices are transposed! => Each column presents either an observation or a kit!
          
              DO i = 1, nrecord
                ind = mTimesInd(i)
                DO j = 1, nft
                    mX(j,i) = in_mX(ind,j)
                END DO
              END DO
                  
              DO i = 1, nrecord
                ind = mTimesInd(i)        
                DO j = 1, 2
                    mY(j,i) = in_mY(ind,j)
                END DO
              END DO    
          
              DO i = 1, nkits
                DO j = 1, nft
                    mK(j,i) = in_mK(i,j)
                END DO
              END DO          
          
              DO i = 1, nkits         
                 mC(i) = in_mC(i)
              END DO
          
              CALL failures()
              
              IF (scale_in_use) THEN
                 CALL scaling()
              END IF               
              
              ! The best beta_solution for Cox's proportional hazard model without regularization/penalization                    
                                       
              CALL DBDC_algorithm( f_solution, x_koe, x_0, 0.0_dp, 0.0_dp, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, 1, 1, user_n)            
                            
              ! Notice: * solution x_koe is obtained by fitting Cox's model to data without regularization
              !         * x_koe is utilized in formation of starting points   
              
               x_ed = 0.01_dp    ! We initialize the previous solution, since we do not have a value for it 
               WRITE(*,*) 'Value of Cox proportional hazard model without regularization:', f_solution
               WRITE(*,*) 
               
               
              !--------------------------------------------------------------------------
              !                 POPULATING AND STORING OF DATA COMPLETED 
              !---------------------------------------------------------------------------           
               
              !--------------------------------------------------------------------------
              !                 SOLVING OF L0-NORM PROBLEM BEGINS 
              !---------------------------------------------------------------------------  



              !---------------------------------------------------------
              ! problems are solved in order k=1,2,...,nkits
              !---------------------------------------------------------               
               IF (start > 0) THEN   ! problem are solved in order k=1,2,...,nkits
               
                 kit_num_ed = 0        ! The number of kits in the previous solution
                 kits_beta_ed = 0      ! The kits in the previous solution
                 
                 IF (start == 1) THEN
                    nstart = 1
                 ELSE IF (start == 2 .OR. start == 3) THEN
                    nstart = nkits
                 ELSE IF (start == 4) THEN 
                    nstart = start_max           ! The number of random starting points
                 END IF 
               
                 DO k = 1, nkits                 ! In this loop we calculate the solution for problem 3 wiht L0-norm such that the number of kits varies from 1 to k
                  CALL set_k(k)                  ! The number of nonzero kits is fixed
                  f_solution = (10.0_dp)**10     ! The best value for this solution is not yet known (therefore a big value is set)
                
                   DO i = 1, nstart              ! Different starting points are looked through to solve the problem 3 with fixed number of nonzero kits
                     x_0 = x_ed                  ! The base of the starting point is the previous solution (k-1 nonzero kits)
                     mukana = .FALSE.           
                     IF (start == 2 .OR. start == 3) THEN          ! Only nkits-k+1 starting points used             
                       DO ii = 1, kit_num_ed       ! Then we check if the kit i belongs to the previous solution
                         IF (i == kits_beta_ed(ii)) THEN   
                           mukana = .TRUE.         ! 'mukana' tells if the kit i is in solution     
                         END IF
                       END DO
                     END IF
                     
                     ! Is a new starting point generated?
                     IF (start == 2) THEN
                        IF(mukana) THEN
                           new_start = .FALSE.
                        ELSE
                           new_start = .TRUE.
                        END IF               
                     ELSE
                        new_start = .TRUE.
                     END IF                  
                                         
                     
                   IF ( new_start ) THEN       ! .TRUE. if a new starting point is generated
                    
                    IF (start == 2 .OR. start == 3) THEN 
                      DO ii = 1, nft 
                        IF (mK(ii,i)==1) THEN
                           x_0(ii) = x_koe(ii)    ! In starting poitnt x_0, we initialize features in kit i with values from solution x_koe
                        END IF   
                      END DO
                    END IF  
                    
                    IF (start == 4) THEN ! A random starting point is generated
                       x_0 = 0.01_dp
                       DO ii = 1, nkits
                         CALL RANDOM_NUMBER(random_num)
                         mRand(ii) = random_num
                         mRandInd(ii) = ii
                       END DO 
                       CALL heapsort_ind(mRand,mRandInd)
                       !WRITE(*,*) mRand
                       !WRITE(*,*) mRandInd                     
                       DO ii = 1, k
                          DO i2 = 1, nft 
                            IF (mK(i2,mRandInd(ii))==1) THEN
                              x_0(i2) = x_koe(i2)    ! In starting point x_0, we initialize features in kit i with values from solution x_koe
                            END IF   
                          END DO     
                       END DO
                    END IF                  
                      
                    run_stop = .FALSE.          ! Initialization of 'run_stop' to .FALSE. since we cannot stop
                    num_rho = 0                 ! The selection of the first value of penalization parameter rho
                  
                    IF (iprint >= 2) THEN
                       WRITE(*,*) '-------------------------------------' 
                    END IF
                  
                    DO WHILE(.NOT. run_stop)    ! The optimization begins for the selected starting point      
                  
                    num_rho = num_rho + 1       ! The update of the parameter rho
                    
                    IF (num_rho < 9) THEN
                     SELECT CASE(num_rho)        ! The selection of the parameter rho
                  
                      CASE(1)
                        rho = mrho(1)
                      CASE(2)
                        rho = mrho(2)                  
                      CASE(3)
                        rho = mrho(3)                  
                      CASE(4)
                        rho = mrho(4)              
                      CASE(5)
                        rho = mrho(5)                  
                      CASE(6)
                        rho = mrho(6)                  
                      CASE(7)
                        rho = mrho(7)                 
                      CASE(8)
                        rho = mrho(8)                             
                     END SELECT 
                    ELSE
                      rho = 10.0_dp * rho
                      
                      IF (num_rho>=24) THEN 
                      run_stop = .TRUE.         ! Forced stop in optimization
                      END IF 
                    END IF
                
                    ! The optimization problem is solved fot the current rho and x_0
     
                    CALL DBDC_algorithm( f_solution, beta_solution, x_0, rho, 0.0_dp, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n) 
            
                    IF (ed_sol_in_pen) THEN 
                      x_0 = beta_solution   ! Starting point for the next round
                    END IF    
                                
                    laskuri = 0           
                    DO j = 1, user_n         ! The number of zero components in beta vector is computed
                      IF ( ABS(beta_solution(j)) < tol_zero) THEN
                        laskuri = laskuri +1
                        beta_solution(j) = 0.0_dp
                      END IF
                    END DO
                                
                    cost = 0.0_dp            ! The initialization of the cost of 'beta_solution'
                    kit_num = 0              ! The initialization of the number of kits in 'beta_solution'
                    kits_beta = 0            ! The initialization of kits in 'beta_solution'
                  
                    !Calculation of kits in solution
                    DO j1 = 1, nkits        ! Each kit is looked through
                      kit_in_use = .FALSE.  
                      DO j2 = 1, nft
                        IF (mk(j2,j1)==1) THEN 
                          IF ( ABS(beta_solution(j2)) >= tol_zero) THEN    ! If the condition is satisfied the kit j1 is in solution 'beta_solution'
                             kit_in_use = .TRUE.                                  
                          END IF
                        END IF
                      END DO
                      IF (kit_in_use) THEN       ! Executed if kit j1 is in solution 'beta_solution'
                        cost = cost + mC(j1)     ! The cost of kit j1 is taken into account
                        kit_num = kit_num + 1    ! The number of kits is updated
                        kits_beta(kit_num) = j1  ! The index of kit j1 is updated to table kits_beta                      
                      END IF 
                    END DO
                              
                    ! We check if the optimization of problem 3 with k kits can be stopped            
                    IF (kit_num <= k .OR. run_stop) THEN 
                       IF (run_stop) THEN 
                         f_solution = (10.0_dp)**10
                       END IF
                       run_stop = .TRUE.     ! The optimization can be stopped
                       DO j = 1, nft
                         points(j,i) = beta_solution(j)        ! We store the solution
                       END DO       
                       f1_current = f1(beta_solution,3,user_n)
                       f2_current = f2(beta_solution,3,user_n)
                       f_points(i) = f1_current-f2_current      ! We store the objective funtion value without penalization
                    END IF               
                                
                    f1_current = f1(beta_solution,3,user_n)
                    f2_current = f2(beta_solution,3,user_n)
                  
                    IF (iprint == 3 .OR. run_stop) THEN                   
                      !WRITE(*,*) 'f',f1_current-f2_current, 'kits', kit_num, '<=', k 
                    END IF
                  
                    !WRITE(45,*) 'rho', rho, 'f',f1_current-f2_current, 'zero features', laskuri   ! Prints the solution to the file 'ratkaisu.txt'
                
                    IF (run_stop) THEN              
                     !WRITE(55,*) rho, f1_current-f2_current, laskuri, 38-laskuri, &
                     !& cost, kit_num, kits_beta(1:kit_num)   ! Prints the solution to the file 'ratkaisu.txt'
                    END IF
                    END DO
                   ELSE
                     !ind = (k-1)*nft       ! Solution with k nonzero features/kits is taken down
                     DO j = 1, nft
                        points(j,i) = x_ed(j)
                     END DO     
                     f_points(i) = pieni
                   END IF           
             
                   END DO
                   
                   IF (start == 1) THEN 
                      pieni = f_points(1)
                      min_ind = 1
                   END IF 
                   
                   IF (start == 2 .OR. start == 3) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     pieni = (10.0_dp)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nkits           ! The determination of the smallest objective function value
                       IF (pieni >= f_points(i)) THEN
                         pieni = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF

                   IF (start == 4) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     pieni = (10.0_dp)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nstart          ! The determination of the smallest objective function value
                       IF (pieni >= f_points(i)) THEN
                         pieni = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF                  
                               
                   ind = (k-1)*nft  
                   DO i = 1, nft
                     x_ed(i) = points(i,min_ind)            ! The beta solution yielding the smallest objective function value is set as the previous solution
                     beta_for_k(ind+i) = points(i,min_ind)  ! The beta solution yielding the smallest objective function value is stored to solution vector 
                   END DO
              
                   IF (iprint >= 1) THEN
                     !WRITE(*,*) 'Objective function value for problem 3 with', k, 'nonzero kits:' , pieni
                   END IF
              
                   laskuri = 0           
                   DO j = 1, user_n         ! The number of zero components in the previous beta vector is computed
                      IF ( ABS(x_ed(j)) < tol_zero ) THEN
                        laskuri = laskuri+1
                      END IF
                   END DO


                   cost = 0.0_dp               ! The initialization of the previous solution
                   kit_num_ed = 0              ! The initialization of the number of kits in the previous solution
                   kits_beta_ed = 0            ! The initialization of kits in the previous solution                    

                   !Calculation of kits in the previous solution
                   DO j1 = 1, nkits    ! Each kit is looked through
                     kit_in_use = .FALSE.
                     DO j2 = 1, nft
                       IF (mk(j2,j1)==1) THEN 
                         IF ( ABS(x_ed(j2)) >= tol_zero ) THEN 
                           kit_in_use = .TRUE.           
                         END IF
                       END IF
                     END DO
                     IF (kit_in_use) THEN              ! Executed if kit j1 is in the previous solution
                       cost = cost + mC(j1)            ! The cost of kit j1 is taken into account
                       kit_num_ed = kit_num_ed + 1     ! The number of kits is updated
                       kits_beta_ed(kit_num_ed) = j1   ! The index of kit j1 is updated to table kits_beta                 
                     END IF 
                   END DO       

                   WRITE(*,*) pieni, laskuri, 38-laskuri, &
                    & cost, kit_num_ed, kits_beta_ed(1:kit_num_ed)                
                    
                 END DO
              
              !---------------------------------------------------------
              ! problems are solved in order k=nkits,nkits-1,...,1
              !---------------------------------------------------------      
              ELSE     

                 kit_num_ed = 0        ! The number of kits in the previous solution
                 kits_beta_ed = 0      ! The kits in the previous solution
                 pieni = 100000.0_dp
                 
                 IF (start == -1) THEN
                    nstart = 1
                 ELSE IF (start == -2 .AND. start == -3 ) THEN
                    nstart = nkits
                 ELSE IF (start == -4) THEN     
                    nstart = start_max
                 END IF 
               
                 DO k = nkits, 1, -1            ! In this loop we calculate the solution for problem 3 wiht L0-norm such that the number of kits varies from 1 to k
                  CALL set_k(k)                 ! The number of nonzero kits is fixed
                  f_solution = 1000000.0_dp     ! The best value for this solution is not yet known (therefore a big value is set)
                
                   DO i = 1, nstart              ! Different starting points are looked through to solve the problem 3 with fixed number of nonzero kits
                     x_0 = x_ed                  ! The base of the starting point is the previous solution 
                     mukana = .FALSE.
                     IF (start == -2 .OR. start == -3 ) THEN              
                       DO ii = 1, kit_num_ed       ! Then we check if the kit i belongs to the previous solution
                         IF (i == kits_beta_ed(ii)) THEN   
                           mukana = .TRUE.         ! 'mukana' tells if the kit i is in solution  
                         END IF
                       END DO
                     END IF
                     
                     IF (start == -2) THEN
                        IF(mukana .OR. (nkits==k .AND. i==1)) THEN
                           new_start = .TRUE.
                        ELSE
                           new_start = .FALSE.
                        END IF
                     ELSE 
                        new_start = .TRUE.
                     END IF
                    
                     
                   IF (new_start) THEN    ! .TRUE. if new starting point is generated
                     
                     IF (k==nkits) THEN   ! All the kits are in the solution. Thus, the solution x_koe is the optimal solution for the problem and a good starting point
                        x_0 = x_koe
                     END IF                  
                     
                     IF (mukana) THEN 
                        DO ii = 1, nft 
                           IF (mK(ii,i)==1) THEN
                             x_0(ii) = 0.0_dp        ! In starting point x_0, we initialize features in kit i with value 0
                           END IF   
                        END DO
                     END IF
                     
                    IF (start == -3 .AND. (.NOT. mukana)) THEN
                       DO ii = 1, nft 
                         IF (mK(ii,i)==1) THEN
                           x_0(ii) = x_koe(ii)        ! In starting poitnt x_0, we initialize features in kit i with values from x_koe belonging to kit i
                         END IF   
                       END DO                    
                    END IF  
                    
                    IF (start == -4 .AND. (.NOT. k==nkits) ) THEN
                       DO ii = 1, nkits
                         CALL RANDOM_NUMBER(random_num)
                         mRand(ii) = random_num
                         mRandInd(ii) = ii
                       END DO 
                       CALL heapsort_ind(mRand,mRandInd)
                       !!!WRITE(*,*) mRandInd                      
                       IF (k > 3) THEN
                         nremoved = 0
                         ii = 1
                         DO WHILE (nremoved < 2 .AND. ii <= nkits) 
                            DO i2 = 1, kit_num_ed       ! Then we check if the kit i belongs to the previous solution
                              IF (mRandInd(ii) == kits_beta_ed(i2)) THEN 
                                 DO iii = 1, nft 
                                   IF (mK(iii,mRandInd(ii))==1) THEN
                                      x_0(iii) = 0.0_dp    ! In starting point x_0, we initialize features in kit mRandInd(ii) with value 0
                                   END IF   
                                 END DO 
                                 nremoved = nremoved + 1
                              END IF
                            END DO                          
                            ii = ii + 1
                         END DO   
                       ELSE 
                         nremoved = 0
                         ii = 1
                         DO WHILE (nremoved < 1 .AND. ii <= nkits) 
                            DO i2 = 1, kit_num_ed       ! Then we check if the kit i belongs to the previous solution
                              IF (mRandInd(ii) == kits_beta_ed(i2)) THEN 
                                 DO iii = 1, nft 
                                   IF (mK(iii,mRandInd(ii))==1) THEN
                                      x_0(iii) = 0.0_dp    ! In starting point x_0, we initialize features in kit mRandInd(ii) with value 0
                                   END IF   
                                 END DO
                                 nremoved = nremoved + 1
                              END IF
                            END DO                          
                            ii = ii + 1
                         END DO                          
                       END IF 
                    END IF
                    
                    !!!WRITE(*,*) 'start', x_0 

                    run_stop = .FALSE.          ! Initialization of 'run_stop' to .FALSE. since we cannot stop
                    num_rho = 0                 ! The selection of the first value of penalization parameter rho
                  
                    IF (iprint >= 2) THEN
                       WRITE(*,*) '-------------------------------------' 
                    END IF
                  
                    DO WHILE(.NOT. run_stop)    ! The optimization begins for the selected starting point      
                  
                    num_rho = num_rho + 1       ! The update of the parameter rho
                    
                    IF (num_rho < 9) THEN
                     SELECT CASE(num_rho)        ! The selection of the parameter rho
                  
                      CASE(1)
                        rho = mrho(1)
                      CASE(2)
                        rho = mrho(2)                  
                      CASE(3)
                        rho = mrho(3)                  
                      CASE(4)
                        rho = mrho(4)              
                      CASE(5)
                        rho = mrho(5)                  
                      CASE(6)
                        rho = mrho(6)                  
                      CASE(7)
                        rho = mrho(7)                 
                      CASE(8)
                        rho = mrho(8)                            
                     END SELECT 
                    ELSE
                      rho = 10.0_dp * rho
                      
                      IF (num_rho>=24) THEN 
                      run_stop = .TRUE.         ! Forced stop in optimization
                      END IF 
                    END IF
                
                    ! The optimization problem is solved fot the current rho and x_0
     
                    CALL DBDC_algorithm( f_solution, beta_solution, x_0, rho, 0.0_dp, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n) 
            
                    IF (ed_sol_in_pen) THEN 
                      x_0 = beta_solution   ! Starting point for the next round
                    END IF    
                                
                    laskuri = 0           
                    DO j = 1, user_n         ! The number of zero components in beta vector is computed
                      IF ( ABS(beta_solution(j)) < tol_zero) THEN
                        laskuri = laskuri +1
                        beta_solution(j) = 0.0_dp
                      END IF
                    END DO
                                
                    cost = 0.0_dp            ! The initialization of the cost of 'beta_solution'
                    kit_num = 0              ! The initialization of the number of kits in 'beta_solution'
                    kits_beta = 0            ! The initialization of kits in 'beta_solution'
                  
                    !Calculation of kits in solution
                    DO j1 = 1, nkits        ! Each kit is looked through
                      kit_in_use = .FALSE.  
                      DO j2 = 1, nft
                        IF (mk(j2,j1)==1) THEN 
                          IF ( ABS(beta_solution(j2)) >= tol_zero) THEN    ! If the condition is satisfied the kit j1 is in solution 'beta_solution'
                             kit_in_use = .TRUE.                                  
                          END IF
                        END IF
                      END DO
                      IF (kit_in_use) THEN       ! Executed if kit j1 is in solution 'beta_solution'
                        cost = cost + mC(j1)     ! The cost of kit j1 is taken into account
                        kit_num = kit_num + 1    ! The number of kits is updated
                        kits_beta(kit_num) = j1  ! The index of kit j1 is updated to table kits_beta                      
                      END IF 
                    END DO
                              
                    ! We check if the optimization of problem 3 with k kits can be stopped            
                    IF (kit_num <= k .OR. run_stop) THEN 
                       IF (run_stop) THEN 
                         f_solution = (10.0_dp)**7
                       END IF
                       run_stop = .TRUE.     ! The optimization can be stopped
                       DO j = 1, nft
                         points(j,i) = beta_solution(j)        ! We store the solution
                       END DO       
                       f1_current = f1(beta_solution,3,user_n)
                       f2_current = f2(beta_solution,3,user_n)
                       f_points(i) = f1_current-f2_current      ! We store the objective funtion value without penalization
                    END IF               
                                
                    f1_current = f1(beta_solution,3,user_n)
                    f2_current = f2(beta_solution,3,user_n)
                  
                    IF (iprint == 3 .OR. run_stop) THEN                   
                      !WRITE(*,*) 'f',f1_current-f2_current, 'kits', kit_num, '=', k 
                    END IF
                  
                    !WRITE(45,*) 'rho', rho, 'f',f1_current-f2_current, 'zero features', laskuri   ! Prints the solution to the file 'ratkaisu.txt'
                
                    IF (run_stop) THEN              
                     !WRITE(55,*) rho, f1_current-f2_current, laskuri, 38-laskuri, &
                     !& cost, kit_num, kits_beta(1:kit_num)   ! Prints the solution to the file 'ratkaisu.txt'
                    END IF
                    END DO
                    
                   ELSE
                     !ind = (k-1)*nft       ! Solution with k nonzero features/kits is taken down
                     DO j = 1, nft
                        points(j,i) = x_ed(j)
                     END DO     
                     f_points(i) = 100000000000.0_dp
                   END IF           
             
                   END DO
                   
                   IF (start == -1) THEN 
                      pieni = f_points(1)
                      min_ind = 1
                   END IF 
                   
                   IF (start == -2 .OR. start == -3) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     pieni = 10000000.0_dp     ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nkits           ! The determination of the smallest objective function value
                       IF (pieni >= f_points(i)) THEN
                         pieni = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF 
                   
                   IF (start == -4) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     pieni = 10000000.0_dp     ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nstart          ! The determination of the smallest objective function value
                       IF (pieni >= f_points(i)) THEN
                         pieni = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF                      
                               
                   ind = (k-1)*nft  
                   DO i = 1, nft
                     x_ed(i) = points(i,min_ind)            ! The beta solution yielding the smallest objective function value is set as the previous solution
                     beta_for_k(ind+i) = points(i,min_ind)  ! The beta solution yielding the smallest objective function value is stored to solution vector 
                   END DO
              
                   IF (iprint >= 1) THEN
                     ! WRITE(*,*) 'Objective function value for problem 3 with', k, 'nonzero kits:' , pieni
                   END IF
              
                   laskuri = 0           
                   DO j = 1, user_n         ! The number of zero components in the previous beta vector is computed
                      IF ( ABS(x_ed(j)) < tol_zero ) THEN
                        laskuri = laskuri+1
                      END IF
                   END DO

                   cost = 0.0_dp               ! The initialization of the previous solution
                   kit_num_ed = 0              ! The initialization of the number of kits in the previous solution
                   kits_beta_ed = 0            ! The initialization of kits in the previous solution                    

                   !Calculation of kits in the previous solution
                   DO j1 = 1, nkits    ! Each kit is looked through
                     kit_in_use = .FALSE.
                     DO j2 = 1, nft
                       IF (mk(j2,j1)==1) THEN 
                         IF ( ABS(x_ed(j2)) >= tol_zero ) THEN 
                           kit_in_use = .TRUE.           
                         END IF
                       END IF
                     END DO
                     IF (kit_in_use) THEN              ! Executed if kit j1 is in the previous solution
                       cost = cost + mC(j1)            ! The cost of kit j1 is taken into account
                       kit_num_ed = kit_num_ed + 1     ! The number of kits is updated
                       kits_beta_ed(kit_num_ed) = j1      ! The index of kit j1 is updated to table kits_beta                 
                     END IF 
                   END DO       

                   WRITE(*,*) pieni, laskuri, 38-laskuri, &
                    & cost, kit_num_ed, kits_beta_ed(1:kit_num_ed)                
                    
                 END DO

              END IF              
              
              !--------------------------------------------------------------------------
              !                 SOLVING OF L0-NORM PROBLEM ENDS 
              !---------------------------------------------------------------------------  

            
              !--------------------------------------------------------------------------
              !                 RESCALING OF SOLUTION VECTORS beta
              !---------------------------------------------------------------------------                

              IF (scale_in_use) THEN     
                 
                 CALL set_k(nkits)               
                 user_lambda = 0.0_dp
                 CALL rescaling()        ! the rescaling of data 
                 
                 DO k = 1, nkits         ! Each solution is rescaled
                   ind = (k-1)*nft
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, nft
                        beta_solution(i) = beta_for_k(ind+i)
                   END DO
                   CALL rescaling_beta(beta_solution)             ! Rescaling of solution
                   f1_current = f1(beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(beta_solution,problem1,user_n) ! The f_2 value
                   f_for_k(k) = f1_current-f2_current             ! The objective function value for problem 3 with k nonzero kits         
                   DO i = 1, nft 
                      beta_for_k(ind+i) = beta_solution(i)        ! the beta vector for problem 3 with k nonzero kits
                   END DO       
                   
                   IF (iprint == 3) THEN
                      !WRITE(*,*) 'Objective function value for problem 3 with', k, 'nonzero kits:' , f_for_k(k)
                   END IF
                 END DO 
            
              END IF            

              !--------------------------------------------------------------------------
              !                 RESCALING OF SOLUTION VECTORS beta COMPLETED
              !---------------------------------------------------------------------------  
               
             CALL cpu_time(f_time)   ! Finish CPU timing    
             cpu = f_time-s_time             
             
             IF (iprint >= 2) THEN               
                 WRITE(*,*) 'Used CPU:', cpu                 
             END IF
       
             CALL deallocate_data() 
                

         END SUBROUTINE coxdbdc_loop_kits        
  
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |        INITIALIZATION OF SEED NUMBER USED IN RANDOM NUMBER GENERATOR           |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        !This subroutine initializes the seed number used in random number generator        
           SUBROUTINE init_random_gen
              IMPLICIT NONE 
              INTEGER(KIND=ip) :: seed_size, state
              INTEGER(KIND=ip), DIMENSION(8) :: t
              INTEGER(KIND=ip), DIMENSION(:), ALLOCATABLE :: seed
 
              CALL RANDOM_SEED(size=seed_size) 
              ALLOCATE(seed(seed_size), stat=state)
              IF(state>0) STOP 'Reservation of space failed!'
              CALL DATE_AND_TIME(values=t)
              seed = 100 * t(7) + t(8)/10
              CALL RANDOM_SEED(put=seed)
           END SUBROUTINE           
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
  
      END MODULE cox_dbdc
