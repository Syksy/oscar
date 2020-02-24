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

    	use, intrinsic :: iso_c_binding
    	
        !USE constants_f, ONLY : dp, ip  ! C type double precision and C type integer (i.e. accuracies)
         
        IMPLICIT NONE   

        !REAL(KIND=c_double), PARAMETER :: dp = c_double ! Changed, was int
        !INTEGER(KIND=c_int), PARAMETER :: ip = c_int
        
        PRIVATE
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
        
        ! TDL Note: need the row splits due to too long line for subroutine ...
           SUBROUTINE coxdbdc_loop_kits(nft, nrecord, nkits, in_vX, in_vY, in_vC, in_vK, &
           	& f_for_k, beta_for_k, iprint, start) &
           	& BIND(C, name = "coxdbdc_loop_kits_f_")
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
               REAL(KIND=c_double), DIMENSION(nft*nkits), INTENT(OUT) :: beta_for_k  ! the solution vectors beta obtained for the problem 3 with different k
               REAL(KIND=c_double), DIMENSION(nkits), INTENT(OUT)     :: f_for_k     ! the objective function values for the problem 3 with different k
               
               REAL(KIND=c_double), DIMENSION(nrecord*nft), INTENT(IN) :: in_vX    ! predictor matrix in vector format
               INTEGER(KIND=c_int), DIMENSION(nrecord*2), INTENT(IN) :: in_vY            ! observed times and labels matrix in vector format
               INTEGER(KIND=c_int), DIMENSION(nkits*nft), INTENT(IN) :: in_vK            ! kit matrix in vector format
               REAL(KIND=c_double), DIMENSION(nkits), INTENT(IN) :: in_vC          ! kit costs in vector format              
               
               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nkits                   ! the number of kits
               
               INTEGER(KIND=c_int), INTENT(IN) :: iprint               ! specifies the print
                                                              !   iprint = 0 : print is suppressed
                                                              !   iprint = 1 : prints for each L0-norm problem the final value of Cox's proportional hazard model
                                                              !   iprint = 2 : prints for each L0-norm problem the final value of Cox's proportional hazard model obtained from each starting point
                                                              !   iprint = 3 : prints for each starting point and L0-norm problem all intermediate results together the final results

               INTEGER(KIND=c_int), INTENT(IN) :: start                ! specifies how starting points are selected when L0-norm problem is solved for the fixed number of kits
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
              
               REAL(KIND=c_double), DIMENSION(nft) :: beta_solution    ! the solution vector beta obtained for the problem
               REAL(KIND=c_double) :: f_solution                       ! the objective function value at the solution 'beta_solution'

               REAL(KIND=c_double), DIMENSION(nft,nkits) :: points     ! the beta_solutions for problem 3 for fixed k ('nkits' different starting points) 
               REAL(KIND=c_double), DIMENSION(nkits) ::     f_points   ! the objective function values for problem 3 for fixed k ('nkits' different starting points)
               REAL(KIND=c_double), DIMENSION(nft) ::       x_koe      ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nft) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller
               
               REAL(KIND=c_double), DIMENSION(nrecord,nft) :: in_mX    ! predictor matrix (row is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrecord,2) :: in_mY   ! observed times and labels matrix (row is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nkits,nft) :: in_mK   ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nkits) :: in_mC          ! kit costs                   
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta     ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'
               REAL(KIND=c_double), DIMENSION(nft) :: x_0                ! the starting point
               
               REAL(KIND=c_double), DIMENSION(nrecord) :: mTimes         ! Times for the observations   
               INTEGER(KIND=c_int), DIMENSION(nrecord) :: mTimesInd   ! Labels of times for the observations (0=alive, 1=death)  
              
               REAL(KIND=c_double), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=c_int) :: nstart                 ! the number of starting point
               INTEGER(KIND=c_int) :: start_max              ! the number of starting point when 'start = 5'
               INTEGER(KIND=c_int) :: termination            ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter      ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n                ! the dimension of the problem
                       
               INTEGER(KIND=c_int) :: mit                   ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds               ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds_clarke        ! the maximum number of rounds during one 'Clarke stationary' algorithm
                                                         ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
               
               INTEGER(KIND=c_int) :: iprint_DBDC  ! variable that specifies print option in DBDC method:
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
       
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
               REAL(KIND=c_double) :: pieni           ! The cost of solution beta
               
               REAL(KIND=c_double) :: s_time, f_time  ! The start and finish times
               REAL(KIND=c_double) :: cpu             ! The cpu time
               
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                        ! Random number
               REAL(KIND=c_double), DIMENSION(nkits) :: mRand           ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=c_int) :: problem1              ! The DC component f1 
               INTEGER(KIND=c_int) :: problem2              ! The DC component f2 
        
               INTEGER(KIND=c_int) :: laskuri
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: nremoved
               INTEGER(KIND=c_int) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=c_int) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii, iprint0
               
               WRITE(*,*) 'Starting subroutine in Fortran'
               
               !iprint0 = iprint
               !CALL cpu_time(s_time)   ! Start CPU timing     
              
               ! The default print is used if user specifieed value is not acceptable
               !IF (iprint0 < 0 .OR. iprint0 > 4) THEN
               !   iprint0 = 1
               !END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5
               
               mrounds = 5000            ! maximum number of rounds during one 'main iterations'
               mit = 5000                ! maximum number of 'main iteration'
               mrounds_clarke = 500      ! maximum number of rounds during one 'Clarke stationary' algorithm
          
               iprint_DBDC = 0           ! basic print of intermediate results and extended print of final results

               WRITE(*,*) 'Checkpoint 1'
          
               agg_used = .TRUE.         ! Aggregation is used
               stepsize_used = .FALSE.   ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.     ! The scaling of data is used
            
               ed_sol_in_pen = .FALSE.
               !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem
               !Values for parameter rho used in penalization problem 
               !mrho = (/0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp /) 
               !mrho = (/0.2_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp /) 

               mrho = (/0.2_c_double, 0.5_c_double, 1.0_c_double, 2.0_c_double, &
               	& 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double /) 

                ! Problem 
                problem1 = 3
                problem2 = 3                
               
                user_n = nft
                
                tol_zero = (10.0_c_double)**(-6)
       
               IF (ed_sol_in_pen) THEN 
                  WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
                              & hyodynnetaan edellisen kierroksen ratkaisua uutena lahtopisteena'
               ELSE
                  WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
                              & ei hyodynneta edellisen kierroksen ratkaisua uutena lahtopisteena'
               END IF
               
               !IF (start > 0) THEN 
               !   WRITE(*,*) 'Tehtavat lapi pienimmasta suurimpaan'
               !ELSE   
               !   WRITE(*,*) 'Tehtavat lapi suurimmasta pienimpaan'
               !END IF
               
               !WRITE(*,*) 'Lahtopisteet generoidaan tavalla:', start 
               WRITE(*,*) 'Penalisointitehtavan ratkaisemiseen kaytetyt rho-arvot:', mrho
               
               WRITE(*,*) 
                           
               WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
       

               WRITE(*,*) 'Checkpoint 2'
          
               
              !--------------------------------------------------------------------------------------
          
              ! The starting points for the test problems presented in MODULE functions.f95 and in the articles [1] and [2]  
        
                  x_0 = 0.0_c_double
                                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------
 
               !! Populate the matrix in_mX
               !ind = 0
               !DO j = 1, nft
               !   DO i = 1, nrecord
               !      ind = ind +1
               !      in_mX(i,j) = in_vX(ind)
               !   END DO
               !END DO
               ! WRITE(*,*) 'Checkpoint 3'          
               !
               !! Populate the matrix in_mY
               !ind = 0
               !DO j = 1, 2
               !   DO i = 1, nrecord
               !      ind = ind +1
               !      in_mY(i,j) = in_vY(ind)
               !   END DO
               !END DO              
               ! WRITE(*,*) 'Checkpoint 4'          
 	       !
               !! Populate the matrix in_mK
               !ind = 0
               !DO j = 1, nft
               !   DO i = 1, nkits
               !      ind = ind +1
               !      in_mK(i,j) = in_vK(ind)
               !   END DO
               !END DO           
               !! Populate the vector in_mC
               ! WRITE(*,*) 'Checkpoint 5'          
               !DO j = 1, nkits
               !   in_mC(j) = in_vC(j)
               !END DO
               ! WRITE(*,*) 'Checkpoint 6'          
            
               WRITE(*,*) 'Matrix X populated:', in_mX
               WRITE(*,*) 'Matrix Y populated:', in_mY
               WRITE(*,*) 'Matrix K populated:', in_mK
               WRITE(*,*) 'Vector C populated:', in_mC
           
              ! Varsinainen koodi tulee seuraavaksi ...
 
                
         END SUBROUTINE coxdbdc_loop_kits        


	!
	! Testing 
	!

         SUBROUTINE coxdbdc_test(nft, nrecord, nkits, in_vX, in_vY, in_vC, in_vK, f_for_k, beta_for_k, &
         	& start, iprint) &
           	& BIND(C, name = "coxdbdc_test_f_")

            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER *************************************           
               REAL(KIND=c_double), DIMENSION(nft*nkits), INTENT(OUT) :: beta_for_k  ! the solution vectors beta obtained for the problem 3 with different k
               REAL(KIND=c_double), DIMENSION(nkits), INTENT(OUT)     :: f_for_k     ! the objective function values for the problem 3 with different k
               
               REAL(KIND=c_double), DIMENSION(nrecord*nft), INTENT(IN) :: in_vX    ! predictor matrix in vector format
               INTEGER(KIND=c_int), DIMENSION(nrecord*2), INTENT(IN) :: in_vY            ! observed times and labels matrix in vector format
               INTEGER(KIND=c_int), DIMENSION(nkits*nft), INTENT(IN) :: in_vK            ! kit matrix in vector format
               REAL(KIND=c_double), DIMENSION(nkits), INTENT(IN) :: in_vC          ! kit costs in vector format              
                                                                      
               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nkits                   ! the number of kits

	       INTEGER(KIND=c_int), INTENT(IN) :: start
	       INTEGER(KIND=c_int), INTENT(IN) :: iprint

           !***************************** LOCAL VARIABLES ************************************      
              
               REAL(KIND=c_double), DIMENSION(nft) :: beta_solution    ! the solution vector beta obtained for the problem
               REAL(KIND=c_double) :: f_solution                       ! the objective function value at the solution 'beta_solution'

               REAL(KIND=c_double), DIMENSION(nft,nkits) :: points     ! the beta_solutions for problem 3 for fixed k ('nkits' different starting points) 
               REAL(KIND=c_double), DIMENSION(nkits) ::     f_points   ! the objective function values for problem 3 for fixed k ('nkits' different starting points)
               REAL(KIND=c_double), DIMENSION(nft) ::       x_koe      ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nft) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller
               
               REAL(KIND=c_double), DIMENSION(nrecord,nft) :: in_mX    ! predictor matrix (row is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrecord,2) :: in_mY   ! observed times and labels matrix (row is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nkits,nft) :: in_mK   ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nkits) :: in_mC          ! kit costs                   
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta     ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'
               REAL(KIND=c_double), DIMENSION(nft) :: x_0                ! the starting point
               
               REAL(KIND=c_double), DIMENSION(nrecord) :: mTimes         ! Times for the observations   
               INTEGER(KIND=c_int), DIMENSION(nrecord) :: mTimesInd   ! Labels of times for the observations (0=alive, 1=death)  
              
               REAL(KIND=c_double), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=c_int) :: nstart                 ! the number of starting point
               INTEGER(KIND=c_int) :: start_max              ! the number of starting point when 'start = 5'
               INTEGER(KIND=c_int) :: termination            ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter      ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n                ! the dimension of the problem
                       
               INTEGER(KIND=c_int) :: mit                   ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds               ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds_clarke        ! the maximum number of rounds during one 'Clarke stationary' algorithm
                                                         ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
               
               INTEGER(KIND=c_int) :: iprint_DBDC  ! variable that specifies print option in DBDC method:
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
       
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
               REAL(KIND=c_double) :: pieni           ! The cost of solution beta
               
               REAL(KIND=c_double) :: s_time, f_time  ! The start and finish times
               REAL(KIND=c_double) :: cpu             ! The cpu time
               
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                        ! Random number
               REAL(KIND=c_double), DIMENSION(nkits) :: mRand           ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=c_int) :: problem1              ! The DC component f1 
               INTEGER(KIND=c_int) :: problem2              ! The DC component f2 
        
               INTEGER(KIND=c_int) :: laskuri
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: nremoved
               INTEGER(KIND=c_int) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=c_int) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii, iprint0


	       ! Outputting successful testing...
	       
	       WRITE(*,*) 'Starting test subroutine'

               !iprint0 = iprint
               !CALL cpu_time(s_time)   ! Start CPU timing     
              
               ! The default print is used if user specifieed value is not acceptable
               !IF (iprint0 < 0 .OR. iprint0 > 4) THEN
               !   iprint0 = 1
               !END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5

	       WRITE(*,*) 'Checkpoint 1'
               
               mrounds = 5000            ! maximum number of rounds during one 'main iterations'
               mit = 5000                ! maximum number of 'main iteration'
               mrounds_clarke = 500      ! maximum number of rounds during one 'Clarke stationary' algorithm
          
               iprint_DBDC = 0           ! basic print of intermediate results and extended print of final results
          
               agg_used = .TRUE.         ! Aggregation is used
               stepsize_used = .FALSE.   ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.     ! The scaling of data is used

	       WRITE(*,*) 'Checkpoint 2'
               
            
               mrho = (/0.2_c_double, 0.5_c_double, 1.0_c_double, 2.0_c_double, &
               	& 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double /) 

                ! Problem 
                problem1 = 3
                problem2 = 3                
                                            
                user_n = nft
                
               ! tol_zero = (10.0_c_double)**(-6)

               ed_sol_in_pen = .FALSE.

	       WRITE(*,*) 'Checkpoint 3'
               
       
               !IF (ed_sol_in_pen) THEN 
               !   WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
               !               & hyodynnetaan edellisen kierroksen ratkaisua uutena lahtopisteena'
               !ELSE
               !   WRITE(*,*) 'Penalisointitehtavaa ratkaistaessa &
               !               & ei hyodynneta edellisen kierroksen ratkaisua uutena lahtopisteena'
               !END IF
   
	       !! Tästä eteenpäin ei käänny, ei saa jostain syystä esim. start printattua tai starting vertailua
   
   
               !IF (start > 0) THEN 
               !   WRITE(*,*) 'Tehtavat lapi pienimmasta suurimpaan'
               !ELSE   
               !   WRITE(*,*) 'Tehtavat lapi suurimmasta pienimpaan'
               !END IF
               
               !WRITE(*,*) 'Lahtopisteet generoidaan tavalla:', start 
               WRITE(*,*) 'Penalisointitehtavan ratkaisemiseen kaytetyt rho-arvot:', mrho
               
               WRITE(*,*) 
                           
               WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
       
       	       WRITE(*,*) 'Checkpoint 4'
                      
       	       !WRITE(*,*) 'in_vX'
	       !WRITE(*,*) in_vX
       	       !WRITE(*,*) 'in_vY'
	       !WRITE(*,*) in_vY
    	       !WRITE(*,*) 'in_vK'
	       !WRITE(*,*) in_vK
       	       !WRITE(*,*) 'in_vC'
	       !WRITE(*,*) in_vC

              !--------------------------------------------------------------------------------------
          
              ! The starting points for the test problems presented in MODULE functions.f95 and in the articles [1] and [2]  
        
              !    x_0 = 0.0_c_double
                                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------
 		
 		!! Seuraavat kohdat kaatavat Fortran subrutiinin
 
               ! Populate the matrix in_mX
               !ind = 0
               !DO j = 1, nft
               !   DO i = 1, nrecord
               !      ind = ind +1
               !      in_mX(i,j) = in_vX(ind)
               !   END DO
               !END DO
               
	       WRITE(*,*) 'Checkpoint 5'
               
               ! Populate the matrix in_mY
               !ind = 0
               !DO j = 1, 2
               !   DO i = 1, nrecord
               !      ind = ind +1
               !      in_mY(i,j) = in_vY(ind)
               !   END DO
               !END DO              
	       WRITE(*,*) 'Checkpoint 6'
               
               ! Populate the matrix in_mK
               !ind = 0
               !DO j = 1, nft
               !   DO i = 1, nkits
               !      ind = ind +1
               !      in_mK(i,j) = in_vK(ind)
               !   END DO
               !END DO           
               ! Populate the vector in_mC
	       WRITE(*,*) 'Checkpoint 7'
               
               !DO j = 1, nkits
               !   in_mC(j) = in_vC(j)
               !END DO
            
               !WRITE(*,*) 'Matrix X populated:', in_mX
               !WRITE(*,*) 'Matrix Y populated:', in_mY
               !WRITE(*,*) 'Matrix K populated:', in_mK
               !WRITE(*,*) 'Vector C populated:', in_mC
           
              ! Varsinainen koodi tulee seuraavaksi ...


         END SUBROUTINE coxdbdc_test        
  
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
              INTEGER(KIND=c_int) :: seed_size, state
              INTEGER(KIND=c_int), DIMENSION(8) :: t
              INTEGER(KIND=c_int), DIMENSION(:), ALLOCATABLE :: seed
 
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
      