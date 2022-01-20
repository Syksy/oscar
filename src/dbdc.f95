        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |      DBDC - THE PROXIMAL DOUBLE BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION      | | 
        !| |                                                                                  | |
        !| |                                for OSCAR with                                    | |
        !| |                                                                                  | |
        !| |                       1) Cox's proportional hazard model                         | |
        !| |                       2) Mean square error model                                 | |
        !| |                       3) Logistic regression model                               | |
        !| |                                                                                  | |
        !| |                       by Kaisa Joki (last modified  May 2021)                    | |
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


      MODULE bundle1  

      
        USE, INTRINSIC :: iso_c_binding

        
        IMPLICIT NONE 
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |           THE BUNDLE ELEMENT AND THE BUNDLE OF THE DC COMPONENT F_1              | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   

        TYPE bundle_element1 ! bundle element of F_1
           PRIVATE
           REAL (KIND=c_double), DIMENSION(:), POINTER  :: subgrad    ! subgradient of the bundle element
           REAL (KIND=c_double) :: lin_error                          ! linearization error of the bundle element
        END TYPE bundle_element1        

        TYPE kimppu1 ! bundle of F_1
           PRIVATE
           TYPE(bundle_element1), DIMENSION(:), POINTER :: b_elements   ! bundle elements (does NOT contain the 'current element' and 'agg_element')
           TYPE(bundle_element1) :: current_element ! bundle element calculated at the current iteration point ('current element') 
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           TYPE(bundle_element1) :: agg_element ! the aggregated bundle element ('agg_element')
           
           INTEGER(KIND=c_int) :: n         ! number of variables (also the length of subgradients)
           INTEGER(KIND=c_int) :: b_maxsize ! 'maximum size of the bundle' - 1, (i.e. b_maxsize=size(b_elements) NOTICE: the 'current element' and 'agg_element' are stored separately)        
           INTEGER(KIND=c_int) :: b_size    ! the current size of the bundle without the 'current element' and 'agg_element' (the actual size of the bundle is 'b_size+1' and 'agg_element is NOT taken into account in this value)
           INTEGER(KIND=c_int) :: indeksi   ! the place where the next bundle element is tried to be added in the bundle element table 'b_elements'  
          
           LOGICAL :: full      ! tells whether the bundle is full or not             
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           LOGICAL :: agg       ! tells whether the aggregated bundle element was inserted into the bundle during the previous round           
        END TYPE kimppu1 


        CONTAINS
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                           CONTAINS SUBROUTINES:                                  | | 
        !| |                                                                                  | |
        !| |    INITIALIZATION         : init_bundle_b1(set, set_size, grad_length)           | |
        !| |    ADD ELEMENT            : add_element_b1(set, grad, alpha)                     | |
        !| |    ADD AGGREGATED ELEMENT : add_agg_element_b1(set, grad, alpha)                 | |
        !| |    ADD 1. CURRENT ELEMENT : add_first_element_b1(set, grad)                      | |
        !| |    UPDATE BUNDLE          : update_b1(set, new_grad, d, value_change)            | |
        !| |    RESETS BUNDLE          : reset_b1(set)                                        | |
        !| |    DEALLOCATION           : deallocate_b1(set)                                   | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |              CONTAINS FUNCTIONS GIVING DIFFERENT VALUES:                         | |   
        !| |                                                                                  | |
        !| |    MATRIX OF SUBGRADIENTS      : grad_matrix(set)                                | |
        !| |    MATRIX OF LIN. ERRORS       : lin_error_matrix(set)                           | |
        !| |    MATRIX OF SUBGRADIENTS +AGG : grad_matrix_agg(set)                            | |
        !| |    MATRIX OF LIN. ERRORS  +AGG : lin_error_matrix_agg(set)                       | |       
        !| |    BUNDLE SIZE                 : give_size_b1(set)                               | |
        !| |    NUMBER OF VARIABLES         : give_n_b1(set)                                  | |
        !| |    IS BUNDLE FULL?             : is_full_b1(set)                                 | |
        !| |    IS AGGREGATION USED?        : is_agg_used(set)                                | |
        !| |    SUBGRADIENT OF ELEMENT i    : give_subgrad_b1(set, i)                         | | 
        !| |    LIN. ERROR OF ELEMENT i     : give_linerr_b1(set, i)                          | | 
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        
        
        

        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                SUBROUTINES                                       | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        
        !**************************************************************************************
        !                                                                                     |
        !                               INITIALIZATION                                        |
        !                                                                                     |
        !************************************************************************************** 
           
           SUBROUTINE init_bundle_b1(set, set_size, grad_length) 
               !
               ! Initializes the bundle 'set'. Now the size of the bundle is 'set_size' and the length of subgradients is 'grad_size'.
               ! 
               ! 
               ! NOTICE: * 'grad_length' >= 1
               !         * IF (set_size < 2 ) THEN the size of the bundle is set to be 1 and only the 'current element' is stored (If aggregation is used, then also the aggregated element 'agg_element' is stored)               
               !         * 'set_size' does NOT include the 'aggregated element'. So if aggregation is used, then the actual size of the bundle is 'set_size+1'. 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set          ! bundle
               INTEGER(KIND=c_int), INTENT(IN):: set_size, grad_length  ! bundle size and the length of subgardients
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, allocstat
               
                           
               IF (set_size < 2) THEN
                    set%b_maxsize = 0                       ! only the 'current element' and the 'agg_element' (if used) are stored 
                    set%full = .TRUE.
                
               ELSE     
                    set%b_maxsize = set_size - 1            ! the biggest possible size of the bundle without the 'current element' (the 'agg_element' is not taken into account here)
                    set%indeksi = 1
                    set%full = .FALSE.
               END IF
               
               set%b_size = 0                    ! the number of stored bundle elements in the table 'b_elements' ( ! without the 'current element' and 'agg_element' ! )  
               set%n = grad_length               ! the number of variables (this is also the length of subgradients)
               set%agg = .FALSE.
               
               ALLOCATE(set%b_elements(set%b_maxsize), STAT=allocstat)  ! initializes the maximum size of the bundle table 'b_elements' 
               ALLOCATE(set%current_element%subgrad(grad_length), &     ! initializes the length of the subgradient in the 'current element'
                         & STAT=allocstat)             
               ALLOCATE(set%agg_element%subgrad(grad_length), &         ! initializes the length of the subgradient in the 'aggregated element'
                         & STAT=allocstat)  

               DO i=1, set%b_maxsize
                    ALLOCATE(set%b_elements(i)%subgrad(grad_length), &  ! initializes the length of subgradients in the table 'b_elements'
                        & STAT=allocstat)     
               END DO 
               
           END SUBROUTINE init_bundle_b1
           
           

        !**************************************************************************************
        !                                                                                     |
        !                     ADD ELEMENT INTO TO THE BUNDLE                                  | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_element_b1(set, grad, alpha)
               !
               ! Adds the element '(grad, alpha)' into the bundle 'set' (i.e. into the bundle element table 'b_elements').
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.
               !         * IF the size of the bundle is 1, THEN nothing is added to the bundle element table 'b_elements'.         
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: grad   ! the subgradient                     
               REAL(KIND=c_double), INTENT(IN) :: alpha                   ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i
               
               IF (set%b_maxsize > 0 ) THEN ! executed if bundle is larger than 1 (i.e. something can be stored into the table 'b_elements')
                   IF ( set%indeksi > set%b_maxsize ) THEN
                       set%indeksi = 1         
                   END IF
               
                   i = set%indeksi
                   set%b_elements(i)%subgrad = grad     ! adds the new subgradient into position i
                   set%b_elements(i)%lin_error = alpha  ! adds the new linearization error into position i
                   set%indeksi = i + 1                  ! the position where the next element is tried to be added
                   
                   IF ( .NOT. set%full ) THEN           ! if the bundle was not full during the previous round, then the size of the bundle is increased with 1
                       set%b_size = set%b_size + 1
                   END IF                  
                   
                   IF(set%b_size == set%b_maxsize) THEN  ! we test: Is the bundle full ?
                       set%full = .TRUE.
                   ELSE
                       set%full = .FALSE.
                   END IF 
                   
               END IF
               
           END SUBROUTINE add_element_b1
                 
           
           
        !**************************************************************************************
        !                                                                                     |
        !                     ADD AGGREGATED ELEMENT INTO TO THE BUNDLE                       | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_agg_element_b1(set, grad, alpha)
               !
               ! Adds the aggregated element '(grad, alpha)' into the bundle 'set'.
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.           
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: grad   ! the aggregated subgradient                      
               REAL(KIND=c_double), INTENT(IN) :: alpha                   ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************
                           
               set%agg_element%subgrad = grad
               set%agg_element%lin_error = alpha 
               set%agg = .TRUE.
               
           END SUBROUTINE add_agg_element_b1           
       
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                  INITIALIZE/ADD THE FIRST CURRENT ELEMENT                           |
        !                                                                                     |     
        !**************************************************************************************            
           
           SUBROUTINE add_first_element_b1(set, grad)
               !
               ! Adds the element '(grad, 0)' calculated at the first iteration point x_0 into the bundle 'set'.
               !
               ! NOTICE: * the dimension of the 'grad' has to be 'set%n'.
               !         * the linearization error of the first current element is always zero.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                   ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: grad    ! the subgradient    
               !**************************** OTHER VARIABLES **************************************            
               
               set%current_element%subgrad = grad
               set%current_element%lin_error = 0.0_c_double ! the linearization error is zero at the iteration point x_0
               
           END SUBROUTINE add_first_element_b1  

           
           
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE update_b1(set, new_grad, d, value_change)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this updates also all the other linearization errors in the bundle 'set' 
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * f1(x_{k+1}) - f1(x_k) is the 'value_change'             
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: new_grad   ! the subgradient calculated at the new iteration point               
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=c_double), INTENT(IN) :: value_change                ! f1(x_{k+1}) - f1(x_k), i.e. the value change in the objective function
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i            
               
               ! The old 'current element' is added into the bundle set 'b_elements' and after that 
               ! the 'current element' can be updated with the new element
               ! (the linearization error of the current element is always zero, thus it is not changed).                  
               CALL add_element_b1(set, set%current_element%subgrad, 0.0_c_double)
               set%current_element%subgrad = new_grad
           
               ! Linearization error update in the bundle set 'b_elements'
                
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error + &
                          & value_change - DOT_PRODUCT(set%b_elements(i)%subgrad, d)
               END DO
               
               IF (set%agg) THEN
                   set%agg_element%lin_error = set%agg_element%lin_error + &
                          & value_change - DOT_PRODUCT(set%agg_element%subgrad, d)           ! update of aggregated element  
               END IF
               
           
           END SUBROUTINE update_b1 
           
           
        !**************************************************************************************
        !                                                                                     |
        !                               RESET THE BUNDLE                                      |
        !                                                                                     |
        !************************************************************************************** 
        
           SUBROUTINE reset_b1(set)
               !
               ! deletes all the elements from the bundle 'set' except the 'current element' (Also the 'aggregated element' is deleted)
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set ! bundle

               IF (set%b_maxsize > 0) THEN      ! Reset is executed if it is possible that we have something in the bundle element table 'b_elements'
                   set%b_size = 0               ! the number of stored subgradient in 'b_elements' after reset
                   set%indeksi = 1              ! the place where the next element is placed
                   set%full = .FALSE.           ! the bundle is not full because elements from 'b_elements' are removed
               END IF

               set%agg = .FALSE.            ! the aggregated element is deleted
               
           END SUBROUTINE reset_b1             
          
        !**************************************************************************************
        !                                                                                     |     
        !                                DEALLOCATION                                         |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE deallocation_b1(set)
               !
               ! Deallocates arrays     
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, allocstat             
               
               DEALLOCATE(set%current_element%subgrad, STAT=allocstat)     ! deallocates the subgradient in the 'current element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"              
               DEALLOCATE(set%agg_element%subgrad, STAT=allocstat)         ! deallocates the subgradient in the 'aggregated element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"                          

               DO i=1, set%b_maxsize
                    DEALLOCATE(set%b_elements(i)%subgrad, STAT=allocstat)  ! deallocates subgradients in the table 'b_elements'
                    IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"             
                    
               END DO            
               DEALLOCATE(set%b_elements, STAT=allocstat)  ! deallocates the bundle table 'b_elements' 
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"              

               
           END SUBROUTINE deallocation_b1         
           
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                        FUNCTIONS GIVING DIFFERENT VALUES                         | |   
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   


        !**************************************************************************************
        !                                                                                     |
        !                             MATRIX OF SUBGRADIENTS                                  |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+1)'.
               !         * The subgradient corresponding to the current element is in the last position.
               !         * The aggregated element is NOT taken into account.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=c_double), DIMENSION(set%n*(set%b_size+1)) :: m  ! the subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i, j, length, start
               
               length = set%n               ! the length of subgradients
                       
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! the place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO
           END FUNCTION grad_matrix
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                           VECTOR OF LINEARIZATION ERRORS                            |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+1'.
               !         * The linearization error corresponding to the current element is in the last position.    
               !         * The aggregated element is NOT taken into account.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=c_double), DIMENSION(set%b_size+1) :: m  ! the linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the last position

           END FUNCTION lin_error_matrix

        !**************************************************************************************
        !                                                                                     |
        !                         MATRIX OF SUBGRADIENTS WITH AGGREGATION                     |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix_agg(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' (with the aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+2)'.
               !         * The subgradient corresponding to the current element is in the 'last but one' position.
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=c_double), DIMENSION(set%n*(set%b_size+2)) :: m  ! subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i, j, length, start
               
               length = set%n               ! the length of subgradients
               
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO

               start = (set%b_size + 1) * length  ! the place where the 'aggregated element' is copied             
               DO j = 1, length 
                   m(start + j) = set%agg_element%subgrad(j)
               END DO
               
           END FUNCTION grad_matrix_agg
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                  VECTOR OF LINEARIZATION ERRORS WITH AGGREGATION                    |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix_agg(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' (with aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction. 
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+2'.
               !         * The linearization error corresponding to the current element is in the 'last but one' position.  
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=c_double), DIMENSION(set%b_size+2) :: m  ! linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the 'last but one' position
               m(set%b_size+2) = set%agg_element%lin_error      ! the lin. error corresponding to the 'aggregated element' is copied into the last position

           END FUNCTION lin_error_matrix_agg           
        


           
        !**************************************************************************************
        !                                                                                     |
        !                                    BUNDLE SIZE                                      |
        !                                                                                     |
        !**************************************************************************************
       
           PURE FUNCTION give_size_b1(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element and the aggregated element ! )
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: bundle_size           ! size of the bundle
               
               bundle_size = set%b_size
               
           END FUNCTION give_size_b1
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              NUMBER OF VARIABLES                                    |
        !                                                                                     |
        !**************************************************************************************        
               
           PURE FUNCTION give_n_b1(set) RESULT(variable_number)
               !
               ! Gives the number of varibles in the minimization problem (this is also the length of subgradients)
               ! Number of varibles is (i.e. has to be) same as in kimppu2 when used in the algorithm.
               !
               IMPLICIT NONE 
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               INTEGER(KIND=c_int) :: variable_number        ! the number of variables
               
               variable_number = set%n
               
           END FUNCTION give_n_b1          
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS BUNDLE FULL?                                     |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_full_b1(set) RESULT(isfull)
               !
               ! Returns .TRUE. if bundle 'set' (i.e. the bundle element table 'b_element') is full otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isfull                 ! tells whether the bundle is full or not
 
               isfull = set%full 
               
            END FUNCTION is_full_b1         

           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS AGGREGATION USED?                                |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_agg_used(set) RESULT(isUsed)
               !
               ! Returns .TRUE. if aggregation is used in the bundle 'set' otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isUsed                 ! tells whether the aggregation is used or not
 
               isUsed = set%agg
               
            END FUNCTION is_agg_used            
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                            SUBGRADIENT OF ELEMENT i                                 |
        !                                                                                     |
        !**************************************************************************************
        
           FUNCTION give_subgrad_b1(set, i) RESULT(grad)
               !
               ! Gives the subgradient of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the subgradient of the 'current element'.    
               !         * If 'i=-1' then gives the subgradient of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)

               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set         ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i                 ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double), DIMENSION(set%n) :: grad  ! subgradient at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! subgradient is from the bundle element table 'b_elements'
                   grad = set%b_elements(i)%subgrad
               ELSE IF (i == 0) THEN                        ! subgradient is from the 'current element'
                   grad = set%current_element%subgrad
               ELSE IF (i == -1) THEN                       ! subgradient is from the 'aggregated element'
                   grad = set%agg_element%subgrad                  
               ELSE                                         ! otherwise we are outside the bundle     
                   WRITE(*,*) 'CANNOT RETURN SUBGRADIENT! index ' &
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_subgrad_b1

           
           
        !**************************************************************************************
        !                                                                                     |
        !                          LINEARIZATION ERROR OF ELEMENT i                           |
        !                                                                                     |
        !**************************************************************************************        
    
           FUNCTION give_linerr_b1(set, i) RESULT(error)
               !
               ! Gives the linearization error of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the linearization error of the 'current element'.    
               !         * If 'i=-1' then gives the linearization error of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i         ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double) :: error           ! the linearization error at the position 'i'
            
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! the linearization error is from the bundle element table 'b_elements'
                   error = set%b_elements(i)%lin_error
               ELSE IF (i==0) THEN                          ! the linearization error is from the 'current element'
                   error = set%current_element%lin_error
               ELSE IF (i==-1) THEN                         ! the linearization error is from the 'aggregated element'      
                   error = set%agg_element%lin_error               
               ELSE                                         ! otherwise we are outside the bundle
                   WRITE(*,*) 'CANNOT RETURN LINEARIZATION ERROR! index '&
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_linerr_b1         

      END MODULE bundle1


      MODULE bundle2
      
        USE, INTRINSIC :: iso_c_binding
      

        IMPLICIT NONE 

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |           THE BUNDLE ELEMENT AND THE BUNDLE OF THE DC COMPONENT F_2              | |  
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        TYPE bundle_element2 ! bundle element of F_2
           PRIVATE
           REAL (KIND=c_double), DIMENSION(:), POINTER  :: subgrad   ! subgradient
           REAL (KIND=c_double), DIMENSION(:), POINTER  :: direction ! search direction for the subproblem           
           REAL (KIND=c_double) :: lin_error     ! linearization error
           REAL (KIND=c_double) :: change        ! value of the predicted decrease  
           REAL (KIND=c_double) :: subprob_value ! value of the subproblem objective
        END TYPE bundle_element2        

        TYPE kimppu2 ! bundle of F_2
           PRIVATE
           TYPE(bundle_element2), DIMENSION(:), POINTER :: b_elements ! bundle elements
           TYPE(bundle_element2) :: current_element ! bundle element at the current iteration point ('current element')
           
           INTEGER(KIND=c_int) :: n         ! number of variables in vector x (also the length of subgradients)
           INTEGER(KIND=c_int) :: b_maxsize ! 'maximum size of the bundle' - 1 ,(i.e. b_maxsize=size(b_elements) NOTICE: the current_element is stored separately)   
           INTEGER(KIND=c_int) :: b_size    ! the current size of the bundle without the 'current element' (the actual size of the bundle is 'b_size+1')
           INTEGER(KIND=c_int) :: glob_ind  ! the position of the bundle element giving the global solution    
           INTEGER(KIND=c_int) :: indeksi   ! the place where the next element is tried to be added in the bundle element table 'b_elements'       
           
           LOGICAL :: full     ! tells whether this bundle is full or not          
        END TYPE kimppu2


        CONTAINS
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAIN SUBROUTINES:                                  | | 
        !| |                                                                                  | |
        !| |    INITIALIZATION               : init_bundle_b2(set, set_size, grad_length)     | |
        !| |    ADD ELEMENT                  : add_element_b2(set, grad, alpha)               | |
        !| |    FIRST CURRENT ELEMENT        : add_first_element_b2(set, grad)                | |
        !| |    UPDATE BUNDLE                : update_b2(set, new_grad, d, value_change)      | |
        !| |    SOLUTION FOR SUBPROBLEM i    : add_solution(set, i , d, delta, obj)           | |   
        !| |    ADD INDEX OF GLOBAL SOLUTION : add_glob_index(set)                            | |
        !| |    RESET BUNDLE:                : reset_b2(set)                                  | |       
        !| |    DEALLOCATION:                : deallocation_b2(set)                           | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |              CONTAIN FUNCTIONS GIVING DIFFERENT VALUES:                          | |   
        !| |                                                                                  | |
        !| |    GLOBAL SOLUTION                   : give_solution(set)                        | |
        !| |    SOLUTION OF SUBPROBLEM i          : give_subprob_solution(set,i)              | |
        !| |    PREDICTED DECREASE                : give_decrease(set)                        | |
        !| |    PREDICT. DECRAESE OF SUBPROBLEM i : give_subprob_decrease(set,i)              | |
        !| |    INDEX OF SOLUTION                 : give_solution_ind(set)                    | |
        !| |    INDEX OF LAST ELEMENT             : give_last_element_ind_b2(set)             | |
        !| |    BUNDLE SIZE                       : give_size_b2(set)                         | |
        !| |    MAX BUNDLE SIZE                   : give_max_size_b2(set)                     | |
        !| |    NUMBER OF VARIABLES               : give_n_b2(set)                            | |
        !| |    IS BUNDLE FULL?                   : is_full_b2(set)                           | |
        !| |    SUBGRADIENT OF ELEMENT i          : give_subgrad_b2(set, i)                   | | 
        !| |    LIN. ERROR OF ELEMENT i           : give_linerr_b2(set, i)                    | | 
        !| |    MAXIMUM NORM VALUE                : max_norm_value(set)                       | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        
        
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                SUBROUTINES                                       | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*               
        

        !**************************************************************************************
        !                                                                                     |
        !                               INITIALIZATION                                        |
        !                                                                                     |
        !**************************************************************************************
           
           SUBROUTINE init_bundle_b2(set, set_size, grad_length) 
               !
               ! Initializes the bundle 'set'. Now the size of the bundle is 'set_size' and the length of subgradients is 'grad_length'.
               !
               ! NOTICE: * 'grad_length' >= 1
               !         * IF (set_size < 2 ) THEN the size of the bundle is set to be 1 and only the 'current element' is stored 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set          ! bundle
               INTEGER(KIND=c_int), INTENT(IN):: set_size, grad_length  ! the size of the bundle and the length of subgradients
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i, allocstat
            
               IF (set_size < 2) THEN
                    set%b_maxsize = 0               ! only the 'current element' is stored 
                    set%full = .TRUE. 
               ELSE     
                    set%b_maxsize = set_size - 1    ! the biggest possible size of the bundle without the current element
                    set%indeksi = 1
                    set%full = .FALSE.
               END IF

               set%b_size = 0           ! the number of stored bundle elements in the table 'b_elements' ( ! without the current element ! )  
               set%n = grad_length      ! the number of variables (this is also the length of subgradients)
            
               ALLOCATE(set%b_elements(set%b_maxsize), STAT=allocstat)   ! initializes the maximum size of the bundle element table 'b_elements'

               DO i=1, set%b_maxsize
                    ALLOCATE(set%b_elements(i)%subgrad(grad_length), &   ! initializes the length of subgradients in the table 'b_elements'
                        & STAT=allocstat) 
                    ALLOCATE(set%b_elements(i)%direction(grad_length), & ! initializes the length of searh directions in the table 'b_elements'
                        & STAT=allocstat)       
               END DO  
               
               ALLOCATE(set%current_element%subgrad(grad_length), &      ! initialize the length of the subgradient in the 'current element'
                        & STAT=allocstat) 
               ALLOCATE(set%current_element%direction(grad_length), &    ! initialize the length of the searh direction in the 'current element'
                        & STAT=allocstat) 
                        
           END SUBROUTINE init_bundle_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                     ADD ELEMENT INTO TO THE BUNDLE                                  | 
        !                                                                                     |
        !**************************************************************************************
        
           SUBROUTINE add_element_b2(set, grad, alpha)
               !
               ! Adds the new element (grad, alpha) into the bundle 'set' (i.e. into the bundle element table 'b_elements').
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.
               !         * IF the size of the whole bundle is 1, THEN nothing is added to the bundle element table 'b_elements'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: grad ! the subgradient                       
               REAL(KIND=c_double), INTENT(IN) :: alpha                 ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i
               
               IF (set%b_maxsize > 0 ) THEN ! executed if bundle size is larger than 0 (i.e. something can be stored into the table 'b_elements')
                   IF ( set%indeksi > set%b_maxsize ) THEN
                        set%indeksi = 1
                   END IF

                   i = set%indeksi  
                   
                   ! In the algorithm we use the case where the 'bundle element' yielding the previous global solution of the search direction problem cannot be replaced
                   
                   IF( set%full .AND. (i == set%glob_ind ) ) THEN   ! the bundle element which yields the previous global solution of 
                       i = i + 1                                    ! the search direction problem cannot be replaced
                       IF ( i > set%b_maxsize ) THEN                ! we make sure that the updated index is not outside the bundle element table 'b_elements'
                            i = 1
                       END IF 
                   END IF              
                   
                   set%b_elements(i)%lin_error = alpha   ! adds a new linearization error into the position i                      
                   set%b_elements(i)%subgrad = grad      ! adds a new subgradient into the position i

                   set%indeksi = i + 1                   ! the position where the next element is tried to be added
                   
                   IF ( .NOT. set%full ) THEN            ! if the bundle wasn't full during the previous round then the size of the bundle is increased with 1
                      set%b_size = set%b_size + 1
                   END IF
                               
                   IF(set%b_size == set%b_maxsize) THEN  ! we test: Is the bundle full ?
                       set%full = .TRUE.
                   ELSE
                       set%full = .FALSE.
                   END IF 
                   
               END IF
           END SUBROUTINE add_element_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                  INITIALIZE/ADD THE FIRST CURRENT ELEMENT                           |
        !                                                                                     |     
        !**************************************************************************************        
           
           SUBROUTINE add_first_element_b2(set, grad)
               !
               ! Adds the element '(grad, 0)' calculated at the first iteration point x_0 into the bundle 'set'.
               !
               ! NOTICE: * the dimension of the 'grad' has to be 'set%n'.
               !         * the linearization error of the first current element is always zero.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: grad ! subgradient at the first iteration point x_0                         
               
               set%current_element%subgrad = grad
               set%current_element%lin_error = 0.0_c_double ! linearization error is zero at the iteration point x_0      
           
           END SUBROUTINE add_first_element_b2      
           
           
  
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************
        
           SUBROUTINE update_b2(set, new_grad, d, value_change)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this also the linearization errors are updated in the bundle
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * f2(x_{k+1}) - f2(x_k) is the 'value_change'
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN):: new_grad   ! subgradient calculated at the new iteration point                   
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=c_double), INTENT(IN) :: value_change                ! f2(x_{k+1}) - f2(x_k), i.e. the value change in the objective function 
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: i            
               
               ! the old 'current element' is added to the bundle set 'b_elements' and after that 
               ! the 'current_element' can be updated with the new element 
               ! (the linearization error of the current element is always zero, thus it is not changed)           
               CALL add_element_b2(set, set%current_element%subgrad, 0.0_c_double)
               set%current_element%subgrad = new_grad
               
           
               !linearization error update on the bundle set 'b_elements'
           
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error + &
                          & value_change - DOT_PRODUCT(set%b_elements(i)%subgrad, d)
               END DO
           
           END SUBROUTINE update_b2
           
           
        !**************************************************************************************
        !                                                                                     |
        !            ADD DIRECTION, PREDICTED DECREASE AND SUBPROBLEM OBJECTIVE VALUE         |
        !                                FOR THE SUBPROBLEM i                                 |
        !                                                                                     |
        !************************************************************************************** 
        
            SUBROUTINE add_solution(set, i , d, delta, obj )
               !
               ! Adds the search direction 'd', the predicted decrease 'delta' (delta1+delta2)
               ! and the objective value 'obj' related to the subproblem 'i'
               !
               ! NOTICE: * 0 <= i <= set%b_size  (other indices are outside the current bundle)
               !         * If i=0 values are added to the current element
               !         * the dimension of 'd' has to be 'set%n'
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set              ! bundle
               INTEGER(KIND=c_int) :: i                                     ! the index of the subproblem
               REAL(KIND=c_double), DIMENSION(set%n), INTENT(IN) :: d ! the search direction
               REAL(KIND=c_double), INTENT(IN) :: delta, obj          ! the predicted decraese and the objective value of the subproblem
               
               IF ( (i > set%b_size) .OR. (i < 0) ) THEN        ! subproblem index is outside (i.e. does not belong to) the bundle  
                   WRITE(*,*) 'CANNOT ADD SOLUTION! index ', i ,&
                         & 'outside the bundle (size:', set%b_size, ')'
               ELSE IF (i == 0) THEN                            ! values are added to the current element
                   set%current_element%direction = d
                   set%current_element%change = delta
                   set%current_element%subprob_value = obj                 
               ELSE                                             ! values are added to the bundle element table 'b_elements' (position is 'i')
                   set%b_elements(i)%direction = d
                   set%b_elements(i)%change = delta
                   set%b_elements(i)%subprob_value = obj
               END IF             
            END SUBROUTINE add_solution
            
            
           
        !**************************************************************************************
        !                                                                                     |
        !               ADD INDEX OF SUBPROBLEM YIELDING THE GLOBAL SOLUTION                  |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE add_glob_index(set) 
               !
               ! Calculates the index of the subproblem which 
               ! gives the global solution of the 'search direction' problem
               !
               ! NOTICE: * 'glob_index' is from the interval from '0' to 'b_size'.
               !         * IF 'glob_index' = 0, THEN the 'current element' gives the solution.
               !         * OTHERWISE the global solution is from the bundle element table 'b_elements'
               !           and 'glob_index' is the position of 'b_elements' which yields the global solution.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set   ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: ind, i                    
               
               IF (set%b_size == 0 ) THEN   ! true if only the 'current element' is in the bundle (i.e. we have only one element in the bundle).
                    set%glob_ind = 0        ! In this case the subproblem related to the 'current element' yields the global solution.
               ELSE            
                    ind = 1 
                    DO i = 2, set%b_size    ! calculation of the index of the element yielding the minimum value of the objective in the table 'b_elements'
                        IF ( set%b_elements(ind)%subprob_value &
                             & >  set%b_elements(i)%subprob_value ) THEN
                              ind = i
                        END IF
                    END DO

                    IF( set%b_elements(ind)%subprob_value &                ! the minimum value of the objective from the table 'b_elements' is compared
                             & > set%current_element%subprob_value ) THEN  ! with the 'current element'
                        ind = 0
                    END IF  
                    
                    set%glob_ind = ind
               END IF                   
           END SUBROUTINE add_glob_index
           

        !**************************************************************************************
        !                                                                                     |
        !                               RESET THE BUNDLE                                      |
        !                                                                                     |
        !************************************************************************************** 
        
           SUBROUTINE reset_b2(set)
               !
               ! Deletes from the bundle 'set' all the elements except the 'current element'
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set ! bundle

               IF (set%b_maxsize > 0) THEN   ! Reset is executed if it is possible that we have something in the bundle element table 'b_elements'
                   set%b_size = 0            ! the number of stored subgradient in 'b_elements' after reset
                   set%indeksi = 1           ! the place where the next element is placed
                   set%full = .FALSE.        ! the bundle is not full because elements from 'b_elements' are removed
               END IF
                       
           END SUBROUTINE reset_b2         
                     

        !**************************************************************************************
        !                                                                                     |     
        !                                  DEALLOCATION                                       |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE deallocation_b2(set)
               !
               ! deallocates arrays        
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(INOUT) :: set                      ! bundle
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, allocstat             
               
               DO i=1, set%b_maxsize
                    DEALLOCATE(set%b_elements(i)%subgrad, STAT=allocstat)    ! deallocates subgradients in the table 'b_elements'
                    IF (allocstat /= 0) STOP "*** Could Not release memory B2 ***"             
                    DEALLOCATE(set%b_elements(i)%direction, STAT=allocstat)  ! deallocates searh directions in the table 'b_elements'
                    IF (allocstat /= 0) STOP "*** Could Not release memory B2 ***"                          
               END DO  

               DEALLOCATE(set%b_elements , STAT=allocstat)   ! deallocates the bundle element table 'b_elements'    
               IF (allocstat /= 0) STOP "*** Could Not release memory B2 ***"                      
 
               DEALLOCATE(set%current_element%subgrad, STAT=allocstat)      ! deallocates the subgradient in the 'current element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B2 ***"              
               
               DEALLOCATE(set%current_element%direction, STAT=allocstat)    ! deallocates the searh direction in the 'current element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B2 ***"                          
               
           END SUBROUTINE deallocation_b2   

           

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                        FUNCTIONS GIVING DIFFERENT VALUES                         | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*                 

        
        !**************************************************************************************
        !                                                                                     |
        !                            GLOBAL SOLUTION VECTOR                                   |
        !                                                                                     |
        !**************************************************************************************     
        
           PURE FUNCTION give_solution(set) RESULT(solution)
               !
               ! Gives the search direction 'solution' which is the global solution.
               !
               ! NOTICE: * 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !           since otherwise the index of global solution is not right.
               !         * the dimension of 'solution' is 'set%n'
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set             ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double), DIMENSION(set%n) :: solution  ! global solution vector
               
               IF (set%glob_ind > 0) THEN    ! global solution is from the bundle element table 'b_elements'
                    solution = set%b_elements(set%glob_ind)%direction
               ELSE                          ! global solution is from the 'current_element'
                    solution = set%current_element%direction
               END IF
           END FUNCTION give_solution
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                           SOLUTION VECTOR OF SUBPROBLEM i                           |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_subprob_solution(set,i) RESULT(solution)
               !
               ! Gives the search direction 'solution' of subproblem 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= set%b_size (other indices are outside the current bundle B_2).
               !         * IF i=0 THEN gives the solution of the subproblem which is get by using the 'current element' of B_2.
               !         * The dimension of 'solution' is 'set%n'.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i         ! index of the subproblem
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double), DIMENSION(set%n) :: solution ! solution vector of the subproblem 'i'
               
               IF ( (i > 0) .AND. (i <= set%b_size) ) THEN   ! solution vector is from the bundle element table 'b_elements'  
                    solution = set%b_elements(i)%direction
               ELSE IF ( i == 0) THEN                        ! solution vector is from the 'current element'
                    solution = set%current_element%direction
               END IF
           END FUNCTION give_subprob_solution   
           


        !**************************************************************************************
        !                                                                                     |
        !                              PREDICTED DECREASE                                     |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_decrease(set) RESULT(dec)
               !
               ! Gives the value of predicted decrease (i.e. delta_1 + delta_2) at the global solution.
               !
               ! NOTICE: 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !         since otherwise the index of global solution is not right.            
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double) :: dec             ! value of the predicted decrease at the global solution 

               IF (set%glob_ind > 0) THEN                       ! predicted decrease is from the bundle element table 'b_elements'
                    dec = set%b_elements(set%glob_ind)%change   
               ELSE                                             ! predicted decrease is from the 'current element'
                    dec = set%current_element%change          
               END IF
           END FUNCTION give_decrease       
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                         PREDICTED DECREASE OF SUBPROBLEM i                          |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_subprob_decrease(set,i) RESULT(dec)
               !
               ! Gives the value of predicted decrease 'dec' in the subproblem 'i' (i.e. delta_1 + delta_2 in the subproblem 'i').
               ! 
               ! NOTICE: * 0 <= 'i' <= set%b_size (other indices are outside the current bundle B_2).
               !         * IF i=0 THEN gives the predicted decrease of the subproblem which is get by using the 'current element' of B_2.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i         ! index of the subproblem   
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double) :: dec             ! value of the predicted decrease of subproblem 'i'
               
               IF ( (i > 0) .AND. (i <= set%b_size) ) THEN  ! index is from the bundle element table 'b_elements'
                    dec = set%b_elements(i)%change 
               ELSE IF (i==0) THEN                          ! index is from the 'current element'
                    dec = set%current_element%change
               END IF
           END FUNCTION give_subprob_decrease              
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              INDEX OF SOLUTION                                      |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION give_solution_ind(set) RESULT(ind)
               !
               ! Gives the index 'ind' of the subproblem which yields the global solution.
               !
               ! NOTICE: * 'CALL add_glob_index(set)' has to be executed before using this FUNCTION, 
               !           since otherwise the index of global solution is not right.                       
               !         * 'ind' is zero if the current element gives the solution.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: ind                   ! index of the global solution 
               
               ind = set%glob_ind
           END FUNCTION give_solution_ind
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                             INDEX OF THE LAST ELEMENT                               |
        !                                                                                     |
        !**************************************************************************************            
         
           PURE FUNCTION give_last_element_ind_b2(set) RESULT(ind)
               !
               ! Gives the index 'ind' of the place where the last bundle element was added in the bundle element table 'b_elements'.
               !
               ! NOTICE: The index 'ind' is zero if there is nothing in the bundle element table 'b_elements'.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: ind                   ! index of the place where the last element was added
               
               IF (set%b_size /= 0) THEN        ! there is something in the bundle element table 'b_elements'
                   ind = set%indeksi - 1 
               ELSE                             ! there in nothing in the bundle element table 'b_elements'
                   ind = 0
               END IF 
           END FUNCTION give_last_element_ind_b2    



        !**************************************************************************************
        !                                                                                     |
        !                                BUNDLE SIZE                                          |
        !                                                                                     |
        !************************************************************************************** 
       
           PURE FUNCTION give_size_b2(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element ! ).
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: bundle_size            ! size of the bundle
               
               bundle_size = set%b_size
           END FUNCTION give_size_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                               MAX BUNDLE SIZE                                       |
        !                                                                                     |
        !************************************************************************************** 
       
           PURE FUNCTION give_max_size_b2(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element ! ).
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: bundle_size            ! size of the bundle
               
               bundle_size = set%b_maxsize
           END FUNCTION give_max_size_b2           
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              NUMBER OF VARIABLES                                    |
        !                                                                                     |
        !**************************************************************************************            
       
           PURE FUNCTION give_n_b2(set) RESULT(variable_number)
               !
               ! Gives the number of varibles in the minimization problem (this is also the length of subgradients).
               ! Number of variables is (i.e. has to be) same as in kimppu1 when used in algorithm.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER(KIND=c_int) :: variable_number        ! number of variables
               
               variable_number = set%n
           END FUNCTION give_n_b2          
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              IS BUNDLE FULL?                                        |
        !                                                                                     |
        !**************************************************************************************            

           PURE FUNCTION is_full_b2(set) RESULT(isfull)
               !
               ! Returns .TRUE. if bundle is full otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************            
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isfull                ! tells whether the bundle is full or not
 
               isfull = set%full 
            END FUNCTION is_full_b2     
            
            
            
        !**************************************************************************************
        !                                                                                     |
        !                           SUBGRADIENT OF SUBPROBLEM i                               |
        !                                                                                     |
        !**************************************************************************************             
            
           FUNCTION give_subgrad_b2(set, i) RESULT(grad)
               ! 
               ! Gives the subgradient of the bundle element at position 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= 'set%b_size' (otherwise we are outside the bundle).
               !         * If 'i'=0 then gives the subgradient of the 'current element'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set         ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i                 ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double), DIMENSION(set%n) :: grad  ! subgradient at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! subgradient is from the bundle element table 'b_elements'
                   grad = set%b_elements(i)%subgrad
               ELSE IF (i == 0) THEN                        ! subgradient is from the 'current element'
                   grad = set%current_element%subgrad
               ELSE                                         ! otherwise we are outside the bundle         
                   WRITE(*,*) 'CANNOT RETURN SUBGRADIENT! index ' &
                        & , i , 'outside the bundle (size without current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_subgrad_b2
                   
           
           
        !**************************************************************************************
        !                                                                                     |
        !                       LINEARIZATION ERROR OF SUBPROBLEM i                           |
        !                                                                                     |
        !**************************************************************************************            
 
           FUNCTION give_linerr_b2(set, i) RESULT(error)
               !
               ! Gives the linearization error of the bundle element at position 'i'.
               !
               ! NOTICE: * 0 <= 'i' <= 'set%b_size' (otherwise we are outside the bundle).
               !         * If 'i'=0 then gives the linearization error of the 'current element'.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set ! bundle
               INTEGER(KIND=c_int), INTENT(IN) :: i         ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double) :: error           ! the linearization error at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! the linearization error is from the bundle element table 'b_elements'
                   error = set%b_elements(i)%lin_error
               ELSE IF (i==0) THEN                          ! the linearization error is from the 'current element'
                   error = set%current_element%lin_error
               ELSE                                         ! otherwise we are outside the bundle 
                   WRITE(*,*) 'CANNOT RETURN LINEARIZATION ERROR! index '&
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_linerr_b2
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              MAXIMUM NORM VALUE                                     |
        !                                                                                     |
        !**************************************************************************************            
        
           PURE FUNCTION max_norm_value(set) RESULT(max_norm)
               !
               ! Gives the value of the maximum subgradient norm.
               !
               IMPLICIT NONE    
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu2), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=c_double) :: max_norm         ! value of the maximum subgradient norm
               REAL(KIND=c_double) :: norm             ! 'help variable'
               INTEGER(KIND=c_int) :: i
               
               max_norm = DOT_PRODUCT(set%current_element%subgrad, &   ! the square of the norm ||\bxi_2(x_k)|| 
                                         & set%current_element%subgrad) 
                                     
               DO i = 1, set%b_size
                   norm = DOT_PRODUCT(set%b_elements(i)%subgrad, &     ! the square of the norm ||\bxi_2(y_i)||
                                  & set%b_elements(i)%subgrad)
                   IF (max_norm < norm) THEN  
                        max_norm = norm
                   END IF
               END DO
               
               max_norm = SQRT(max_norm)  ! the maximum norm ||\bxi_{2,max}|| (WITHOUT square ! )
           END FUNCTION max_norm_value
           

      END MODULE bundle2


      MODULE functions    


        USE, INTRINSIC :: iso_c_binding
        
        IMPLICIT NONE
        
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                 INFORMATION SUPPLIED BY THE USER:                        | | 
        !| |                                                                          | |
        !| |    * PROBLEM specification:                                              | |
        !| |                                                                          | |       
        !| |        - The DC componen f1:           'problem1'                        | |
        !| |        - The DC componen f2:           'problem2'                        | |
        !| |        - the number of variables:      'user_n'                          | |               
        !| |                                                                          | |
        !| |    * Different PARAMETERS:                                               | |
        !| |                                                                          | |       
        !| |        - the stopping tolerance:      'user_crit_tol'                    | |
        !| |                                                                          | |       
        !| |        MAIN ITERATION:                                                   | |       
        !| |        - the size of bundle B_1:      'user_size_b1'                     | |       
        !| |        - the size of bundle B_2:      'user_size_b2'                     | |       
        !| |        - the descent parameter:       'user_m'                           | |       
        !| |        - the decrease parameter:      'user_c'                           | |       
        !| |        - the decrease parameter:      'user_r_dec'                       | |       
        !| |        - the increase parameter:      'user_r_inc'                       | |           
        !| |        - the enlargement parameter:   'user_eps_1'                       | |       
        !| |                                                                          | |       
        !| |        CLARKE STATIONARY ALGORITHM:                                      | |       
        !| |        - the size of bundle B:        'user_size'                        | | 
        !| |        - the proximity measure:       'user_eps'                         | |               
        !| |        - the descent parameter:       'user_m_clarke'                    | |       
        !| |                                                                          | |       
        !| |                                                                          | |       
        !| |    * Computation of the value of the DC functions f_1 and f_2:           | |
        !| |        - f1(y, problem1)   the value of DC component f_1 at a point y    | |
        !| |        - f2(y, problem2)   the value of DC component f_2 at a point y    | |           
        !| |                                                                          | |               
        !| |                                                                          | |               
        !| |    * Computation of the subgradient of the DC components f_1 and f_2:    | |
        !| |        - subgradient_f1(y, problem1)    the subgradient of f_1 at y      | |
        !| |        - subgradient_f2(y, problem2)    the subgradient of f_2 at y      | |       
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  INFORMATION ABOUT PARAMETERS:                            | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------       
        
    
        !****************** GLOBAL PARAMETERS *******************************************

        ! Cox's proportional hazard model + L1-norm:  
        ! Test problem 1  : problem1 = 1   &  problem2 = 1

        ! Cox's proportional hazard model + L0-norm for features (parametric):  
        ! Test problem 2  : problem1 = 2   &  problem2 = 2      
        
        ! Cox's proportional hazard model + L0-norm for kits (parametric):  
        ! Test problem 3  : problem1 = 3   &  problem2 = 3    
        
        ! Mean square error + L0-norm for kits (parametric):  
        ! Test problem 4  : problem1 = 4   &  problem2 = 4  
        
        ! Logistic regression + L0-norm for kits (parametric):  
        ! Test problem 5  : problem1 = 5   &  problem2 = 5          
        
        
        !-------------------------------------------------------------------------------------------------      
        !__________________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
        !****************** PARAMETRES NEEDED ONLY IN MAIN ITERATION AlGORITHM ********************



        TYPE set_info ! all information in the bundle
            
          INTEGER(KIND=c_int) :: user_size_b1                     ! The biggest possible size of the bundle B_1 
                                                                      ! If user_size_b1 <= 0 then DEFAULT value MIN(user_n+5,1000) is used 
                                                                   
          INTEGER(KIND=c_int) :: user_size_b2                     ! The biggest possible size of the bundle B_2 
                                                                      ! If user_size_b2 <= 0 then DEFAULT value 3 is used  
              
          REAL(KIND=c_double) :: user_m                           ! The descent parameter  If user_m <= 0.0_c_double .OR. user_m >= 1.0_c_double 
                                                                               !               then DEFAULT value 0.2_c_double is used

          REAL(KIND=c_double) :: user_c                           ! The decrease parameter c in DBDC  
                                                                      ! If user_c <= 0.0_c_double or user_c > 1.0_c_double then DEFAULT value 0.1_c_double is used 

                                                                 
          REAL(KIND=c_double) :: user_r_dec                       ! The decrease parameter r in DBDC 
        
        !If user_r_dec <= 0.0_c_double .OR. user_r_dec >= 1.0_c_double then DEFAULT value is used.
        !                               
        !   DEFAULT value:                          
        !     If user_n < 10:           user_r_dec = 0.75_c_double    
        !     If 10 <= user_n < 300:    user_r_dec = the first two decimals of n/(n+5)
        !     If user_n >= 300:         user_r_dec = 0.99_c_double
        !
        !   Some examples of the DEFAULT value of the parameter 'user_r_dec':
        !     If user_n=10:     user_r_dec = 0.66_c_double                          
        !     If user_n=20:     user_r_dec = 0.80_c_double                         
        !     If user_n=25:     user_r_dec = 0.83_c_double                         
        !     If user_n=50:     user_r_dec = 0.90_c_double                         
        !     If user_n=100:    user_r_dec = 0.95_c_double                         
        !     If user_n=150:    user_r_dec = 0.96_c_double     
        !     If user_n=200:    user_r_dec = 0.97_c_double                      
        !     If user_n=250:    user_r_dec = 0.98_c_double    
        !
        
          REAL(KIND=c_double) :: user_r_inc                       ! The increase parameter R: If user_r_inc <= 1.0_c_double 
                                                                      !                           then DEFAULT value (10.0_c_double)**7 is used
 
          REAL(KIND=c_double) :: user_eps_1                       ! The enlargement parameter: If user_eps_1 <= 0.0_c_double .OR. user_eps_1 > 1.0_c_double 
                                                                      !                            then DEFAULT value 5*(10.0_c_double)**(-5) is used

 
        !____________________________________________________________________________________________                                                       
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETRES NEEDED ONLY IN CLARKE STATIONARY AlGORITHM *******************
        
        
           INTEGER(KIND=c_int) :: user_size                    ! The biggest possible bundle size for the bundle used in 'Clarke stationary' algorithm

           REAL(KIND=c_double) :: user_m_clarke                ! The descent parameter: If user_m_clarke <= 0.0_c_double .OR. user_m_clarke >= 1.0_c_double 
                                                                   !                                    then DEFAULT value 0.01_c_double is used
                                                                    
          
           REAL(KIND=c_double) :: user_eps                                 ! The proximity measure: If user_eps <= 0.0_c_double 
                                                                               !                        then DEFAULT value (10.0_c_double)**(-6) is used when n <= 50
                                                                               !                                           (10.0_c_double)**(-5) is used when n > 50
                                                                    
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETER NEEDED IN STOPPING CONDITIONS ********************
              
        ! ** The stopping tolerance **      
          REAL(KIND=c_double) :: user_crit_tol                                  ! If user_crit_tol <= 0.0_c_double then DEFAULT value (10.0_c_double)**(-5) is used when n <=200
                                                                                    !                                                     (10.0_c_double)**(-4) is used when n > 20     
       
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        
        
        !-------------------------------------------------------------------------------
        !                             DATA MATRICES
        !-------------------------------------------------------------------------------
         
          INTEGER(KIND=c_int) :: nft0              ! the maximum number of features in a predictor
          INTEGER(KIND=c_int) :: nrecord0          ! the maximum number of observations (data points)
          INTEGER(KIND=c_int) :: nkits0            ! the maximum number of kits
          INTEGER(KIND=c_int) :: nk0               ! the k-norm used

          REAL(KIND=c_double), DIMENSION(:,:), ALLOCATABLE :: mX          ! predictor matrix (column is an observation)
          REAL(KIND=c_double), DIMENSION(:), ALLOCATABLE :: mC            ! kit costs

          INTEGER(KIND=c_int), DIMENSION(:,:), ALLOCATABLE :: mY          ! observed times and labels matrix (column is an observation) in Cox's proportional hazard model  
          INTEGER(KIND=c_int), DIMENSION(:), ALLOCATABLE :: mY_log        ! vector of outputs in the logistic regression model
          
          REAL(KIND=c_double), DIMENSION(:), ALLOCATABLE :: mY_mse        ! vector of outputs in the mean square error model  
        
          INTEGER(KIND=c_int), DIMENSION(:,:), ALLOCATABLE :: mK          ! kit matrix (column is a kit)
        
          INTEGER(KIND=c_int), DIMENSION(:), ALLOCATABLE :: k_norm_ind    ! The last k places have the indices defining the k-norm
          INTEGER(KIND=c_int), DIMENSION(:), ALLOCATABLE :: k_norm_ind_k  ! The last k places have the indices defining the k-norm in kit structure
          INTEGER(KIND=c_int), DIMENSION(:,:), ALLOCATABLE :: mFail       ! The matrix containing information about failures (time and index)
          INTEGER(KIND=c_int), DIMENSION(:,:), ALLOCATABLE :: mUnique     ! The matrix containing first index for each failure time and the number of failures at that time
         
          REAL(KIND=c_double), DIMENSION(:,:), ALLOCATABLE :: hav         ! The matrix containing means and standard deviations for features
        
          REAL(KIND=c_double) :: user_rho                              ! the penalization parameter in L0-norm     
          REAL(KIND=c_double) :: user_lambda                           ! the penalization parameter in L1-norm     
          REAL(KIND=c_double) :: user_a                                ! the penalization parameter        
          INTEGER(KIND=c_int) :: nfail                                 ! The number of failures
          INTEGER(KIND=c_int) :: nfailunique                           ! The number of unique failures

        !REAL(KIND=c_double), DIMENSION(nft,nrecord), SAVE :: mX       ! predictor matrix (column is an observation)
        !INTEGER(KIND=c_int), DIMENSION(2,nrecord), SAVE :: mY         ! observed times and labels matrix (column is an observation)  
        !REAL(KIND=c_double), DIMENSION(nrecord), SAVE :: mY_mse       ! outputs in mean square error model  
        !INTEGER(KIND=c_int), DIMENSION(nrecord), SAVE :: mY_log       ! outputs in logistic regression model  
        !INTEGER(KIND=c_int), DIMENSION(nft,nkits), SAVE :: mK         ! kit matrix (column is a kit)
        !REAL(KIND=c_double), DIMENSION(nkits), SAVE :: mC             ! kit costs
         
        !INTEGER(KIND=c_int), DIMENSION(nft), SAVE :: k_norm_ind       ! The first k places have the indices defining the k-norm
        
          END TYPE set_info                 
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                          | |
        !| |                                                                                          | |
        !| |    SUBROUTINES:                                                                          | |
        !| |                                                                                          | |
        !| |       allocate_parameters((in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc,  &              | |
        !| |                           & in_eps1, in_b, in_m_clarke, in_eps, in_crit_tol)             | |
        !| |                                                                                          | |
        !| |       allocate_data_cox(set, nft,nrecord,nkits,nk)  - allocates size of data matrices    | |
        !| |       allocate_matrices_cox(set, ...)               - allocates data matrices            | |
        !| |       deallocate_data_cox(set)                      - deallocates data matrices          | |
        !| |                                                                                          | |
        !| |       allocate_data_mse(set,nft,nrecord,nkits,nk)   - allocates size of data matrices    | |
        !| |       allocate_matrices_mse(set, ...)               - allocates data matrices            | |
        !| |       deallocate_data_mse(set)                      - deallocates data matrices          | |
        !| |                                                                                          | |
        !| |       allocate_data_log(set,nft,nrecord,nkits,nk)   - allocates size of data matrices    | |
        !| |       allocate_matrices_mse(set, ...)               - allocates data matrices            | |
        !| |       deallocate_data_log(set)                      - deallocates data matrices          | |
        !| |                                                                                          | |
        !| |       set_k(nk)                               - defines the L0-norm used                 | |
        !| |       scaling_cox(set)                        - scaling of data                          | |
        !| |       scaling_matrix_cox(set, in_mX)          - scaled data matrix is stored to in_mX    | |
        !| |       rescaling_cox(set)                      - rescaling of data                        | |
        !| |       rescaling_beta_cox(y)                   - rescaling of solution                    | |
        !| |       scaling_mse(set)                        - scaling of data                          | |
        !| |       scaling_matrix_mse(set, in_mX)          - scaled data matrix is stored to in_mX    | |
        !| |       rescaling_mse(set)                      - rescaling of data                        | |
        !| |       rescaling_beta_mse(y)                   - rescaling of solution                    | |
        !| |       scaling_log(set)                        - scaling of data                          | |
        !| |       scaling_matrix_log(set, in_mX)          - scaled data matrix is stored to in_mX    | |
        !| |       rescaling_log(set)                      - rescaling of data                        | |
        !| |       rescaling_beta_log(y)                   - rescaling of solution                    | |
        !| |       failures(set)                           - Allocates the failure matrix and         | |
        !| |                                                 initializes its values                   | |       
        !| |                                                                                          | |
        !| |       f1_sub(y,problem1,user_n,f1,grad1)  - value and subgradient of f1 at y             | |
        !| |                                                                                          | |
        !| |       f1(y,problem1,user_n)               - value of f1 at y                             | |
        !| |       f2(y,problem1,user_n)               - value of f2 at y                             | |
        !| |                                                                                          | |
        !| |       subgradient_f1(y,problem1,user_n)   - subgradient of f1 at y                       | |
        !| |       subgradient_f2(y,problem2,user_n)   - subgradient of f2 at y                       | |
        !| |                                                                                          | |
        !| |       SUBROUTINES to use headsort algorithm                                              | |
        !| |       heapsort(a)                                                                        | |
        !| |       heapsort_ind(a,b)                                                                  | |
        !| |       heapsort_k(a,b,k)                                                                  | |
        !| |                                                                                          | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
    
        
        CONTAINS
        
        
           !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_parameters(set, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)
               !
               ! Allocates the parameters use in DBDC method
               !       
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information 
               
               INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance 
               
               !**************************** OTHER VARIABLES **************************************
               
               set%user_size_b1 = in_b1
               set%user_size_b2 = in_b2
               set%user_size = in_b
               set%user_m = in_m
               set%user_m_clarke = in_m_clarke
               set%user_c = in_c
               set%user_r_dec = in_r_dec                   
               set%user_r_inc = in_r_inc    
               set%user_eps_1 = in_eps1
               set%user_eps = in_eps
               set%user_crit_tol = in_crit_tol
               
               set%user_a = 100.0_c_double
           
           END SUBROUTINE allocate_parameters
           !------------------------------------------------------------------------------------------      
        
           !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_data_cox(set, nft, nrecord, nkits, nk)
               !
               ! Allocates the sizes of data matrices for the Cox's proportiona hazard model.
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information 
               
               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nkits                   ! the number of kits
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! defines the k-norm used
               !**************************** OTHER VARIABLES **************************************
               
               set%nft0 = nft
               set%nrecord0 = nrecord
               set%nkits0 = nkits
               set%nk0 = nk
               
               ALLOCATE(set%mX(nft,nrecord),set%mY(2,nrecord),set%mK(nft,nkits), & 
                        & set%mC(nkits),set%k_norm_ind(nft), &
                        & set%k_norm_ind_k(nkits),set%hav(2,nft))            
           
           END SUBROUTINE allocate_data_cox
           !------------------------------------------------------------------------------------------


           !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_data_mse(set, nft, nrecord, nkits, nk)
               !
               ! Allocates the sizes of data matrices for the mean square error model.
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information  
               
               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nkits                   ! the number of kits
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! defines the k-norm used
               !**************************** OTHER VARIABLES **************************************
               
               set%nft0 = nft
               set%nrecord0 = nrecord
               set%nkits0 = nkits
               set%nk0 = nk
               
               ALLOCATE(set%mX(nft,nrecord),set%mY_mse(nrecord),set%mK(nft,nkits), & 
                        & set%mC(nkits),set%k_norm_ind(nft), &
                        & set%k_norm_ind_k(nkits),set%hav(2,nft+1))            
           
           END SUBROUTINE allocate_data_mse
           !------------------------------------------------------------------------------------------           
          

           !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_data_log(set, nft, nrecord, nkits, nk)
               !
               ! Allocates the sizes of data matrices for the logistic regression model.
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information               
               
               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nkits                   ! the number of kits
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! defines the k-norm used
               !**************************** OTHER VARIABLES **************************************
               
               set%nft0 = nft
               set%nrecord0 = nrecord
               set%nkits0 = nkits
               set%nk0 = nk
               
               ALLOCATE(set%mX(nft,nrecord),set%mY_log(nrecord),set%mK(nft,nkits), & 
                        & set%mC(nkits),set%k_norm_ind(nft), &
                        & set%k_norm_ind_k(nkits),set%hav(2,nft+1))            
                   
           END SUBROUTINE allocate_data_log
           !------------------------------------------------------------------------------------------  


 
            !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_matrices_cox(set, in_mX, in_mY, in_mK, in_mC,  & 
                                            & nrecord, nft, nk)
               !
               ! Allocates the data matrices in Cox's proportiona hazard model
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set                       ! The set of information        

               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(2,nrecord) :: in_mY         ! observed times and labels matrix (column is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          

               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! the number of kits
                
               !**************************** OTHER VARIABLES **************************************           
               INTEGER :: i, j               

              
                  DO j = 1, nrecord
                    DO i = 1, nft
                        set%mX(i,j) = in_mX(i,j)
                    END DO
                  END DO
                      
                  DO j = 1, nrecord
                    DO i = 1, 2
                        set%mY(i,j) = in_mY(i,j)
                    END DO
                  END DO    
              
                  DO i = 1, nk
                    DO j = 1, nft
                        set%mK(j,i) = in_mK(i,j)
                    END DO
                  END DO          
              
                  DO i = 1, nk         
                     set%mC(i) = in_mC(i)
                  END DO
              
                  CALL failures(set)
                    
           
           END SUBROUTINE allocate_matrices_cox
           !------------------------------------------------------------------------------------------
           
    
           !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_matrices_mse(set, in_mX, in_mY, in_mK, in_mC,  & 
                                            & nrecord, nft, nk)
               !
               ! Allocates the data matrices in the mean square error model.
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set                       ! The set of information        

               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               REAL(KIND=c_double), DIMENSION(nrecord) :: in_mY           ! responses 
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          

               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! the number of kits
                
               !**************************** OTHER VARIABLES **************************************           
               INTEGER :: i, j               

              
                  DO j = 1, nrecord
                    DO i = 1, nft
                        set%mX(i,j) = in_mX(i,j)
                    END DO
                  END DO
                      
                  DO j = 1, nrecord
                        set%mY_mse(j) = in_mY(j)
                  END DO    
              
                  DO i = 1, nk
                    DO j = 1, nft
                        set%mK(j,i) = in_mK(i,j)
                    END DO
                  END DO          
              
                  DO i = 1, nk         
                     set%mC(i) = in_mC(i)
                  END DO
              
           
           END SUBROUTINE allocate_matrices_mse
           !------------------------------------------------------------------------------------------
 
 
 
            !------------------------------------------------------------------------------------------
           SUBROUTINE allocate_matrices_log(set, in_mX, in_mY, in_mK, in_mC,  & 
                                            & nrecord, nft, nk)
               !
               ! Allocates the data matrices in the logistic regression model.
               ! 
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0
               !         * 'nrecord' > 0               
               !         * 'nkits' > 0               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set                       ! The set of information        

               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrecord) :: in_mY           ! responses 
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          

               INTEGER(KIND=c_int), INTENT(IN) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int), INTENT(IN) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int), INTENT(IN) :: nk                      ! the number of kits
                
               !**************************** OTHER VARIABLES **************************************           
               INTEGER :: i, j               

              
                  DO j = 1, nrecord
                    DO i = 1, nft
                        set%mX(i,j) = in_mX(i,j)
                    END DO
                  END DO
                      
                  DO j = 1, nrecord
                        set%mY_log(j) = in_mY(j)
                  END DO    
              
                  DO i = 1, nk
                    DO j = 1, nft
                        set%mK(j,i) = in_mK(i,j)
                    END DO
                  END DO          
              
                  DO i = 1, nk         
                     set%mC(i) = in_mC(i)
                  END DO
             
           
           END SUBROUTINE allocate_matrices_log
           !------------------------------------------------------------------------------------------
 

 
           !------------------------------------------------------------------------------------------
           SUBROUTINE set_k(set, nk)
               !
               ! Changes the value of k (defines the L0-norm used)
               ! 
               ! NOTICE: * 'nft' >= 'nk' > 0              
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information               
               
               INTEGER(KIND=c_int), INTENT(IN) :: nk         ! defines the k-norm used
               !**************************** OTHER VARIABLES **************************************
               
                  set%nk0 = nk
                       
           END SUBROUTINE set_k
           !------------------------------------------------------------------------------------------
           
           
           !------------------------------------------------------------------------------------------
           SUBROUTINE deallocate_data_cox(set)
               !
               ! Deallocates the data matrices for the Cox's proportional hazard model.           
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               
               DEALLOCATE(set%mX,set%mY,set%mK,set%mC,set%k_norm_ind,& 
                         & set%k_norm_ind_k,set%mFail,set%mUnique,set%hav)
           
           END SUBROUTINE deallocate_data_cox
           !------------------------------------------------------------------------------------------
           
                    
           !------------------------------------------------------------------------------------------
           SUBROUTINE deallocate_data_mse(set)
               !
               ! Deallocates the data matrices for the mean square error model.              
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               
               DEALLOCATE(set%mX,set%mY_mse,set%mK,set%mC, &
                         & set%k_norm_ind,set%k_norm_ind_k,set%hav)
           
           END SUBROUTINE deallocate_data_mse
           !------------------------------------------------------------------------------------------
           
                    
           !------------------------------------------------------------------------------------------
           SUBROUTINE deallocate_data_log(set)
               !
               ! Deallocates the data matrices for the logistic regression model.              
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               
               DEALLOCATE(set%mX,set%mY_log,set%mK,set%mC,& 
                           & set%k_norm_ind,set%k_norm_ind_k,set%hav)
           
           END SUBROUTINE deallocate_data_log
           !------------------------------------------------------------------------------------------
           

           
           !******************************************************************************************
           !                                                                                         *
           !                 **** COX's proportional hazards model - SCALING ****                    *
           !                                                                                         *
           !******************************************************************************************
           !
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_cox(set)
               !
               ! Scales the used data in the Cox's proportional hazards model.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
                INTEGER(KIND=c_int) :: i, j         
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double) :: b                              ! help variable
               
               ! mean input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + set%mX(i,j)
                  END DO
                  b = 1.0_c_double / set%nrecord0
                  set%hav(1,i) = a * b
               
               END DO
               
               ! deviaition input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + (set%mX(i,j)-set%hav(1,i))**2
                  END DO
                  b = 1.0_c_double / (set%nrecord0)
                  a = a * b
                  set%hav(2,i) = SQRT(a)               
                  
               END DO       
               
               ! Scaling of input
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = (set%mX(i,j)-set%hav(1,i))/set%hav(2,i)
                  END DO
               END DO
               
           END SUBROUTINE scaling_cox   
           !------------------------------------------------------------------------------------------        

           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_matrix_cox(set, in_mX)
               !
               ! Stores scaled matrix of info_set 'set' to vector 'in_mX'.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(IN) :: set          ! The set of information            
                               
               REAL(KIND=c_double), DIMENSION(set%nft0,set%nrecord0), INTENT(INOUT) :: in_mX       ! scaled predictor matrix (column is an observation)
               
               !**************************** OTHER VARIABLES **************************************
               INTEGER ::  i, j               

                  DO j = 1, set%nrecord0
                    DO i = 1, set%nft0
                       in_mX(i,j) = set%mX(i,j) 
                    END DO
                  END DO
                  
               
           END SUBROUTINE scaling_matrix_cox   
           !------------------------------------------------------------------------------------------    
           
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_cox(set)
               !
               ! Rescales the data to the original data in the Cox's proportional hazards model.                
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, j                      ! help variables
   
               ! Rescaling of the data
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = set%mX(i,j)*set%hav(2,i)+set%hav(1,i)
                  END DO
               END DO
                       
           END SUBROUTINE rescaling_cox     
           !------------------------------------------------------------------------------------------                      
                     
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_beta_cox(set,y)
               !
               ! Rescales the solution point y to the original data in the Cox's proportional hazards model.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: y      ! a point which is rescaled
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i                     ! help variables
           
               ! Rescaling of the solution y
               DO i = 1, set%nft0
                  y(i) = y(i) / set%hav(2,i)
               END DO
               
           END SUBROUTINE rescaling_beta_cox               
           !------------------------------------------------------------------------------------------    
           



           !******************************************************************************************
           !                                                                                         *
           !                     **** Mean square error model - SCALING ****                         *
           !                                                                                         *
           !******************************************************************************************
           !
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_mse(set)
               !
               ! Scales the used data in the mean square error model.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, j         
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double) :: b                              ! help variable
               
               ! means input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + set%mX(i,j)
                  END DO
                  b = 1.0_c_double / set%nrecord0
                  set%hav(1,i) = a * b
               
               END DO
               
               ! deviations input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + (set%mX(i,j)-set%hav(1,i))**2
                  END DO
                  b = 1.0_c_double / (set%nrecord0)
                  a = a * b
                  set%hav(2,i) = SQRT(a)               
                  
               END DO       
               
               ! Scaling of input
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = (set%mX(i,j)-set%hav(1,i))/set%hav(2,i)
                  END DO
               END DO
               
               ! mean output
             !  a = 0.0_c_double
             !  DO j = 1, nrecord0
             !     a = a + mY_mse(j)
             !  END DO
             !  b = 1.0_c_double / nrecord0
             !  hav(1,nft0+1) = a * b
             
               ! deviations output
             !  a = 0.0_c_double
             !  DO j = 1, nrecord0
             !      a = a + (mY_mse(j)-hav(1,nft0+1))**2
             !  END DO
             !  b = 1.0_c_double / (nrecord0)
             !  a = a * b
             !  hav(2,nft0+1) = SQRT(a)               
                                 
               ! Scaling of output
             !  DO j = 1, nrecord0
             !     mY_mse(j) = (mY_mse(j)-hav(1,nft0+1))/hav(2,nft0+1)
             !  END DO
               
           END SUBROUTINE scaling_mse  
           !------------------------------------------------------------------------------------------   

           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_matrix_mse(set, in_mX)
               !
               ! Stores scaled matrix of info_set 'set' to vector 'in_mX'.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(IN) :: set          ! The set of information            
                               
               REAL(KIND=c_double), DIMENSION(set%nft0,set%nrecord0), INTENT(INOUT) :: in_mX       ! scaled predictor matrix (column is an observation)
               
               !**************************** OTHER VARIABLES **************************************
               INTEGER ::  i, j               

                  DO j = 1, set%nrecord0
                    DO i = 1, set%nft0
                       in_mX(i,j) = set%mX(i,j) 
                    END DO
                  END DO
                   
               
           END SUBROUTINE scaling_matrix_mse   
           !------------------------------------------------------------------------------------------     

           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_mse(set)
               !
               ! Rescales the data to the original data in the mean square error model.                
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, j                      ! help variables
   
               ! Rescaling of the input data
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = set%mX(i,j)*set%hav(2,i)+set%hav(1,i)
                  END DO
               END DO
                                 
           END SUBROUTINE rescaling_mse     
           !------------------------------------------------------------------------------------------             
           
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_beta_mse(set,y)
               !
               ! Rescales the solution point y to the original data in the mean square error model.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: y      ! a point which is rescaled
               !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double), DIMENSION(set%nft0+1) :: x                ! a point which is rescaled             
               INTEGER(KIND=c_int) :: i, num                               ! help variables
                
               x = 0.0_c_double
               
               ! Rescaling of the solution y
               DO i = 1, set%nft0
                  x(i) = y(i) / set%hav(2,i)
                  !x(i) = y(i)*hav(2,nft0+1)/ hav(2,i)
               END DO

               num = set%nft0
               DO i = 1, num
                  x(num+1) = x(num+1) - set%hav(1,i)*y(i)/ set%hav(2,i)
               END DO
               x(num+1) = x(num+1) + y(num+1)         
               y = x

               !DO i = 1, nft0
               !   x(nft0+1) = x(nft0+1) - hav(1,i)*hav(2,nft0+1)*y(i)/ hav(2,i)
               !END DO
               !x(nft0+1) = x(nft0+1) + hav(2,nft0+1)*y(nft0+1) + hav(1,nft0+1)             
               !y = x
               
           END SUBROUTINE rescaling_beta_mse               
           !------------------------------------------------------------------------------------------
           

           !******************************************************************************************
           !                                                                                         *
           !                   **** Logistic regression model - SCALING  ****                        *
           !                                                                                         *
           !******************************************************************************************
           !
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_log(set)
               !
               ! Scales the used data in the logistic regression model.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, j         
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double) :: b                              ! help variable
               
               ! means input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + set%mX(i,j)
                  END DO
                  b = 1.0_c_double / set%nrecord0
                  set%hav(1,i) = a * b
               
               END DO
               
               ! deviations input
               DO i = 1, set%nft0
                  a = 0.0_c_double
                  DO j = 1, set%nrecord0
                      a = a + (set%mX(i,j)-set%hav(1,i))**2
                  END DO
                  b = 1.0_c_double / (set%nrecord0)
                  a = a * b
                  set%hav(2,i) = SQRT(a)               
                  
               END DO       
               
               ! Scaling of input
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = (set%mX(i,j)-set%hav(1,i))/set%hav(2,i)
                  END DO
               END DO
               
               ! mean output
             !  a = 0.0_c_double
             !  DO j = 1, nrecord0
             !     a = a + mY_mse(j)
             !  END DO
             !  b = 1.0_c_double / nrecord0
             !  hav(1,nft0+1) = a * b
             
               ! deviations output
             !  a = 0.0_c_double
             !  DO j = 1, nrecord0
             !      a = a + (mY_mse(j)-hav(1,nft0+1))**2
             !  END DO
             !  b = 1.0_c_double / (nrecord0)
             !  a = a * b
             !  hav(2,nft0+1) = SQRT(a)               
                                 
               ! Scaling of output
             !  DO j = 1, nrecord0
             !     mY_mse(j) = (mY_mse(j)-hav(1,nft0+1))/hav(2,nft0+1)
             !  END DO
               
           END SUBROUTINE scaling_log  
           !------------------------------------------------------------------------------------------                      

           !------------------------------------------------------------------------------------------                      
           SUBROUTINE scaling_matrix_log(set, in_mX)
               !
               ! Stores scaled matrixs of info_set 'set' to vector 'in_mX'.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(IN) :: set          ! The set of information            
                               
               REAL(KIND=c_double), DIMENSION(set%nft0,set%nrecord0), INTENT(INOUT) :: in_mX       ! scaled predictor matrix (column is an observation)
               
               !**************************** OTHER VARIABLES **************************************
               INTEGER ::  i, j               

                  DO j = 1, set%nrecord0
                    DO i = 1, set%nft0
                       in_mX(i,j) = set%mX(i,j) 
                    END DO
                  END DO
                   
               
           END SUBROUTINE scaling_matrix_log   
           !------------------------------------------------------------------------------------------ 

           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_log(set)
               !
               ! Rescales the data to the original data in the logistic regression model.                
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: i, j                      ! help variables
   
               ! Rescaling of the input data
               DO i = 1, set%nft0
                  DO j = 1, set%nrecord0
                      set%mX(i,j) = set%mX(i,j)*set%hav(2,i)+set%hav(1,i)
                  END DO
               END DO
               
               ! Rescaling of the output
              ! DO j = 1, nrecord0
              !    mY_mse(j) = mY_mse(j)*hav(2,nft0+1)+hav(1,nft0+1)
              ! END DO 
                       
           END SUBROUTINE rescaling_log     
           !------------------------------------------------------------------------------------------                      
           
           !------------------------------------------------------------------------------------------                      
           SUBROUTINE rescaling_beta_log(set, y)
               !
               ! Rescales the solution point y to the original data in the logistic regression model.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: y      ! a point which is rescaled
               !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double), DIMENSION(set%nft0+1) :: x                ! a point which is rescaled             
               INTEGER(KIND=c_int) :: i, num                               ! help variables
                
               x = 0.0_c_double
               
               num = set%nft0
               ! Rescaling of the solution y
               DO i = 1, num
                  x(i) = y(i) / set%hav(2,i)
                  !x(i) = y(i)*hav(2,nft0+1)/ hav(2,i)
               END DO

               DO i = 1, num
                  x(num+1) = x(num+1) - set%hav(1,i)*y(i)/ set%hav(2,i)
               END DO
               x(num+1) = x(num+1) + y(num+1)         
               y = x

               !DO i = 1, nft0
               !   x(nft0+1) = x(nft0+1) - hav(1,i)*hav(2,nft0+1)*y(i)/ hav(2,i)
               !END DO
               !x(nft0+1) = x(nft0+1) + hav(2,nft0+1)*y(nft0+1) + hav(1,nft0+1)             
               !y = x
               
           END SUBROUTINE rescaling_beta_log              
           !------------------------------------------------------------------------------------------   

           

           !------------------------------------------------------------------------------------------
           SUBROUTINE failures(set)
               !
               ! Allocates the failure matrix and initializes its values in the Cox's proportional hazards model.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
               !**************************** OTHER VARIABLES **************************************
               INTEGER(KIND=c_int) :: a1, a2, m, eka, fsum, d, place           
               INTEGER(KIND=c_int) :: i, unique, ind 
               LOGICAL :: cont             
               
               m = 0
               DO i = set%nrecord0,1,-1
                 IF (set%mY(2,i)==1) THEN 
                   m = m+1
                   eka = i
                 END IF
               END DO
               
               set%nfail = m
               
               ALLOCATE(set%mFail(2,set%nfail))         
               
               set%mFail = 0
               
               m = 0
               a1 = set%mY(1,eka)            ! Time for the first failure
               eka = 1                   ! The place of the first failure
               fsum = 0                  ! One failure at this time
               unique = 1
               DO i = 1, set%nrecord0
                  IF (set%mY(2,i)==1) THEN
                      m = m + 1
                      set%mFail(1,m) = i     ! Index for the observation i                
                      a2 = set%mY(1,i)       ! Time for the failure at observation i
                      IF (a1 < a2 ) THEN 
                         set%mFail(2,eka) = fsum
                         a1 = a2
                         fsum = 1
                         eka = m
                         unique = unique + 1  ! one more unique failure time
                      ELSE
                        fsum = fsum + 1                   
                      END IF
                   END IF
               END DO
               set%mFail(2,eka)=fsum
               set%nfailunique = unique
               
               ALLOCATE(set%mUnique(2,set%nfailunique))
 
               ind = 1
               DO i = 1, set%nfailunique
                  
                  place = set%mFail(1,ind)
                  d = set%mFail(2,ind)
                  cont = .TRUE. 
                  
                  DO WHILE (cont)
                    IF (place > 1) THEN
                      IF (set%mY(1,place-1)==set%mY(1,place)) THEN
                         place = place - 1
                      ELSE
                         set%mUnique(1,i) = place 
                         set%mUnique(2,i) = d
                         cont = .FALSE.                  
                      END IF
                    ELSE
                     set%mUnique(1,i) = place
                     set%mUnique(2,i) = d
                     cont = .FALSE.                                      
                    END IF                    
                  END DO

                  ind = ind + d 
               END DO              
           
           END SUBROUTINE failures
          !------------------------------------------------------------------------------------------

        
        !********************************************************************************
        !                                                                               |
        !         FUNCTION VALUE and SUBGRADIENT OF THE DC COMPONENT f_1                |
        !                                                                               |
        !********************************************************************************

           SUBROUTINE f1_sub(set, y, problem1, user_n, f, grad)        
                !
                ! Calculates the function value and subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                ! The dimension of the problem is 'user_n'. 
                ! The function value is stored to 'f'
                ! The subgradient is stored to 'grad'
                !
                ! NOTICE: The dimension of 'y' and 'grad' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y      ! a point where the function value of the DC component f_1 is calculated
                REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: grad  ! a subgradient of the DC component f_1 is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: problem1             ! the objective function f_1 for which the value is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: user_n               ! the dimension of the problem
                REAL(KIND=c_double), INTENT(OUT) :: f                   ! the function value of the DC component f_1 at a point 'y'               
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double), DIMENSION(set%nft0) :: mG            ! help variable                        
                REAL(KIND=c_double) :: sum_r, apu, div                ! help variables
                REAL(KIND=c_double) :: apu1, apu2                       ! help variables
                REAL(KIND=c_double) :: exp_term                       ! help variables
                REAL(KIND=c_double) :: largest, reaali                ! help variables
                INTEGER(KIND=c_int) :: time1, time2                   ! help variables                
                INTEGER(KIND=c_int) :: i, j, k, ind, d, place         ! help variable
                INTEGER(KIND=c_int) :: ind1, ind2                     ! help variable
                LOGICAL :: use_log
                
                SELECT CASE(problem1)

                   !-------------------------------------
                   !           Problem   1
                   !        COX model + L1  (Lasso)
                   !-------------------------------------
                   CASE(1)  
                      f = 0.0_c_double
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
                      
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
            
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                      
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
                       
                            f = f + d * Log(sum_r)
                      
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) =  set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
    
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_lambda * y(i)
                         grad(i) = grad(i) + set%user_lambda
                       ELSE
                         f = f - set%user_lambda * y(i)
                         grad(i) = grad(i) - set%user_lambda
                      END IF
                     END DO    
                     
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   2
                   !    COX model + L0 for features (parametric)                   
                   !-------------------------------------
                   CASE(2)  
                      f = 0.0_c_double
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
     
                    ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
         
                            f = f + d * Log(sum_r)
                      
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) = set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
 
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         f = f - set%user_rho * y(i)
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
    
                   !-------------------------------------    

                   !-------------------------------------
                   !           Problem   3
                   !    COX model + L0 for kits (parametric)                   
                   !-------------------------------------
                   CASE(3)  
                      f = 0.0_c_double
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
             
                    ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
                 
                            f = f + d * Log(sum_r)
                      
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) =  set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
                 
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         f = f - set%user_rho * y(i)
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
 
                   !-------------------------------------               

 
                   !-------------------------------------
                   !           Problem   4
                   !    MSE + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(4)  
                      f = 0.0_c_double
                      grad = 0.0_c_double
                      
                      ! mean square error
                      DO i = 1, set%nrecord0
                      
                         apu = 0.0_c_double
                         apu = set%mY_mse(i)-y(user_n)
                         DO j = 1, set%nft0
                            apu = apu - y(j)*set%mX(j,i)
                         END DO
                         
                         f = f + apu**2         
                         grad(user_n) = grad(user_n) - 2.0_c_double * apu
                         DO j = 1, set%nft0
                            grad(j) = grad(j) - 2.0_c_double*set%mX(j,i)*apu
                         END DO
                         
                      END DO 
                      
                     div = 0.5_c_double/set%nrecord0
                     f = f * div
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         f = f - set%user_rho * y(i)
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
                 
                   !-------------------------------------  
                   
                    
                   !-------------------------------------
                   !           Problem   5
                   !    Logistic + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(5)  
                      f = 0.0_c_double
                      grad = 0.0_c_double
                      
                      apu1 = 0.0_c_double
                      apu2 = 0.0_c_double
                      
                      ! The logistic regression model
                      DO i = 1, set%nrecord0
                         
                       ! The linear term
                         apu = y(user_n)
                         DO j = 1, set%nft0
                            apu = apu + y(j)*set%mX(j,i)
                         END DO
                         
                         f = f + set%mY_log(i)*apu
                         f = f - Log(1.0_c_double+Exp(apu))
                         
                         grad(user_n) = grad(user_n) + set%mY_log(i)
                         exp_term = Exp(apu)/(1.0_c_double+Exp(apu))
                         
                         IF (exp_term > 100000.0_c_double) THEN
                            WRITE(*,*) 'exp_term', exp_term
                         END IF 
                         
                         grad(user_n) = grad(user_n) - exp_term
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + set%mX(j,i)*set%mY_log(i)
                            grad(j) = grad(j) - set%mX(j,i)*exp_term
                         END DO
                         
                      END DO 

                     div = 1.0_c_double/set%nrecord0
                     f = - f * div
                     grad = - grad * div
                        
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         f = f - set%user_rho * y(i)
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
                 
                   !-------------------------------------              
 
                END SELECT
                
           END SUBROUTINE f1_sub      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

        !********************************************************************************
        !                                                                               |
        !                   FUNCTION VALUE OF THE DC COMPONENT f_1                      |
        !                                                                               |
        !******************************************************************************** 
           FUNCTION f1(set, y, problem1, user_n) RESULT(f)           
                !
                ! Calculates the function value of DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                ! The dimension of the problem is 'user_n'.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: problem1                 ! the objective function f_1 for which the value is calculated   
                INTEGER(KIND=c_int), INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double) :: f                              ! the function value of the DC component f_1 at a point 'y'
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double) :: sum_r, apu, div                ! help variables
                REAL(KIND=c_double) :: largest, reaali                ! help variables
                REAL(KIND=c_double) :: testaus                        ! help variables
                INTEGER(KIND=c_int) :: time1, time2                         ! help variables                
                INTEGER(KIND=c_int) :: i, j, k, ind, d, place               ! help variable
                INTEGER(KIND=c_int) :: ind1, ind2                           ! help variable
                LOGICAL :: use_log
                
                SELECT CASE(problem1)

                   !-------------------------------------
                   !           Problem   1
                   !        COX model + L1  (Lasso)
                   !-------------------------------------
                   CASE(1)  
                      f = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
      
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu
                        
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              sum_r = sum_r - apu
                   
                            END DO
                            f = f + d * Log(sum_r)              
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                            END IF                          
                         END DO
        
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_lambda * y(i)
                       ELSE
                         f = f - set%user_lambda * y(i)
                      END IF
                     END DO 
 
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   2
                   !    COX model + L0 for features (parametric)                                   
                   !-------------------------------------
                   CASE(2)  
                      f = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
   
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu
                        
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              sum_r = sum_r - apu
                                  
                            END DO
                            f = f + d * Log(sum_r)              
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                            END IF                          
                         END DO
             
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                       ELSE
                         f = f - set%user_rho * y(i)
                      END IF
                     END DO 
                   
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Problem   3
                   !    COX model + L0 for kits (parametric)                                   
                   !-------------------------------------
                   CASE(3)  
                      f = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               f = f - set%mX(k,place)*y(k)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
                       
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu
                        
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              sum_r = sum_r - apu
                                  
                            END DO
                            f = f + d * Log(sum_r)              
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                            f = f + d * Log(sum_r)
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                            END IF                          
                         END DO
                                                        
                         reaali = 1.0_c_double * set%nft0                                     
                         f = f + d * largest + Log(reaali)
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     f = f * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                       ELSE
                         f = f - set%user_rho * y(i)
                      END IF
                     END DO 
                                
                   !-------------------------------------                  

                   !-------------------------------------
                   !           Problem   4
                   !    MSE + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(4)  
                      f = 0.0_c_double
  
                      ! mean square error
                      DO i = 1, set%nrecord0
                      
                         apu = 0.0_c_double
                         apu = set%mY_mse(i)-y(user_n)
                         DO j = 1, set%nft0
                            apu = apu - y(j)*set%mX(j,i)
                         END DO
                         
                         f = f + apu**2         
                         
                      END DO                      
                      
                     div = 0.5_c_double/set%nrecord0
                     f = f * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                       ELSE
                         f = f - set%user_rho * y(i)
                      END IF
                     END DO 
                     
                 
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   5
                   !    Logistic + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(5)  
                      f = 0.0_c_double
                      
                      ! The logistic regression model
                      DO i = 1, set%nrecord0
                         
                       ! The linear term
                         apu = y(user_n)
                         DO j = 1, set%nft0
                            apu = apu + y(j)*set%mX(j,i)
                         END DO
                         
                         f = f + set%mY_log(i)*apu
                         f = f - Log(1.0_c_double+Exp(apu))
                        
                      END DO 

                     div = 1.0_c_double/set%nrecord0
                     f = - f * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         f = f + set%user_rho * y(i)
                       ELSE
                         f = f - set%user_rho * y(i)
                      END IF
                     END DO 
                 
                   !-------------------------------------   
                
                END SELECT              
            

           END FUNCTION f1
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        !********************************************************************************
        !                                                                               |
        !                   FUNCTION VALUE OF THE DC COMPONENT f_2                      |
        !                                                                               |
        !********************************************************************************            
           FUNCTION f2(set, y, problem2, user_n) RESULT(f)           
                !
                ! Calculates the function value of DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                ! The dimension of the problem is 'user_n'.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: problem2                 ! the objective function f_2 for which the value is calculated   
                INTEGER(KIND=c_int), INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double) :: f                             ! the function value of the DC component f_2 at a point 'y'
                REAL(KIND=c_double), DIMENSION(user_n) :: absterm    ! the absolute values of the beta vector
                REAL(KIND=c_double), DIMENSION(set%nkits0) :: absterm_k  ! the absolute values of the beta vector
                INTEGER(KIND=c_int), DIMENSION(user_n) :: absind           ! the indices of the absolute values
                INTEGER(KIND=c_int), DIMENSION(set%nkits0) :: absind_k         ! the indices of the absolute values
                INTEGER(KIND=c_int) :: i, j                                ! help variables
                

                
                SELECT CASE(problem2)
                
                 
                   !-------------------------------------
                   !           Problem   1
                   !        COX model + L1  (Lasso)                
                   !-------------------------------------
                   CASE(1)  
                      f = 0.0_c_double

                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   2
                   !    COX model + L0 for features (parametric)                                                   
                   !-------------------------------------
                   CASE(2)  
                      f = 0.0_c_double
                      
                      DO i = 1, set%nft0
                         absterm(i) = ABS(y(i))
                         absind(i) = i
                      END DO
                      
                      CALL heapsort_k(absterm,absind,set%nk0)
                      
                      DO i = user_n-set%nk0+1, user_n
                        f = f + absterm(i)
                      END DO
                      f = f * set%user_rho
                      
                      ! the indices of the k largest absolute values of beta
                      set%k_norm_ind = absind
                      
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Problem   3
                   !    COX model + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(3)  
                      f = 0.0_c_double
                      absterm_k = 0.0_c_double
                      
                      DO i = 1, set%nkits0
                        DO j = 1, set%nft0
                          absterm_k(i) = absterm_k(i) + set%mK(j,i) * ABS(y(j))
                        END DO  
                        absind_k(i) = i
                      END DO
                      
                      CALL heapsort_k(absterm_k,absind_k,set%nk0)
                      
                      DO i = set%nkits0-set%nk0+1, set%nkits0
                        f = f + absterm_k(i)
                      END DO
                      f = f * set%user_rho
                      
                      ! the indices of the k largest absolute values of beta for kit structure
                      set%k_norm_ind_k = absind_k
                      
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   4
                   !    MSE + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(4)  
                      f = 0.0_c_double
                      absterm_k = 0.0_c_double
                      
                      DO i = 1, set%nkits0
                        DO j = 1, set%nft0
                          absterm_k(i) = absterm_k(i) + set%mK(j,i) * ABS(y(j))
                        END DO  
                        absind_k(i) = i
                      END DO
                      
                      CALL heapsort_k(absterm_k,absind_k,set%nk0)
                      
                      DO i = set%nkits0-set%nk0+1, set%nkits0
                        f = f + absterm_k(i)
                      END DO
                      f = f * set%user_rho
                      
                      ! the indices of the k largest absolute values of beta for kit structure
                      set%k_norm_ind_k = absind_k
                      
                   !------------------------------------- 

                   !-------------------------------------
                   !           Problem   5
                   !    Logistic + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(5)  
                      f = 0.0_c_double
                      absterm_k = 0.0_c_double
                      
                      DO i = 1, set%nkits0
                        DO j = 1, set%nft0
                          absterm_k(i) = absterm_k(i) + set%mK(j,i) * ABS(y(j))
                        END DO  
                        absind_k(i) = i
                      END DO
                      
                      CALL heapsort_k(absterm_k,absind_k,set%nk0)
                      
                      DO i = set%nkits0-set%nk0+1, set%nkits0
                        f = f + absterm_k(i)
                      END DO
                      f = f * set%user_rho
                      
                      ! the indices of the k largest absolute values of beta for kit structure
                      set%k_norm_ind_k = absind_k
                      
                   !-------------------------------------                  
                   
                
                END SELECT              
            


           END FUNCTION f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           

        !********************************************************************************
        !                                                                               |
        !                    SUBGRADIENT OF THE DC COMPONENT f_1                        |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_f1(set, y, problem1, user_n) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                ! The dimension of the problem is 'user_n'.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_1 is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: problem1                 ! the objective function f_1 for which the subgradient is calculated      
                INTEGER(KIND=c_int), INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'
                REAL(KIND=c_double) :: a                              ! help variable
                REAL(KIND=c_double), DIMENSION(set%nft0) :: mG            ! help variable                        
                REAL(KIND=c_double) :: sum_r, apu, div                ! help variables
                REAL(KIND=c_double) :: exp_term                       ! help variables
                REAL(KIND=c_double) :: largest                        ! help variables
                INTEGER(KIND=c_int) :: time1, time2                   ! help variables                
                INTEGER(KIND=c_int) :: i, j, k, ind, d, place         ! help variable
                INTEGER(KIND=c_int) :: ind1, ind2                     ! help variable
                LOGICAL :: use_log
                
                SELECT CASE(problem1)

                   !-------------------------------------
                   !           Problem   1
                   !        COX model + L1  (Lasso)
                   !-------------------------------------
                   CASE(1)  
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
          
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
                            
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) =  set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
                         
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         grad(i) = grad(i) + set%user_lambda
                       ELSE
                         grad(i) = grad(i) - set%user_lambda
                      END IF
                     END DO 
                     
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   2
                   !    COX model + L0 for features (parametric)                                                   
                   !-------------------------------------                    
                   CASE(2)              
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
 
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
                   
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) =  set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
 
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
                
                   !-------------------------------------             

                   !-------------------------------------
                   !           Problem   3
                   !    COX model + L0 for kits (parametric)                                                   
                   !-------------------------------------                    
                   CASE(3)              
                      grad = 0.0_c_double
                      
                      ! linear term in Cox model
                      ind = 1                       ! Position of the first failure at time t_1 in mFail
                      DO i = 1, set%nfailunique
                         d = set%mFail(2,ind)           ! the number of failures at time t_ind
                         DO j = ind, ind+d-1
                            place = set%mFail(1,j)      ! Position of the failure in mX
                            DO k = 1, set%nft0
                               grad(k) = grad(k) - set%mX(k,place)
                            END DO                          
                         END DO
                         ind = ind + d
                      END DO
                      
                      use_log = .TRUE.   
                      apu = 0.0_c_double
                      DO j = 1, set%nft0
                          apu = apu + set%mX(j,1)*y(j)
                      END DO
                       
                      a = set%user_a 
                      IF (ABS(apu)>=set%user_a) THEN 
                        IF (ABS(apu)<(set%user_a+1.0_c_double)) THEN
                           a = set%user_a + 1.0_c_double
                        END IF
                      END IF    
                         
                      IF (ABS(apu)>=a) THEN 
                         use_log = .FALSE.            ! ln-exp-term cannot be used, instead maximum is used
                      END IF 
                    
                      ! ln-exp-term in Cox model                      
                      IF (use_log) THEN
                      
                      mG = 0.0_c_double
                      sum_r = 0.0_c_double
                      
                      DO i = 1, set%nrecord0
                         apu = 0.0_c_double
                         DO j = 1, set%nft0
                           apu = apu + set%mX(j,i)*y(j) 
                         END DO     

                         apu = Exp(apu)                  
                         sum_r = sum_r + apu

                         DO j = 1, set%nft0
                            mG(j) = mG(j) + set%mX(j,i) * apu                    
                         END DO                         
                      END DO
                      
                      ind = 1                              ! The first failure happens at time t_1
                      time1 = 1       
                     
                      DO k = 1, set%nfailunique
                         time2 = set%mUnique(1,k)              ! The first failure at time t_ind 
                         d = set%mFail(2,ind)                  ! The number of failures for time t_ind
                         IF (time1<time2) THEN 
                            DO i = time1, time2-1
                              apu = 0.0_c_double
                              DO j = 1, set%nft0
                                apu = apu + set%mX(j,i)*y(j) 
                              END DO       
                              apu = Exp(apu)
                              
                              sum_r = sum_r - apu
                              
                              DO j = 1, set%nft0
                                 mG(j) = mG(j) - set%mX(j,i) * apu                   
                              END DO                                  
                            END DO
                            
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2 
                         ELSE 
                  
                            div = 1.0_c_double / sum_r
                            DO j = 1, set%nft0
                               grad(j) = grad(j) + d * mG(j) * div
                            END DO
                        
                            ind = ind + d
                            time1 = time2                         
                         END IF 
                      END DO
                      ELSE
                      ! maximum term is used instead of ln-exp-term
                      
                      ind2 = set%nrecord0
                      largest = -(10.0_c_double)**10
                      DO k = set%nfailunique, 1, -1
                         ind1 = set%mUnique(1,k)               ! The first failure at time t_k 
                         d = set%mUnique(2,k)                  ! The number of failures for time t_k
                         DO i = ind1, ind2
                            apu = 0.0_c_double
                            DO j = 1, set%nft0
                              apu = apu + set%mX(j,i)*y(j) 
                            END DO  
                            IF (apu > largest) THEN         
                                largest = apu
                                DO j = 1, set%nft0
                                 mG(j) =  set%mX(j,i)                
                                END DO      
                            END IF                          
                         END DO
                       
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + d * mG(j) 
                         END DO
                        
                         ind2 = ind1-1
                     END DO
                      
                     END IF
                      
                     div = 2.0_c_double/set%nrecord0
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
               
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Problem   4
                   !    MSE + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(4)  
                      grad = 0.0_c_double
                      
                      ! mean square error
                      DO i = 1, set%nrecord0
                      
                         apu = 0.0_c_double
                         apu = set%mY_mse(i)-y(user_n)
                         DO j = 1, set%nft0
                            apu = apu - y(j)*set%mX(j,i)
                         END DO
                         
                         grad(user_n) = grad(user_n) - 2.0_c_double * apu
                         DO j = 1, set%nft0
                            grad(j) = grad(j) - 2.0_c_double*set%mX(j,i)*apu
                         END DO
                         
                      END DO 
                      
                     div = 0.5_c_double/set%nrecord0
                     grad = grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
                 
                   !-------------------------------------  

                   !-------------------------------------
                   !           Problem   5
                   !    Logistic + L0 for kits (parametric)                   
                   !-------------------------------------
                   ! beta_0 = y(user_n) and beta_1,...,beta_nft = y(1),...,y(nft=user_n-1)
                   CASE(5)  
                      grad = 0.0_c_double
                      
                      ! The logistic regression model
                      DO i = 1, set%nrecord0
                         
                       ! The linear term
                         apu = y(user_n)
                         DO j = 1, set%nft0
                            apu = apu + y(j)*set%mX(j,i)
                         END DO
              
                         grad(user_n) = grad(user_n) + set%mY_log(i)
                         exp_term = Exp(apu)/(1.0_c_double+Exp(apu))
                         grad(user_n) = grad(user_n) - exp_term
                         DO j = 1, set%nft0
                            grad(j) = grad(j) + set%mX(j,i)*set%mY_log(i)
                            grad(j) = grad(j) - set%mX(j,i)*exp_term
                         END DO
                         
                      END DO 

                     div = 1.0_c_double/set%nrecord0
                     grad = - grad * div
                      
                     DO i = 1, set%nft0
                       IF (y(i)>=0.0_c_double) THEN 
                         grad(i) = grad(i) + set%user_rho
                       ELSE
                         grad(i) = grad(i) - set%user_rho
                      END IF
                     END DO 
                 
                   !-------------------------------------                       
 
                END SELECT                      
                
           END FUNCTION subgradient_f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
        !********************************************************************************
        !                                                                               |
        !                    SUBGRADIENT OF THE DC COMPONENT f_2                        |
        !                                                                               |       
        !********************************************************************************    
           FUNCTION subgradient_f2(set, y, problem2, user_n) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                ! The dimension of the problem is 'user_n'.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                TYPE(set_info), INTENT(INOUT) :: set          ! The set of information              
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_2 is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: problem2           ! the objective function f_2 for which the subgradient is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: user_n             ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'
                INTEGER(KIND=c_int) :: i, j, ind                      ! help variables
                
               
                SELECT CASE(problem2)
                
                   
                   !-------------------------------------
                   !           Problem   1
                   !       COX model + L1 (Lasso)                                                  
                   !-------------------------------------                    
                   CASE(1)              
                     grad = 0.0_c_double           

                   !-------------------------------------      
                   
                   !-------------------------------------
                   !           Problem   2
                   !    COX model + L0 for features (parametric)                                                   
                   !-------------------------------------                    
                   CASE(2)              
                      grad = 0.0_c_double    

                      DO i = user_n-set%nk0+1, user_n
                        ind = set%k_norm_ind(i)
                        IF (y(ind)>=0.0_c_double) THEN 
                          grad(ind) = grad(ind) + set%user_rho
                        ELSE
                          grad(ind) = grad(ind) - set%user_rho                      
                        END IF
                      END DO

                   !-------------------------------------     

                   !-------------------------------------
                   !           Problem   3
                   !    COX model + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(3)  
                      grad = 0.0_c_double    

                      DO i = set%nkits0-set%nk0+1, set%nkits0
                        ind = set%k_norm_ind_k(i)
                        DO j = 1, set%nft0
                          IF(set%mK(j,ind)==1) THEN
                            IF (y(j)>=0.0_c_double) THEN 
                               grad(j) = grad(j) + set%user_rho
                            ELSE
                               grad(j) = grad(j) - set%user_rho                     
                            END IF                            
                          END IF  
                        END DO
                      END DO

                   !-------------------------------------        
                   
                   !-------------------------------------
                   !           Problem   4
                   !    MSE + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(4)  
                      grad = 0.0_c_double    

                      DO i = set%nkits0-set%nk0+1, set%nkits0
                        ind = set%k_norm_ind_k(i)
                        DO j = 1, set%nft0
                          IF(set%mK(j,ind)==1) THEN
                            IF (y(j)>=0.0_c_double) THEN 
                               grad(j) = grad(j) + set%user_rho
                            ELSE
                               grad(j) = grad(j) - set%user_rho                     
                            END IF                            
                          END IF  
                        END DO
                      END DO

                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem 5 
                   !    Logistic + L0 for kits (parametric)                                                   
                   !-------------------------------------
                   CASE(5)  
                      grad = 0.0_c_double    

                      DO i = set%nkits0-set%nk0+1,set%nkits0
                        ind = set%k_norm_ind_k(i)
                        DO j = 1, set%nft0
                          IF(set%mK(j,ind)==1) THEN
                            IF (y(j)>=0.0_c_double) THEN 
                               grad(j) = grad(j) + set%user_rho
                            ELSE
                               grad(j) = grad(j) - set%user_rho                     
                            END IF                            
                          END IF  
                        END DO
                      END DO

                   !-------------------------------------                   
                END SELECT  

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
           
           
          ! Orders the table a in increasing order
           SUBROUTINE heapsort(a)
 
                REAL(KIND=c_double), INTENT(INOUT) :: a(0:)
                INTEGER(KIND=c_int) :: start, n, bottom
                REAL(KIND=c_double) :: temp
 
                n = size(a)
                DO start = (n - 2) / 2, 0, -1
                   call siftdown(a, start, n)
                END DO
 
                DO bottom = n - 1, 1, -1
                   temp = a(0)
                   a(0) = a(bottom)
                   a(bottom) = temp
                   call siftdown(a, 0, bottom)
                END DO
 
           END SUBROUTINE heapsort
 
           SUBROUTINE siftdown(a, start, bottom)
 
                REAL(KIND=c_double), INTENT(INOUT) :: a(0:)
                INTEGER(KIND=c_int), INTENT(IN) :: start, bottom
                INTEGER(KIND=c_int) :: child, root
                REAL(KIND=c_double) :: temp
 
                root = start
                DO WHILE(root*2 + 1 < bottom)
                   child = root * 2 + 1
    
                   IF (child + 1 < bottom) THEN
                      IF (a(child) < a(child+1)) child = child + 1
                   END IF
 
                   IF (a(root) < a(child)) THEN
                      temp = a(child)
                      a(child) = a(root)
                      a(root) = temp
                      root = child
                   ELSE
                      RETURN
                   END IF
                END DO      
 
           END SUBROUTINE siftdown

           ! Orders the whole table a in increasing order.
           ! Table b will give the original indices of elements in the ordered table a.
           SUBROUTINE heapsort_ind(a,b)
 
                REAL(KIND=c_double), INTENT(INOUT) :: a(0:)
                INTEGER(KIND=c_int), INTENT(INOUT) :: b(0:)
                INTEGER(KIND=c_int) :: start, n, bottom
                REAL(KIND=c_double) :: temp
                INTEGER(KIND=c_int) :: apu
 
                n = size(a)
                DO start = (n - 2) / 2, 0, -1
                   CALL siftdown_ind(a, b, start, n)
                END DO
 
                DO bottom = n - 1, 1, -1
                   !swap in a
                   temp = a(0)
                   a(0) = a(bottom)
                   a(bottom) = temp
                   !swap in b 
                   apu = b(0)
                   b(0) = b(bottom)
                   b(bottom) = apu  
                   CALL siftdown_ind(a, b, 0, bottom)
                END DO
 
           END SUBROUTINE heapsort_ind
 
           SUBROUTINE siftdown_ind(a, b, start, bottom)
 
                REAL(KIND=c_double), INTENT(INOUT) :: a(0:)
                INTEGER(KIND=c_int), INTENT(INOUT) :: b(0:)
                INTEGER(KIND=c_int), INTENT(IN) :: start, bottom
                INTEGER(KIND=c_int) :: child, root
                REAL(KIND=c_double) :: temp
                INTEGER(KIND=c_int) :: apu
 
                root = start
                DO WHILE(root*2 + 1 < bottom)
                    child = root * 2 + 1
  
                    IF (child + 1 < bottom) THEN
                       IF (a(child) < a(child+1)) child = child + 1
                    END IF
 
                    IF (a(root) < a(child)) THEN
                        ! swap in a 
                        temp = a(child)
                        a(child) = a(root)
                        a(root) = temp
                        ! swap in b
                        apu = b(child)
                        b(child) = b(root)
                        b(root) = apu

                        root = child
                    ELSE
                        RETURN
                    END IF  
                END DO      
 
           END SUBROUTINE siftdown_ind

           ! Orders the table a in increasing order such that only k largest values are calculated. 
           ! They can be found from positions n-k+1, ... , n of table a (they are in incresing order). 
           ! This means that the beginning of table a is not ordered, i.e. positions from 1 to n-k.
           ! Table b will give the original indices of elements of table a after it is ordered. 
           SUBROUTINE heapsort_k(a,b,k)
 
                REAL(KIND=c_double), INTENT(INOUT) :: a(0:)
                INTEGER(KIND=c_int), INTENT(INOUT) :: b(0:)
                INTEGER(KIND=c_int), INTENT(IN) :: k 
                INTEGER(KIND=c_int) :: start, n, bottom, pituus
                REAL(KIND=c_double) :: temp
                INTEGER(KIND=c_int) :: apu
 
                n = size(a)
   
                IF (k > n) THEN 
                    WRITE(*,*) 'WARNING: k is larger than the dimension in heapsort -> k is reduced'   
                    pituus = n
                ELSE
                    pituus = k   
                END IF 
   
                DO start = (n - 2) / 2, 0, -1
                    CALL siftdown_ind(a, b, start, n)
                END DO
 
                DO bottom = n - 1, n-pituus+1, -1
                    !swap in a
                    temp = a(0)
                    a(0) = a(bottom)
                    a(bottom) = temp
                    !swap in b 
                    apu = b(0)
                    b(0) = b(bottom)
                    b(bottom) = apu  
                    CALL siftdown_ind(a, b, 0, bottom)
                END DO
   
                IF (pituus < n) THEN 
                    bottom = n-pituus
                    !swap in a
                    temp = a(0)
                    a(0) = a(bottom)
                    a(bottom) = temp
                    !swap in b 
                    apu = b(0)
                    b(0) = b(bottom)
                    b(bottom) = apu    
                END IF
 
           END SUBROUTINE heapsort_k          

!**          
           
        !********************************************************************************
        !                                                                               |
        !                FUNCTION VALUES OF THE FUNCTION f AT A POINT x                 |
        !                                                                               |
        !********************************************************************************

           SUBROUTINE func(user_n, y, f, nproblem, set)     
                !
                ! Calculates the function at a point 'y' for LMBM.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                !               
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value is calculated
                INTEGER(KIND=c_int), INTENT(IN) :: user_n             ! the dimension of the problem
                INTEGER(KIND=c_int), INTENT(IN) :: nproblem           ! the solved problem
                REAL(KIND=c_double), INTENT(OUT) :: f                 ! the function value of f at a point y  
                TYPE(set_info), INTENT(INOUT) :: set                  ! The set of information                  
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=c_double) :: fv1, fv2                       ! Function values of f1 and f2 at y
                 
                 fv1 = f1(set, y, nproblem, user_n)    
                 fv2 = f2(set, y, nproblem, user_n)
                 f = fv1 - fv2
                               
           END SUBROUTINE func   
           
         
        !********************************************************************************
        !                                                                               |
        !                  SUBGRADIENT OF THE FUNCTION f AT A POINT x                   |
        !                                                                               |
        !********************************************************************************
        
           SUBROUTINE subgra(user_n, y, grad, nproblem, set)          
                !
                ! Calculates the subgradient of function f at a point 'x' for LMBM.
                !
                ! NOTICE: The dimension of 'x' and 'grad' has to be 'user_n'.

                
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: y         ! a point where the subgradient of the function f is calculated
                REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: grad   ! subgradient
                INTEGER(KIND=c_int), INTENT(IN) :: user_n         ! the dimension of the problem
                INTEGER(KIND=c_int), INTENT(IN) :: nproblem       ! the solved problem
                TYPE(set_info), INTENT(INOUT) :: set                   ! The set of information                                 
                !**************************** OTHER VARIABLES **************************************             
                REAL(KIND=c_double), DIMENSION(user_n) :: grad1, grad2     ! the subgradients of f1 and f2
      
                 grad1 = subgradient_f1(set, y, nproblem, user_n)           
                 grad2 = subgradient_f2(set, y, nproblem, user_n)           
                 grad = grad1 - grad2 

           END SUBROUTINE subgra
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                   
!**          
           
      END MODULE functions     
       

      MODULE dbdc
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |         THE PROXIMAL DOUBLE BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION          | | 
        !| |                                 (version 2)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |      Features :                                                                  | |
        !| |                                                                                  | |
        !| |           * Possibility to use simple stepsize determination after               | |
        !| |             each 'main iteration'.                                               | |
        !| |                                                                                  | |        
        !| |                                                                                  | |                
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|    Utilizes the new version of PLQDF1 by Ladislav Luksan as a quadratic solver.      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*


        USE, INTRINSIC :: iso_c_binding
       
        USE bundle1                    ! The BUNDLE of the DC component f_1
        USE bundle2                    ! The BUNDLE of the DC component f_2
        USE functions                  ! Contains INFORMATION from the USER
        
        IMPLICIT NONE    
        
        EXTERNAL PLQDF1             ! The QUADRATIC SOLVER by Ladislav Luksan
        
        CONTAINS
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      DBDC_algorithm      : The double bundle method for DC optimization.         | | 
        !| |      main_iteration      : The main iteration algorithm needed in the double     | |
        !| |                            bundle method for DC optimization.                    | |
        !| |      guaranteeing_clarke : The algorithm guaranteeing Clarke stationarity        | |       
        !| |                            for a solution obtained.                              | |
        !| |                                                                                  | |       
        !| |      quadratic_solver    : The solver for the quadratic norm minimization        | |       
        !| |                            problem. (Needed in the Clarke stationary algorithm)  | |       
        !| |      subproblem_solver   : The solver for subproblems in the search direction    | |
        !| |                            problem. (Needed in the main iteration algorithm)     | |   
        !| |                                                                                  | |   
        !| |                                                                                  | |       
        !| |                            CONTAINS FUNCTIONS:                                   | |
        !| |                                                                                  | |
        !| |        select_value_t : Selects the parameter t from interval [t_min,t_max].     | |               
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
        !  |           THE DOUBLE BUNDLE ALGORITHM FOR NONSMOOTH DC OPTIMIZATION            |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE DBDC_algorithm( f_solution, x_solution, x_0, rho, lambda, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUTime,  &
                            & agg_used, stepsize_used, iprint, problem1, problem2, user_n, &
                            & set)
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'problem1'       : The f1 in problem
            !         * 'problem2'       : The f2 in problem     
            !         * 'user_n'         : The dimension of the problem   
            !         * 'mit'            : The maximum number of 'main iterations'      
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration'
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Clarke stationary' algorithm
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm          
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm            
            !         * 'iprint'         : Specifies the print      
            !               
            !         * 'rho'            : the penalization parameter for L0-norm
            !         * 'lambda'         : the penalization parameter for L1-norm
            !
            !         * 'x_0'            : the starting point 
            !         
            !
            ! OUTPUT: * 'f_solution' : The objective function value at the solution 'x_solution'  
            !         * 'x_solution' : The solution obtained to the minimizatin problem
            !         * 'termination': The cause of the termination in the algorithm
            !         * 'counter'    : Gives the values of different counters
            !         * 'CPUtime'    : the CPU time
            !
            ! INOUT:
            !
            !         * 'set'            : the set of information 
            !
            ! NOTICE: * The dimension of vectors 'x_0' and 'x_solution' has to be 'nft' ('nft' is defined by USER in MODULE functions)
            !         * The dimension of the vector 'counter' has to be 8.
            !         * 'mit', 'mrounds' and 'mrounds_clarke' have to be integers.
            !         * IF ('mit' <= 0) THEN DEFAULT value 1000 is used
            !         * IF ('mrounds' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * IF ('mrounds_clarke' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * 'iprint' has to be -4, -3, -2, -1, 0, 1, 2, 3 or 4 . If it is NOT then DEFAULT value 1 is used. 
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER ************************************* 
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information               
                
               REAL(KIND=c_double), DIMENSION(user_n), INTENT(OUT) :: x_solution  ! the solution obtained to the problem

               REAL(KIND=c_double), DIMENSION(user_n), INTENT(IN) :: x_0          ! the starting point (the dimension 'user_n' is the number of variables)
 
               REAL(KIND=c_double), INTENT(OUT) :: f_solution             ! the objective function value at the solution 'x_solution'
               REAL(KIND=c_double), INTENT(OUT) :: CPUtime                ! the CPU time 
               
               REAL(KIND=c_double), INTENT(IN) :: rho                     ! the penalization parameter
               REAL(KIND=c_double), INTENT(IN) :: lambda                  ! the penalization parameter
               
               INTEGER(KIND=c_int), INTENT(IN) :: problem1                ! the f1 in problem
               INTEGER(KIND=c_int), INTENT(IN) :: problem2                ! the f2 in problem
               
               INTEGER(KIND=c_int), INTENT(IN) :: user_n                  ! the dimension of the problem

               INTEGER(KIND=c_int), INTENT(INOUT) :: mit                  ! the maximum number of 'main iterations'
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds              ! the maximum number of rounds during one 'main iteration'
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds_clarke       ! the maximum number of rounds during one 'Clarke stationary' algorithm
               
               INTEGER(KIND=c_int), INTENT(OUT) :: termination        ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8), INTENT(OUT) :: counter  ! contains the values of different counteres: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'
                                                              
               LOGICAL, INTENT(IN) :: agg_used          ! .TRUE. if aggregation is used in the algorithm. Otherwise .FALSE.                                                           
               LOGICAL, INTENT(IN) :: stepsize_used     ! .TRUE. if simple stepsize determination is used in the algorithm. Otherwise .FALSE.                                                             
                        
               INTEGER(KIND=c_int), INTENT(INOUT) :: iprint ! variable that specifies print option:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result 
                                                !   iprint = -1: basic print of final result (without the solution vector)
                                                !   iprint = 2 : extended print of final result 
                                                !   iprint = -2: extended print of final result (without the solution vector)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results
                                                !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                                !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                                !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                       
                      
           !***************************** LOCAL VARIABLES ************************************  

               TYPE(kimppu1) :: B1             ! The bundle B_1 for the DC component f_1
               TYPE(kimppu2) :: B2             ! The bundle B_2 for the DC component f_2
               
               REAL(KIND=c_double), DIMENSION(user_n) :: x_current     ! the current iteration point (the dimension 'user_n' is the number of variables)
               REAL(KIND=c_double), DIMENSION(user_n) :: x_new         ! the new iteration point (obtained from the previous 'main iteration')              
           
               REAL(KIND=c_double), DIMENSION(user_n) :: grad1, grad2  ! subgradients (the dimenison 'user_n' is the length of subgradient)       
               REAL(KIND=c_double), DIMENSION(user_n) :: dd            ! the search direction obtained from the 'main iteration': dd=x_new -x_current
                                                                 ! (the dimension 'user_n' is the length of subgradient)               

               REAL(KIND=c_double) :: crit_tol, eps  ! 'crit_tol'=the stopping tolerance and 'eps'=the enlargement parameter
               REAL(KIND=c_double) :: m              ! the descent parameter used in 'main iteration'
               REAL(KIND=c_double) :: r_dec, r_inc   ! 'r_dec'=the decrease parameter and 'R_inc'=the increase parameter  
               REAL(KIND=c_double) :: c              ! 'c'=the decrease parameter
               
               REAL(KIND=c_double) :: step_tol       ! Step-length tolerance used in 'Clarke stationary' algorithm (same as the proximity measure)     
               REAL(KIND=c_double) :: m_clarke       ! the descent parameter used in 'Clarke stationary' algorithm
           
               REAL(KIND=c_double) :: f_0                         ! the value of the objective function f=f_1-f_2 at the starting point x_0
               REAL(KIND=c_double) :: f1_current, f2_current      ! the value of f_1 and f_2 at the current solution  
               REAL(KIND=c_double) :: f1_new, f2_new              ! the value of f_1 and f_2 at the new iteration point x_new (which is obtained from the previous main iteration)        
               REAL(KIND=c_double) :: change                      ! 'change' = f(x_new) - f(x_current)  (i.e. the change in the objective function value)
               REAL(KIND=c_double) :: change1                     ! 'change1' = f1(x_new) - f1(x_current)  (i.e. the change in the DC component f1)
               REAL(KIND=c_double) :: change2                     ! 'change2' = f2(x_new) - f2(x_current)  (i.e. the change in the DC component f2)
               
               REAL(KIND=c_double) :: norm                        ! the distance between two consecutive iteration points        
               
               REAL(KIND=c_double) :: delta                       ! variable which can beused to relax the descent condition          
    
               REAL(KIND=c_double) :: start_time, finish_time                   ! start and finish CPU time
               REAL(KIND=c_double) :: start_time_main_it, finish_time_main_it   ! start and finish CPU time in one 'main iteration'
               
               REAL(KIND=c_double) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER(KIND=c_int) :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time   

               INTEGER(KIND=c_int) :: size_b1 , size_b2    ! The biggest possible size of the bundles B_1 and B_2 
               
               INTEGER(KIND=c_int)  :: iter_counter        ! the number of main iterations executed 
               INTEGER(KIND=c_int)  :: subprob_counter     ! the number of subproblems solved during the execution of the bundle algorithm
               INTEGER(KIND=c_int)  :: stop_cond_counter   ! the number of times 'Clarke stationary' algorithm is used during the algorithm
               INTEGER(KIND=c_int)  :: f_counter           ! the number of function values evaluated for a DC component during 
                                               ! the execution of the bundle algorithm (same for the DC component f_1 and f_2) (without Clarke stationary algorithm)

               INTEGER(KIND=c_int) :: clarke_f_counter     ! the number of function values evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               INTEGER(KIND=c_int) :: clarke_sub_counter   ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               
               INTEGER(KIND=c_int)  :: subgrad1_counter    ! the number of subgradients calculated for f_1 during the bundle algorithm (without Clarke stationary algorithm)
               INTEGER(KIND=c_int)  :: subgrad2_counter    ! the number of subgradients calculated for f_2 during the bundle algorithm (without Clarke stationary algorithm)                 
               
               INTEGER(KIND=c_int) :: help_mit_iter_counter    ! the number of iteration rounds executed in one 'main iteration'
               INTEGER(KIND=c_int) :: help_subprob_counter     ! the number of subproblems solved in one 'main iteration'
               INTEGER(KIND=c_int) :: help_f_counter           ! the number of function values evaluated for a DC component in one 'main iteration' (same for f_1 and f_2)
               INTEGER(KIND=c_int) :: help_subgrad1_counter    ! the number of subgradients calculated for f_1 in one 'main iteration' (without Clarke stationary algorithm)
               INTEGER(KIND=c_int) :: help_subgrad2_counter    ! the number of subgradients calculated for f_2 in one 'main iteration' (without Clarke stationary algorithm)
               INTEGER(KIND=c_int) :: help_stop_cond_counter   ! the number of times CLarke stationarity was tested during the 'main iteration'
               
               INTEGER(KIND=c_int) :: help_clarke_f_counter     ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               INTEGER(KIND=c_int) :: help_clarke_sub_counter   ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               
               INTEGER(KIND=c_int) :: reason_for_stop       ! the reason for stop during the 'main iteration'
                                                ! 0 - a new iteration point found
                                                ! 1 - stopping condition satisfied (Clarke stationarity)
                                                ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)  
                                                ! 3 - the biggest possible number of rounds executed in the 'main iteration'
                                                ! 5 - the biggest possible number of rounds executed in the 'Clarke stationary' algorithm
               
               INTEGER(KIND=c_int) :: i         ! help variables 
                 
               LOGICAL :: stop_alg              ! .TRUE. if the proximal double bundle algorithm DBDC can be terminated                                             

               REAL(KIND=c_double) :: apuf  ! variable which can beused to relax the descent condition          
        
               
               CALL cpu_time(start_time)                ! Start CPU timing     

               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time     
 
               
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
           
               !*************** PARAMETER VALUES FROM USER ARE CHECKED **********************
               !
               !     IF VALUES ARE INCORRECT THEN DEFAULT VALUES ARE USED FOR PARAMETERS
               !
 
                ! ***** penalization parameter 'lambda' for L1-norm *****
               IF ( (lambda<0.0_c_double) ) THEN
                   set%user_lambda = 0.0_c_double
               ELSE       
                   set%user_lambda = lambda 
               END IF
 
               ! ***** penalization parameter 'rho' for L0-norm *****
               IF ( (rho<0.0_c_double) ) THEN
                   set%user_rho = 0.0_c_double
               ELSE       
                   set%user_rho = rho 
               END IF
               
                                     
               ! ***** descent parameter 'm' *****
               IF ( (set%user_m<=0.0_c_double) .OR. (set%user_m>=1.0_c_double) ) THEN
                   m = 0.2_c_double
               ELSE       
                   m = set%user_m 
               END IF
              
               ! ***** stopping tolerance 'crit_tol' *****
               IF (set%user_crit_tol <= 0.0_c_double) THEN
               
                   IF (user_n < 50) THEN
                       crit_tol = (10.0_c_double)**(-5)
                   ELSE IF ((50 <= user_n) .AND. (user_n <=200) ) THEN
                        crit_tol = (10.0_c_double)**(-5)
                   ELSE IF (user_n > 200 ) THEN 
                        crit_tol = (10.0_c_double)**(-4)                 
                   END IF
                   
               ELSE 
                   crit_tol = set%user_crit_tol
               END IF
              
               ! ***** enlargement parameter 'eps' *****
               IF (set%user_eps_1 <= 0.0_c_double) THEN
                   eps = 5*(10.0_c_double)**(-5)
               ELSE
                   eps = set%user_eps_1
               END IF
              
               ! ***** decrease parameter 'r_dec' *****
               IF ((set%user_r_dec <= 0.0_c_double) .OR. (set%user_r_dec >= 1.0_c_double) ) THEN
                   
                   IF( user_n < 10) THEN 
                      r_dec = 0.75_c_double
                   ELSE IF ( user_n == 10) THEN 
                      r_dec = 0.66_c_double
                   ELSE IF ( user_n == 11) THEN 
                      r_dec = 0.68_c_double   
                   ELSE IF ( user_n == 12) THEN 
                      r_dec = 0.70_c_double
                   ELSE IF ( user_n == 13) THEN 
                      r_dec = 0.72_c_double
                   ELSE IF ( user_n == 14) THEN 
                      r_dec = 0.73_c_double
                   ELSE IF ( user_n == 15) THEN 
                      r_dec = 0.75_c_double
                   ELSE IF ( user_n == 16) THEN 
                      r_dec = 0.76_c_double
                   ELSE IF ( user_n == 17) THEN 
                      r_dec = 0.77_c_double
                   ELSE IF ( user_n == 18) THEN 
                      r_dec = 0.78_c_double
                   ELSE IF ( user_n == 19) THEN 
                      r_dec = 0.79_c_double                     
                   ELSE IF ( user_n == 20 .OR. user_n == 21 ) THEN 
                      r_dec = 0.80_c_double    
                   ELSE IF ( user_n == 22 ) THEN 
                      r_dec = 0.81_c_double    
                   ELSE IF ( user_n == 23 .OR. user_n == 24 ) THEN 
                      r_dec = 0.82_c_double    
                   ELSE IF ( user_n == 25 .OR. user_n == 26 ) THEN 
                      r_dec = 0.83_c_double           
                   ELSE IF ( user_n == 27 .OR. user_n == 28 ) THEN 
                      r_dec = 0.84_c_double                         
                   ELSE IF ( user_n == 29 .OR. user_n == 30 ) THEN 
                      r_dec = 0.85_c_double                   
                   ELSE IF ( user_n >= 31 .AND. user_n <= 33 ) THEN 
                      r_dec = 0.86_c_double       
                   ELSE IF ( user_n >= 34 .AND. user_n <= 36 ) THEN 
                      r_dec = 0.87_c_double  
                   ELSE IF ( user_n >= 37 .AND. user_n <= 40 ) THEN 
                      r_dec = 0.88_c_double       
                   ELSE IF ( user_n >= 41 .AND. user_n <= 44 ) THEN 
                      r_dec = 0.89_c_double               
                   ELSE IF ( user_n >= 45 .AND. user_n <= 50 ) THEN 
                      r_dec = 0.90_c_double                             
                   ELSE IF ( user_n >= 51 .AND. user_n <= 57 ) THEN 
                      r_dec = 0.91_c_double  
                   ELSE IF ( user_n >= 58 .AND. user_n <= 66 ) THEN 
                      r_dec = 0.92_c_double   
                   ELSE IF ( user_n >= 67 .AND. user_n <= 78 ) THEN 
                      r_dec = 0.93_c_double                 
                   ELSE IF ( user_n >= 79 .AND. user_n <= 94 ) THEN 
                      r_dec = 0.94_c_double                     
                   ELSE IF ( user_n >= 95 .AND. user_n <= 119 ) THEN 
                      r_dec = 0.95_c_double
                   ELSE IF ( user_n >= 120 .AND. user_n <= 161 ) THEN 
                      r_dec = 0.96_c_double
                   ELSE IF ( user_n >= 162 .AND. user_n <= 244 ) THEN 
                      r_dec = 0.97
                   ELSE IF ( user_n >= 245 .AND. user_n <= 299 ) THEN 
                      r_dec = 0.98
                   ELSE IF ( user_n >= 300 ) THEN 
                      r_dec = 0.99
                   END IF 
                   
               ELSE
                   r_dec = set%user_r_dec
               END IF

               ! ***** increase parameter 'r_inc' *****
               IF ( set%user_r_inc <= 1.0_c_double ) THEN
                   r_inc = (10.0_c_double)**(7) 
               ELSE
                   r_inc = set%user_r_inc
               END IF

               ! ***** bundle B1 size *****
               IF ( set%user_size_b1 <= 0 ) THEN
                   size_b1 = MIN(user_n+5,1000) 
               ELSE
                   size_b1 = set%user_size_b1
               END IF

               ! ***** bundle B2 size *****
               IF ( set%user_size_b2 <= 0 ) THEN
                   size_b2 = 3
               ELSE
                   size_b2 = set%user_size_b2
               END IF      

               ! ***** decrease parameter 'c' *****               
               IF ((set%user_c<=0.0_c_double) .OR. (set%user_c>1.0_c_double)) THEN
                   c = 0.1_c_double
               ELSE
                   c = set%user_c
               END IF              
               
               !--------------------------------------------------------------------------
               !                      CLARKE STATIONARY ALGORITHM
               !--------------------------------------------------------------------------
               ! ***** descent paramter 'm_clarke' ***** 
               IF ((set%user_m_clarke<=0.0_c_double) .OR. (set%user_m_clarke>=1.0_c_double)) THEN
                   m_clarke = 0.01_c_double
               ELSE
                   m_clarke = set%user_m_clarke
               END IF
               
               ! ***** step-lenght tolerance 'step_tol' ***** 
               IF (set%user_eps <= 0.0_c_double) THEN
                   IF (user_n < 50) THEN
                       step_tol = 0.000001_c_double
                   ELSE
                       step_tol = 0.00001_c_double               
                   END IF

               ELSE
                   step_tol = set%user_eps
               END IF 
               
               !----------------------------------------------------------------------------
               ! ***** maximum number of main iterations 'mit' *****
               IF ( mit <= 0 ) THEN
                   mit = 5000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds' in the 'main iteration' *****
               IF ( mrounds <= 0 ) THEN
                   mrounds = 5000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds_clarke' in the Clarke stationary algorithm *****
               IF ( mrounds_clarke <= 0 ) THEN
                   mrounds_clarke = 5000
               END IF                  
               
               ! ***** print option 'iprint' *****
               IF ( (iprint < -4 ) .OR. (iprint > 4) ) THEN   !executed if print value is wrong
                   iprint = 1
               END IF              
              !----------------------------------------------------------------------------
              
           !_______________________________________________________________________________
           !************************ STEP 0: PARAMETER INITIALIZATION *********************                
           
               x_current = x_0                                ! the current iteration point is the starting point x_0
               f2_current = f2(set, x_0, problem2, user_n)    ! the value of the DC component f_2 at x_0
               f_counter = 1                                  ! one function value evaluated for both f_1 and f_2
              
               CALL f1_sub(set, x_0, problem1, user_n, f1_current, grad1)   ! the subgradient and function value at x_0
               subgrad1_counter = 1                                         ! one subgradient calculated for f_1

               f_0 = f1_current - f2_current     ! the value of the objective function at the starting point x_0 
               
               grad2 = subgradient_f2(set, x_0, problem2, user_n)    ! the subgradient of f_2 at x_0
               subgrad2_counter = 1                                  ! one subgradient calculated for f_1
 
 
               ! The bundles B_1 and B_2 are initialized
               CALL init_bundle_b1(B1, size_b1, user_n)    
               CALL init_bundle_b2(B2, size_b2, user_n)
               
               ! The first bundle element is added into the bundle B_1 and B_2 (i.e. the one corresponding to the starting point)
               CALL add_first_element_b1(B1, grad1)
               CALL add_first_element_b2(B2, grad2)
               
               iter_counter = 0             ! the number of 'main iterations' executed so far is zero
               subprob_counter = 0          ! the number of 'subproblems' solved so far is also zero
               stop_cond_counter = 0        ! the number of times approximate stopping condition is tested is zero
               stop_alg = .FALSE.           ! we cannot stop the proximal bundle algorithm
               
               clarke_f_counter = 0         ! the initialization of function value counter for 'Clarke stationary' algorithm
               clarke_sub_counter = 0       ! the initialization of subgradient counter for 'Clarke stationary' algorithm
               
               delta = 0.0_c_double              

              
               IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 4) ) THEN        ! the basic/extended print of indermediate results     

                    WRITE(*,*) 'Main Iter:', 0, 'f(x):', f1_current - f2_current  

               END IF              
           !_______________________________________________________________________________
           !************************ STEP 0: END ******************************************    
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------  
           
            DO WHILE ( ( .NOT. stop_alg) .AND. (iter_counter < mit )  )  ! is repeated until the bundle algorithm DBDC can be terminated 
             
               !_______________________________________________________________________________
               !************************ STEP 1: MAIN ITERATION ******************************* 

                iter_counter = iter_counter + 1                     ! a new 'main iteration' is executed
                
                CALL cpu_time(start_time_main_it)
               
                ! The execution of main iteration
                CALL main_iteration( x_current, f1_current, f2_current, delta, user_n, &   
                            & x_new, f1_new, f2_new, f_0, B1, B2,  &
                            & crit_tol, eps, m, c, r_dec, r_inc, m_clarke, step_tol, mrounds,& 
                            & mrounds_clarke, reason_for_stop, help_mit_iter_counter, &
                            & help_subprob_counter, help_f_counter, help_subgrad1_counter, &
                            & help_subgrad2_counter, help_stop_cond_counter, &
                            & help_clarke_f_counter, help_clarke_sub_counter, &
                            & agg_used, stepsize_used, problem1, problem2, set)  


                CALL cpu_time(finish_time_main_it)
                
                ! Different counters are updated
                subprob_counter = subprob_counter + help_subprob_counter
                f_counter = f_counter + help_f_counter
                stop_cond_counter = stop_cond_counter + help_stop_cond_counter
                subgrad1_counter = subgrad1_counter + help_subgrad1_counter
                subgrad2_counter = subgrad2_counter + help_subgrad2_counter 
                
                clarke_f_counter = clarke_f_counter + help_clarke_f_counter
                clarke_sub_counter = clarke_sub_counter + help_clarke_sub_counter

               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               IF (reason_for_stop == 0) THEN   ! in this case, a new iteration point is found in the 'main iteration'
               !_______________________________________________________________________________
               !************************ STEP 2: BUNDLE UPDATE*******************************
                
                    ! Subgradients of f_1 and f_2 at the new point x_new
                    grad1 = subgradient_f1(set, x_new, problem1, user_n)
                    subgrad1_counter = subgrad1_counter + 1
                    change1 = f1_new - f1_current
                    
                    IF (problem2==2) THEN                                  ! L0-norm (parametric) 
                      apuf = f2(set, x_new, problem2, user_n)                   ! Calculated in order to have correct order in k-norm
                    END IF
                    
                    grad2 = subgradient_f2(set, x_new, problem2, user_n)
                    subgrad2_counter = subgrad2_counter + 1
                    change2 = f2_new - f2_current
                    
                    dd = x_new - x_current                              ! the search direction
                    change = f1_new - f2_new - (f1_current -f2_current) ! the change in the objective function value
                
                    CALL update_b1(B1, grad1, dd, change1)              ! bundle update for B_1
                    CALL update_b2(B2, grad2, dd, change2)              ! bundle update for B_2
                                        
                    x_current = x_new                                   ! update of the current iteration point
        
                    
                    f1_current = f1_new                                 ! update of the function value f_1
                    f2_current = f2_new                                 ! update of the function value f_2
                    
                    norm = SQRT(DOT_PRODUCT(dd,dd))                     ! Distance between two consecutive iteration points
                               
                    IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 4) ) THEN        ! the basic/extended print of indermediate results
                        
                        WRITE(*,*) 'Main Iter:', iter_counter, 'f(x):', f1_current-f2_current, '  change:', change, 'norm: ', norm
                        IF (iprint > 3) THEN    ! IF iprint =  4 then the iteration point is printes
                            WRITE(*,*) 'The new iteration point:'
                            DO i = 1, user_n
                                WRITE(*,*) 'x(', i ,')*=',  x_current(i)
                            END DO
                            WRITE(*,*) ' '
                        END IF 
                        
                        IF ( ABS(iprint) == 4 ) THEN
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                            WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * ' 
                            WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                            WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                            WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                            WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                            WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                            IF (help_stop_cond_counter > 0) THEN
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                                WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                            END IF
                            
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) ' '                          
                        END IF
                        

                        
                    END IF  
                            
               !_______________________________________________________________________________
               !************************ STEP 2: END ****************************************       
               
               ELSE                     ! one of the stopping conditions is fulfilled
               
                    stop_alg = .TRUE.   ! the minimization algorithm can be STOPPED
                    
                    IF(ABS(iprint) == 4) THEN ! the extended print of indermediate results
                        WRITE(*,*) ' '
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * '
                        WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                        WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                        WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                        WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                        WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                        IF (help_stop_cond_counter > 0) THEN    
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                            WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                        END IF                          
                            
                        WRITE(*,*) '--------------------------------------------------------------------------- '       
                        WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        WRITE(*,*) ' '
                        WRITE(*,*) 'During Main iteration round:', iter_counter
                        WRITE(*,*) 'Some of the stopping conditions is fulfilled and the algorithm is stopped.'                     
                        WRITE(*,*) ' '
                    END IF 
                    
               END IF
               
           END DO
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------         
           

           termination = reason_for_stop   ! the cause of the termination
           
           IF ( iter_counter >= mit ) THEN  ! the maximum number of 'main iterations' is executed and this causes the termination
               termination = 4 
           END IF
           
           x_solution = x_current               ! the solution to the minimization problem
           f_solution = f1_current - f2_current ! the objective function value at the solution
           
           ! the values of different counters
           counter(1) = iter_counter 
           counter(2) = subprob_counter
           counter(3) = f_counter
           counter(4) = subgrad1_counter
           counter(5) = subgrad2_counter 
           counter(6) = stop_cond_counter
           counter(7) = clarke_f_counter
           counter(8) = clarke_sub_counter
           
            
           CALL cpu_time(finish_time)         ! Stop CPU timing
           CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
          ! Calculate the elapsed 'clock' time in seconds:
          elapsed_time=(1.0_c_double*clock_end-clock_start)/clock_rate       
          
          CPUtime = finish_time-start_time        

           IF ( (ABS(iprint) == 1) ) THEN ! basic print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '           * * * * * BASIC PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '               
               WRITE(*,*) 'f*=', f_solution     
               WRITE(*,*) '---------------------------------------------------------------------------' 
               WRITE(*,*) 'CPU time used:', finish_time-start_time                              
               WRITE(*,*) '---------------------------------------------------------------------------'
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'            
           END IF
           
           IF ((ABS(iprint) /= 1) .AND. (iprint /= 0)) THEN  ! extended print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '        * * * * * EXTENDED PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '               
               WRITE(*,*) 'f*=', f_solution                    
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '   
               WRITE(*,*) 'The number of main iterations:', iter_counter
               WRITE(*,*) 'The number of subproblems solved:', subprob_counter             
               WRITE(*,*) 'The number of function values evaluated for DC component f_1:', f_counter 
               WRITE(*,*) 'The number of function values evaluated for DC component f_2:', f_counter 
               WRITE(*,*) 'The number of subgradients calculated for DC component f_1:', subgrad1_counter 
               WRITE(*,*) 'The number of subgradients calculated for DC component f_2:', subgrad2_counter 
               WRITE(*,*) '---------------------------------------------------------------------------'
               WRITE(*,*) 'The number of times Clarke stationary algorithm (CSA) was used:', stop_cond_counter  
               WRITE(*,*) 'The number of function values computed for f in CSA:', clarke_f_counter
               WRITE(*,*) 'The number of subgradients calculated for f in CSA: ', clarke_sub_counter
               WRITE(*,*) '---------------------------------------------------------------------------'     
               WRITE(*,*) 'CPU time used:', finish_time-start_time                          
               WRITE(*,*) '--------------------------------------------------------------------------- '
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'                
           END IF 
           


           CALL deallocation_b1(B1)        
           CALL deallocation_b2(B2)        

           
           END SUBROUTINE DBDC_algorithm      
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|


        
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             THE MAIN ITERATION                                 |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE main_iteration( x_k, f1_k, f2_k, f_val, user_n, &
                            & x_new , f1_new, f2_new, f_0, B1, B2, &
                            & crit_tol, eps, m, c, r_dec, r_inc, m_clarke, step_tol, mrounds, &
                            & mrounds_clarke, reason_for_stop, iter_counter, &
                            & subprob_counter, fi_counter, subgrad1_counter, subgrad2_counter,&
                            & stop_cond_counter, clarke_f_counter, clarke_sub_counter, &
                            & agg_in_use, stepsize_used, problem1, problem2, set)                    
            !
            ! Executes the 'main iteration' algorithm. It is needed to find the search direction at the current iteration point. 
            !
            ! INPUT: * 'x_k'               : The current iteration point
            !        * 'user_n'            : The dimension of the problem
            !        * 'problem1'          : The f1 of the problem
            !        * 'problem2'          : The f2 of the problem
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !        * 'f_val'             : The value of objective function used in descent condition
            !        * 'f_0'               : The value of the objective function f at a starting point x_0
            !
            !        * 'crit_tol'          : The stopping tolerance 
            !        * 'eps'               : The enlargement parameter
            !        * 'm'                 : The descent parameter used in main iteration
            !        * 'c'                 : The decrease parameter
            !        * 'r_dec' and 'r_inc' : The decrease and increase parameters   
            !        * 'm_clarke'          : The descent parameter used in 'Clarke stationary' algorithm 
            !        * 'step_tol'          : The step-length tolerance used in 'Clarke stationary' algorithm (i.e. proximity measure)
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'main iteration'
            !        * 'mrounds_clarke'    : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !
            !        * 'agg_in_use'        : If .TRUE. then aggregation is used in the algorithm
            !        * 'stepsize_used'     : If .TRUE. then simple stepsize determination is used in the algorithm after each main iteration
            !
            !
            ! OUTPUT: * 'x_new'              : the new iteration point
            !         * 'f1_new' and 'f2_new': the new value of DC components f_1 and f_2 at 'x_new'
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'main iteration' 
            !
            !         * 'iter_counter'       : the number of rounds needed in 'main iteration'
            !         * 'subprob_counter'    : the number of subproblems solved in 'main iteration'
            !         * 'fi_counter'         : the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2) (wihtout the ones used in Clarke stationary algorithm)
            !         * 'subgrad1_counter'   : the number of subgradients calculated for f_1 in 'main iteration' (without the ones used in Clarke stationary algorithm)
            !         * 'subgrad2_counter'   : the number of subgradients calculated for f_2 in 'main iteration' (without the ones used in Clarke stationary algorithm)
            !         * 'stop_cond_counter'  : the number of times approximate stopping condition is tested during the 'main iteration' 
            !         * 'clarke_f_counter'   : the number of function values evaluated for f in 'Clarke stationary' algorithm
            !         * 'clarke_sub_counter' : the number of subgradients evaluated for f in 'Clarke stationary' algorithm
            !
            ! INOUT: * 'B_1' and B_2'        : the bundles of the DC components f_1 and f_2
            !        * 'set'                 : The set of information
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE main_iteration is used in SUBROUTINE bundle_method.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER *********************************  
               
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information               
            
               TYPE(kimppu1), INTENT(INOUT) :: B1                   ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                   ! the bundle B_2 for the DC component f_2
               REAL(KIND=c_double), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'main iteration'

               REAL(KIND=c_double), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=c_double), INTENT(IN)  :: f_val           ! the value used in descent condition
               REAL(KIND=c_double), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'main iteration'
               REAL(KIND=c_double), INTENT(IN)  :: f_0             ! the value of f at the starting point x_0
               
               REAL(KIND=c_double), INTENT(IN) :: crit_tol       ! crit_tol=stopping tolerance 
               REAL(KIND=c_double), INTENT(IN) :: eps            ! eps=enlargemant parameter
               REAL(KIND=c_double), INTENT(IN) :: m              ! m=descent parameter
               REAL(KIND=c_double), INTENT(IN) :: c              ! c=decrease parameter 
               REAL(KIND=c_double), INTENT(IN) :: r_dec, r_inc   ! r_dec=decrease parameter and R_inc=increase parameter
               REAL(KIND=c_double), INTENT(IN) :: m_clarke       ! The descent parameter used in 'Clarke stationary' algorithm 
               REAL(KIND=c_double), INTENT(IN) :: step_tol       ! The step-length tolerance used in 'Clarke stationary' algorithm (i.e. proximity measure)
               
               INTEGER(KIND=c_int) :: user_n                           ! the dimension of the problem
               INTEGER(KIND=c_int) :: problem1                         ! the f1 of the problem
               INTEGER(KIND=c_int) :: problem2                         ! the f2 of the problem
               
               INTEGER(KIND=c_int) :: mrounds                          ! the maximum number of possible rounds in the 'main iteration'
                                                           ! If mrounds<=0, then DEFAULT value 500 is used (This is also done in the SUBROUTINE bundle_method()
               INTEGER(KIND=c_int) :: mrounds_clarke                   ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm

               INTEGER(KIND=c_int), INTENT(OUT) :: reason_for_stop     ! 0 - a new iteration point found
                                                           ! 1 - stopping condition satisfied (i.e. Clarke stationarity)
                                                           ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)
                                                           ! 3 - maximum number of rounds executed in 'main iteration'
                                                           ! 5 - maximum number of rounds executed in 'Clarke stationary' algorithm
               
               INTEGER(KIND=c_int), INTENT(OUT) :: iter_counter        ! the number of rounds used in 'main iteration'
               INTEGER(KIND=c_int), INTENT(OUT) :: subprob_counter     ! the number of subproblems solved in 'main iteration'
               INTEGER(KIND=c_int), INTENT(OUT) :: fi_counter          ! the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2) (without the ones used in Clarke stationary algorithm)
               INTEGER(KIND=c_int), INTENT(OUT) :: subgrad1_counter    ! the number of subgradients calculated for f_1 in 'main iteration' (without the ones used in Clarke stationary algorithm)
               INTEGER(KIND=c_int), INTENT(OUT) :: subgrad2_counter    ! the number of subgradients calculated for f_2 in 'main iteration' (without the ones used in Clarke stationary algorithm)
               INTEGER(KIND=c_int), INTENT(OUT) :: stop_cond_counter   ! the number of times Clarke stationary algortihm is entered during the 'main iteration'  

               INTEGER(KIND=c_int), INTENT(OUT) :: clarke_f_counter    ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               INTEGER(KIND=c_int), INTENT(OUT) :: clarke_sub_counter  ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               LOGICAL :: agg_in_use              ! If .TRUE. then aggregation is used in the algorithm 
               LOGICAL :: stepsize_used           ! If .TRUE. then step-size determination is used in the algorithm 
               

           !***************************** LOCAL VARIABLES ************************************
               
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: d_t                    ! the search direction
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: y                      ! the new auxilary point
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at y
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: agg_grad               ! the new aggregated subgradient of f_1
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: apu_grad               ! 'help' subgradient
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: vect                   ! 'help' vector               
               
               REAL(KIND=c_double) :: norm1              ! ||\bxi_1(x_k)||
               REAL(KIND=c_double) :: max_norm           ! the value of the maximum subgradient norm ||\bxi_{2,max}||
               REAL(KIND=c_double) :: t                  ! the proximity parameter t
               REAL(KIND=c_double) :: div_t              ! 1.0_c_double divided by the proximity parameter t
               REAL(KIND=c_double) :: t_min, t_max       ! the bounds for proximity parameter t
               
               REAL(KIND=c_double) :: d_norm             ! the norm of the direction vector d_t
               REAL(KIND=c_double) :: delta1, delta2     ! the predicted changes of f_1 and f_2, respectively
               REAL(KIND=c_double) :: f1_y, f2_y         ! the function values of f_1 and f_2 at y
               REAL(KIND=c_double) :: f_k                ! the function value of f at x_k
               REAL(KIND=c_double) :: real_decrease      ! the real decrease of the objective function f
               REAL(KIND=c_double) :: lin_err1, lin_err2 ! the linearization errors of f_1 and f_2 (calculated using y)
               REAL(KIND=c_double) :: agg_lin_err        ! the new aggregated linearization error of f_1
               
               REAL(KIND=c_double) :: help               ! 'help' variable  

               INTEGER(KIND=c_int) :: f_counter                ! help counter for the number of function values
               INTEGER(KIND=c_int) :: sub_counter              ! help counter for the number of subgradients 
               INTEGER(KIND=c_int) :: clarke_iter_counter      ! help counter for the number  of iterations in Clarke stationary algorithm 
               
               INTEGER(KIND=c_int) :: s_counter                ! the number of subproblems solved during the current iteration round of the 'main iteration' algorithm
               INTEGER(KIND=c_int) :: i, ind                   ! help variables
                           
               LOGICAL :: t_changed                ! .TRUE. if t has been changed during the previous round (.TRUE. also when initialization is done)
               LOGICAL :: was_b1_full              ! .TRUE. if during the previous round something was added to B_1 and B_1 was full before this insertion
               LOGICAL :: agg_used                 ! .TRUE. if aggregation element is in use
               LOGICAL :: stop_main_it             ! .TRUE. if the current 'main iteration' can be stopped
               
               LOGICAL :: stop_linesearch          ! .TRUE. if the current 'line_search' can be stopped
               REAL(KIND=c_double) :: stepsize           ! stepsize determined during 'line search' into descent direction
               REAL(KIND=c_double) :: test_f1            ! the value of f1 at the point tested in 'line search'
               REAL(KIND=c_double) :: test_f2            ! the value of f2 at the point tested in 'line search'
               REAL(KIND=c_double), DIMENSION(give_n_b1(B1)) :: test_y  ! a point tested during 'line search'
                       

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           
           !__________________________________________________________________________________     
           !************************ STEP 0: STOPPING CONDITION ******************************
                
               iter_counter = 1                   ! the first round is executed
               subprob_counter = 0                ! the number of subproblems solved
               fi_counter = 0                     ! the number of function values evaluated
               subgrad1_counter = 0               ! the number of subgradients calculated for f_1
               subgrad2_counter = 0               ! the number of subgradients calculated for f_2
               stop_cond_counter = 0              ! the number of times approximate stopping condition is tested
               clarke_f_counter = 0               ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               clarke_sub_counter = 0             ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
                   
               agg_used = .FALSE.                 ! Aggregation is initialized to be .FALSE.
           
               stop_main_it = .FALSE.             ! We are NOT ready to stop
               
               vect = give_subgrad_b1(B1,0) - give_subgrad_b2(B2,0)       ! vector: ' \bxi_1(x_k) - \bxi_2(x_k) '
               
               start: IF ( SQRT(DOT_PRODUCT(vect,vect)) < (crit_tol)) THEN  ! ||\bxi_1(x_k) - \bxi_2(x_k)|| <  crit_tol 
               !**>>**>>**>> APPROXIMATE CRITICALITY OBTAINED <<**<<**<<**
                   
                   d_t = 1.0_c_double
                   
                   ! Clarke stationary algorithm is executed
                   CALL guaranteeing_clarke( x_k, f1_k, f2_k, reason_for_stop, &
                            & x_new , f1_new, f2_new, d_t, &
                            & crit_tol, m_clarke, step_tol, mrounds_clarke, &
                            & clarke_iter_counter, f_counter, sub_counter, &
                            & problem1, problem2, user_n, set)
                    
                   ! Counters are updated   
                   clarke_f_counter = clarke_f_counter + f_counter                  
                   clarke_sub_counter = clarke_sub_counter + sub_counter                    
           
                   stop_cond_counter = stop_cond_counter + 1    ! the update of the stopping condition counter
                   stop_main_it = .TRUE.                        ! the 'main iteration' can be stopped
                  
           !__________________________________________________________________________________         
           !************************ STEP 0: END *********************************************         
               ELSE start   
               !_______________________________________________________________________________
               !************************ STEP 1: PARAMETER INITIALIZATION *********************
           
                   f_k = f1_k - f2_k                       ! the function value of f at x_k
                   
                   vect = give_subgrad_b1(B1,0)            ! the vector \bxi_1(x_k)
                   norm1 = SQRT(DOT_PRODUCT(vect,vect))    ! the value of the norm ||\bxi_1(x_k)||
                   max_norm = max_norm_value(B2)           ! the value of the maximum subgradient norm in the bundle B_2
                                    
                   t_min = (0.5_c_double * r_dec * eps) / ( norm1 + max_norm )  ! the lower bound for t
                   t_max = r_inc * t_min                                  ! the upper bound for t
                   
                   t = select_value_t(t_min,t_max)    ! the parameter t is selected         
            
                   t_changed = .TRUE.                 ! the parameter t was changed
                   was_b1_full = .FALSE.              ! B_1 was not full during the previous round since this is the first round of 'main iteration'
                       
                                   
                   IF ( mrounds <= 0 ) THEN           
                       mrounds = 5000                 ! the DEFAULT value for 'mrounds' is used
                   END IF
               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               END IF start

               
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------
           
               DO WHILE ( (.NOT. stop_main_it ) .AND. (iter_counter <= mrounds ))   ! is repeated until the 'main iteration' can be terminated         
               !______________________________________________________________________________
               !************************ STEP 2: SEARCH DIRECTION ****************************
                   
   
                   IF (agg_in_use) THEN             ! 'agg_in_use' tells whether we use aggregation or not. If 'agg_in_use=.FALSE. then also 'agg_used'=.FALSE.
                       agg_used = is_agg_used(B1)   ! tells whether aggregation element is used or not? 
                   END IF   
                   
                   IF ( agg_used ) THEN             ! IF we have .FALSE. here then the algorithm does NOT use AGGREGATION

                        CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+2, &  ! subproblems are solved with aggregation
                                          &  B1, B2, t, s_counter)
                   
                   ELSE
                   
                        CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+1, &  ! subproblems are solved without aggregation
                                          &  B1, B2, t, s_counter)
                   END IF   
                                        
                   subprob_counter = subprob_counter + s_counter    ! the number of subproblems solved so far
                                   
                   CALL add_glob_index(B2)              ! calculates the index of the subproblem yielding the global solution        
                   d_t = give_solution(B2)              ! the global solution d_t is selected   
                   d_norm = SQRT(DOT_PRODUCT(d_t,d_t))  ! the norm ||d_t|| of the solution d_t is calculated
                   
                   !WRITE(*,*) 'd_norm', d_norm

                   ! Aggregated element is calculated. 
                   ! NOTICE: The aggregated element (if used) is now always added into bundle B_1 (also in those cases when the bundle B_1 is not full)

                   ! It is also possible to add the aggregated element into bundle B_1 only when the bundle B_1 is full
                   ! IF (is_full_b1(B1)) THEN

                   IF (agg_in_use) THEN                     ! Is done only when 'agg_in_use'=.TRUE.
                       ind = give_solution_ind(B2)
                       apu_grad = give_subgrad_b2(B2,ind)
                       div_t = 1.0_c_double / t
                       DO i = 1, give_n_b2(B2)
                          agg_grad(i) = (-div_t) * d_t(i) + apu_grad(i)
                       END DO
                       agg_lin_err = - give_decrease(B2) - (div_t) * DOT_PRODUCT(d_t,d_t)
                       agg_lin_err = agg_lin_err + give_linerr_b2(B2,ind)
                       CALL add_agg_element_b1(B1,agg_grad,agg_lin_err)
                   END IF

                   !END IF  
                                    
               !______________________________________________________________________________
               !************************ STEP 2: END *****************************************     
                   
                   
               !->->->->->-> EXECUTION OF BRANCH BEGINS (2 POSSIBLE BRANCHES) <-<-<-<-<-<-<-<-
                   branches: IF (d_norm < crit_tol) THEN                
               !->->->->->->->->->->->->->-> BRANCH 1 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-             
                   !__________________________________________________________________________
                   !******************** STEP 3: CLARKE STATIONARITY CHECK *******************
                                           
                        ! Clarke stationary algorithm is executed
                        CALL guaranteeing_clarke( x_k, f1_k, f2_k, reason_for_stop, &
                            & x_new , f1_new, f2_new, d_t, &
                            & crit_tol, m_clarke, step_tol, mrounds_clarke, &
                            & clarke_iter_counter, f_counter, sub_counter, &
                            & problem1, problem2, user_n, set)
                       
                        !Counters are updated
                        clarke_f_counter = clarke_f_counter + f_counter                 
                        clarke_sub_counter = clarke_sub_counter + sub_counter   
                       
                        stop_cond_counter = stop_cond_counter + 1   ! the update of the stopping condition counter
                        
                        stop_main_it = .TRUE.                       ! Main iteration can be stopped

                   !__________________________________________________________________________
                   !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************  
                
               !->->->->->->->->->->->->->-> BRANCH 1 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                        
                   ELSE branches
               !->->->->->->->->->->->->->-> BRANCH 2 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                   !__________________________________________________________________________
                   !******************** STEP 4: DESCENT TEST ********************************
 
                       IF (SQRT(DOT_PRODUCT(d_t,d_t))>1.0_c_double) THEN 
                          d_norm = SQRT(DOT_PRODUCT(d_t,d_t))
                          d_t = d_t / d_norm
                           
                       END IF     
                          
 
                       y = x_k + d_t                      ! a new auxilary point y
                       f1_y = f1(set, y, problem1,user_n)      ! the value of f_1 at y
                       f2_y = f2(set, y, problem2,user_n)      ! the value of f_2 at y             
                       
                       fi_counter = fi_counter + 1        ! one more objective function value calculated for a DC component
               
                       real_decrease = ( f1_y - f2_y ) -  f_k  ! the real decrease in the value of f    

                       
                       IF (give_decrease(B2)>0.0_c_double) THEN 
                         WRITE(*,*) "Problems! Positive value of predicted decrease!", give_decrease(B2)
                         real_decrease = 100000_c_double
                       END IF
                   
                       branch2: IF ( real_decrease <= (m * give_decrease(B2) + f_val ) ) THEN   ! the descent condition is dependent on the value 'f_val' 
                       !**>>**>>**>> NEW ITERATION POINT FOUND <<**<<**<<**  
                           
                           IF(stepsize_used) THEN 
                           !-**--**-  STEPSIZE DETERMINATION  -**--**-
                           IF (real_decrease < -0.1_c_double) THEN
                              stepsize = 1.0_c_double
                              stop_linesearch = .FALSE.
                              DO WHILE((.NOT. stop_linesearch))
                                 test_y = y + d_t
                                 test_f1 = f1(set, test_y, problem1,user_n)     ! the value of f_1 at test_y
                                 test_f2 = f2(set, test_y, problem2,user_n)     ! the value of f_2 at test_y   
                                 fi_counter = fi_counter + 1    
                                 IF ( (test_f1 - test_f2 - f_k) <= ( 0.7_c_double * m * give_decrease(B2))) THEN 
                                    stepsize = stepsize + 1.0_c_double
                                    y = test_y
                                    f1_y = test_f1
                                    f2_y = test_f2                                   
                                 ELSE
                                    stop_linesearch = .TRUE.  
                                 END IF
                              END DO
                           END IF 
                           !-**--**- STEPSIZE DETERMINATION END -**--**-
                           END IF
                           
                           x_new = y               ! the new iteration point is the current auxilary point
                           f1_new = f1_y           ! the value of f_1 at x_new
                           f2_new = f2_y           ! the value of f_2 at x_new
                           stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                           reason_for_stop = 0     ! the reason for stopping is the new iteration point
    
                           
                   !__________________________________________________________________________
                   !******************** STEP 4: END *****************************************         
                   
                       ELSE branch2 
                       !----------- SOME UPDATES IN VARIABLES BEGINS -------------------------  
                       
                   !__________________________________________________________________________     
                   !******************** STEP 5: BUNDLE UPDATE *******************************
                          
                          update: IF ( ((f1_y - f2_y - f_0)>0 ) .AND. & 
                                                      & (d_norm > eps)) THEN 
                          !----------- PARAMETER t REDUCTION --------------------------------
                              t = t - r_dec * ( t - t_min )  ! the parameter t is updated
                              t_changed = .TRUE.             ! the parameter t was changed
 

                          ELSE update
                          !----------- BUNDLE AND PARAMETER UPDATE --------------------------
                              
                              i = give_solution_ind(B2)              ! the index of the subproblem yielding the global solution
                              vect = give_subgrad_b2(B2, i)   ! the subgradient of f_2 in the subproblem yielding the global solution
                              delta2 = - DOT_PRODUCT(vect, d_t) + give_linerr_b2(B2,i)  ! the value of delta_2
                              delta1 = give_decrease(B2) - delta2                       ! the value of delta_1
                              
                              new_grad1 = subgradient_f1(set, y, problem1,user_n)              ! a subgradient of f_1 at y
                              subgrad1_counter = subgrad1_counter + 1              ! a new subgradient was evaluated for f_1            

                              lin_err1 = f1_k - f1_y + DOT_PRODUCT(d_t,new_grad1)  ! a new linearization error for f_1

                              IF (real_decrease > -m * give_decrease(B2))  THEN    ! adjustment of the proximity parameter if necessary
                                  t =  t - c * ( t - t_min )
                              END IF                          
                              
                            !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                              bundle1: IF (is_full_b1(B1)) THEN   
                              ! bundle B_1 is full and a new element is added to the bundle B_1
                              ! the overwritten element is written/taken down 
                                 
                                  CALL add_element_b1(B1, new_grad1, lin_err1)
                                  was_b1_full = .TRUE.
                                
                              ELSE bundle1
                              ! bundle B_1 is NOT full and a new element is added to the bundle B_1

                                  CALL add_element_b1(B1, new_grad1, lin_err1)
                                  was_b1_full = .FALSE.

                              END IF bundle1                             
                            !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-


                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE BEGINS <<<<<<<<<<<-             

                            bundle2: IF( delta2 >= 0 .AND. (give_max_size_b2(B2)>0) ) THEN   ! bundle B_2 is updated (happens only when delta2 >= 0 and bundle size is larger than 1) (SOME OTHER INSERTION RULES ALSO POSSIBLE TO APPLY!)

 
                                  new_grad2 = subgradient_f2(set, y, problem2,user_n)       ! a subgradient of f_2 at y  
                                  subgrad2_counter = subgrad2_counter + 1              ! a new subgradient was evaluated for f_2
                                  lin_err2 = f2_k - f2_y + DOT_PRODUCT(d_t,new_grad2)  ! a new linearization error for f_2  
                                
                                  CALL add_element_b2(B2, new_grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                                  !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.
                                  
                   !__________________________________________________________________________     
                   !******************** STEP 5: END *****************************************
                   
                   !__________________________________________________________________________     
                   !******************** STEP 6: PARAMETER UPDATE ****************************
                   
                                  help = SQRT(DOT_PRODUCT(new_grad2,new_grad2)) ! the norm of the new subgradient of f_2 (subgradient is calculated 
                                                                                ! at the new auxilary point y)
                                  IF(help > max_norm) THEN 
                                      max_norm = help            ! the updated value of the maximum subgradient norm in the bundle B_2
                                      t_min = (0.5_c_double * r_dec * eps) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                  END IF            
     
                   !__________________________________________________________________________     
                   !******************** STEP 6: END *****************************************
                   
                            END IF bundle2

                            
                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE ENDS <<<<<<<<<<<<<-
                               
                              t_changed = .FALSE.   ! the parameter t was NOT changed during this ELSE branch of 'update'
                                  
                          END IF update        
                        !------ PARAMETER t REDUCTION & BUNDLE AND PARAMETER UPDATE ENDS -----
                        
                          iter_counter = iter_counter + 1  ! update of the iteration counter                        
                              
                       END IF branch2   
                     !--------- NEW ITERATION POINT & SOME UPDATES IN VARIABLES ENDS --------- 
                     
               !->->->->->->->->->->->->->-> BRANCH 2 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-   
                   END IF branches
               !->->->->->->-> EXECUTION OF BRANCH ENDS  (2 BRANCHES) <-<-<-<-<-<-<-<-<-<-<-<-
              END DO 
              
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------           
              
              IF ( (.NOT. stop_main_it ) ) THEN 
                  reason_for_stop = 3            ! the maximum number of rounds have been executed 
                  iter_counter = iter_counter -1 
              END IF        
              
           END SUBROUTINE main_iteration
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        

        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>       
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           PARAMETER t SELECTION                                |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           
           
           FUNCTION select_value_t(t_min,t_max) RESULT(t)
                IMPLICIT NONE
                REAL(KIND=c_double), INTENT(IN) :: t_min, t_max
                REAL(KIND=c_double) :: t
            
                t = 0.8_c_double * (t_min + t_max) ! selection of t
                
                !OTHER OPTIONS:
                
!               t = 0.5_c_double * (t_min + t_max)    ! selects the middle point from interval [t_min, t_max]
!               t = t_min + 0.8*(t_max - t_min)  

!               t = t_min                   ! selects the lower bound
!               t = t_max                   ! selects the upper bound
                
           END FUNCTION select_value_t
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|    


  
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                      Guaranteeing Clarke stationarity                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE guaranteeing_clarke( x_k, f1_k, f2_k, reason_for_stop,&
                            & x_new , f1_new, f2_new, d_cause,&
                            & crit_tol, m, step_tol, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, &
                            & problem1, problem2, user_n, set)
               
            !
            ! Executes the 'Clarke stationary' algorithm. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary, then this algorithm generates a descent direction.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !
            !        * 'd_cause'           : The search direction causing the execution of the Clarke stationary algorithm
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance (i.e. proximity measure)
            !
            !        * 'user_n'            : The dimension of the problem
            !        * 'problem1'          : The f1 of the problem
            !        * 'problem2'          : The f2 of the problem
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !
            !
            ! OUTPUT: * 'x_new'              : the new iteration point (obtained if 'reason_for_stop' = 0)
            !         * 'f1_new' and 'f2_new': the new values of DC components f_1 and f_2 (obtained if 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'Clarke stationarity' algorithm  
            !
            !         * 'iter_counter'       : the number of rounds needed in 'Clarke stationary' algorithm
            !         * 'f_counter'          : the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !         * 'subgrad_counter'    : the number of subgradients calculated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !
            !
            ! INOUT   * 'set'                : The set of information
            !
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE guaranteeing_clarke is used in SUBROUTINE main_iteration.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
        
               TYPE(set_info), INTENT(INOUT) :: set          ! The set of information
               
               REAL(KIND=c_double), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=c_double), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Clarke stationary' algorithm

               REAL(KIND=c_double), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=c_double), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               
               REAL(KIND=c_double), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance 
               REAL(KIND=c_double), INTENT(IN) :: m             ! 'm'=descent parameter
               REAL(KIND=c_double), INTENT(IN) :: step_tol      ! 'step_tol'=step-length tolerance (i.e. proximity measure)
               
               INTEGER(KIND=c_int), INTENT(IN):: user_n                 ! the dimension of the problem
               INTEGER(KIND=c_int), INTENT(IN):: problem1               ! the f1 of the problem
               INTEGER(KIND=c_int), INTENT(IN):: problem2               ! the f2 of the problem
               
               INTEGER(KIND=c_int) :: mrounds                         ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm
                                                          
               INTEGER(KIND=c_int), INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)
                                                          ! 5 - the maximum number of rounds executed in 'Clarke stationary' algorithm
                                                          
               INTEGER(KIND=c_int), INTENT(OUT) :: iter_counter       ! the number of rounds used in 'Clarke stationary' algorithm
               INTEGER(KIND=c_int), INTENT(OUT) :: f_counter          ! the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
               INTEGER(KIND=c_int), INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for f in 'Clarke stationary' algorithm 
               

           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                ! the bundle B of f containing subgradients from '\partial f(x)'
               
               REAL(KIND=c_double), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=c_double), DIMENSION(user_n) :: x_d                    ! an auxiliary point
               REAL(KIND=c_double), DIMENSION(user_n) :: new_grad               ! the subgradient of f at x
               REAL(KIND=c_double), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at x_d
               REAL(KIND=c_double), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=c_double), DIMENSION(user_n) :: y                      ! a new auxiliary point            
               REAL(KIND=c_double), DIMENSION(user_n) :: test_y                 ! another auxiliary point            
               
               REAL(KIND=c_double) :: f_k              ! the value of the objective function at the current iteration point   
               REAL(KIND=c_double) :: apuf             ! the help objective value 
               REAL(KIND=c_double) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=c_double) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=c_double) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=c_double) :: div_eps          ! 1.0_c_double / eps     
               REAL(KIND=c_double) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=c_double) :: f1_y, f2_y       ! the function values of f_1 and f_2 at y
               REAL(KIND=c_double) :: test_f1, test_f2 ! the function values of f_1 and f_2 at test_y
               REAL(KIND=c_double) :: direc_der        ! directional derivative of f 
               REAL(KIND=c_double) :: stepsize         ! stepsize 
               REAL(KIND=c_double) :: descent_app      ! approximation of descent in the objective function value
                            
               INTEGER(KIND=c_int) :: size_b                 ! the biggest possible bundle size for B
               INTEGER(KIND=c_int) :: N                      ! the current bundle size for B with the current element
                            
               LOGICAL :: stop_alg               ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_step_det          ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************              
                        
            ! the size of the bundle B
               IF (set%user_size < 2) THEN
                   !size_b = user_n + 3 
                   size_b = user_n *2 
               ELSE
                   size_b = set%user_size
               END IF
               
               eps = (10.0_c_double)**(-4)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_c_double / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               
               f_k = f1_k - f2_k                      ! The value of f at the current iteration point
               
               d = -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))    ! The search direction used to determine subgradient for f         
               
            ! The bundles B is initialized
               CALL init_bundle_b1(B, size_b, user_n) 
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.          ! Algorithm cannot be stopped
               

           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************         

               
               IF (iter_counter > 100) THEN  ! If over 100 iterations are executed then stepsize 'eps' is modified
                   eps = (10.0_c_double)**(-6)
               END IF
           
               x_d = x_k + (eps) * d                      ! Auxialiary point 'x_d' is calculated
               
               new_grad1 = subgradient_f1(set, x_d, problem1, user_n)  ! Subgradient of f_1 at 'x_d'
               
               IF (problem2==2) THEN                                  ! L0-norm (parametric) 
                 apuf = f2(set, x_new, problem2, user_n)                   ! Calculated in order to have correct order in k-norm
               END IF
               new_grad2 = subgradient_f2(set, x_d, problem2, user_n)  ! Subgradient of f_2 at 'x_d'
               
               ! new subgradient for f at x_k 
               new_grad = new_grad1 - new_grad2
               subgrad_counter = subgrad_counter +1 
               
               ! New subgradien is added into the bundle 'B'
               IF (N == 0) THEN 
                  CALL add_first_element_b1(B, new_grad)    
                  N = N + 1
               ELSE
                  CALL add_element_b1(B, new_grad, 0.0_c_double)
               END IF 
           
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CALRKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
               
               ! Minimum norm element 'u_k' is determined in the bundle 'B' (i.e. aggregated subgradient)
               CALL quadratic_solver(x_k, user_n, N, B, 1.0_c_double, u_k, obj_u) 
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
    
               norm_u = SQRT( DOT_PRODUCT(u_k,u_k))
               
            
               IF (norm_u < crit_tol) THEN  
                  
                  reason_for_stop = 1     ! Approximate Clarke stationarity is achieved
                  stop_alg = .TRUE.       ! Algorihtm can be stopped
                
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SERACH DIRECTION *****************************            
               
               d =  u_k / norm_u          ! New search direction
               y = x_k + eps * d          ! New auxialiary point    
               
               f1_y = f1(set,y,problem1, user_n)      ! f_1 at 'y'
               f2_y = f2(set, y,problem2, user_n)      ! f_2 at 'y'
               f_counter = f_counter +1 
               
               real_decrease = f1_y - f2_y - f_k      ! Real decraese in the objective function f

               direc_der = div_eps * real_decrease    ! Directional derivative at 'x_k' into the direction 'd'       
               
               IF (direc_der <= (-m * norm_u)) THEN   ! If .TRUE. decrease in the objective f is sufficient
                     reason_for_stop = 10             ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.                ! algorithm can be stopped  
               END IF
               
               iter_counter = iter_counter +1         ! one more iteration is executed
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1) .AND. stop_alg ) THEN 
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  stepsize = 1.0_c_double                      ! initial stepsize
                  stop_step_det = .FALSE.                ! stepsize determination cannot be stopped
                  DO WHILE((.NOT. stop_step_det))
                      test_y = x_k + stepsize * d        ! auxiliary test point
                      test_f1 = f1(set,test_y, problem1, user_n)     ! the value of f_1 at test_y
                      test_f2 = f2(set,test_y, problem2, user_n)     ! the value of f_2 at test_y   
                      f_counter = f_counter + 1          
                      descent_app = -(10.0_c_double)**(-4)     ! sufficient descent 
                      IF ( (test_f1 - test_f2 - f_k) > MIN(-m * stepsize * norm_u,descent_app )) THEN 
                         stepsize = 0.5_c_double * stepsize    ! stepsize is decreased
                         IF (stepsize < step_tol) THEN  
                               stop_step_det = .TRUE.    
                         END IF
                      ELSE
                         stop_step_det = .TRUE.      
                         y = test_y
                         f1_y = test_f1
                         f2_y = test_f2
                      END IF
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-
                            

                  
                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                  
                     x_new = y               ! the new iteration point is the auxilary point y
                     f1_new = f1_y           ! the value of f_1 at x_new
                     f2_new = f2_y           ! the value of f_2 at x_new
                     reason_for_stop = 0     ! the reason for stop is the new iteration point             
                     
                  ELSE
                     reason_for_stop = 2     ! the reason for stop is that the approximate stopping condition is satisfied (i.e. the step-length beta* < step_tol)
                     
                  END IF 
               
                  iter_counter = iter_counter -1    ! one extra iteration in the counter needs to be removed 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                               
              END IF           
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         
          
           
           END SUBROUTINE guaranteeing_clarke
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|    


        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                          QUADRATIC NORM SOLVER                                 |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE quadratic_solver(x, NF, NA, B, t, d, obj) 
            ! Solves the quadratic norm problem in the 'Clarke stationary' algorithm
            ! Calls PLQDF1 by Ladislav Luksan
            !
            ! INPUT: * 'x'               : The current iteratioin point
            !        * 'NF' and 'NA'     : The number of variables and the size of the bundle 'B', respectively    
            !        * 'B'               : The bundle 'B'
            !        * 't'               : The proximity parameter
            !
            ! OUTPUT: * 'd'              : the vector giving minimum norm
            !         * 'obj'            : the optimal value of the quadratic norm problem
            !
            ! NOTICE: The dimensions of vectors 'x' and 'd' has to be same. (i.e. the dimensio of 'x' and 'd' has to be
            !         'user_n' when SUBROUTINE quadratic_solver is used in SUBROUTINE guaranteeing_clarke.
            
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN) :: B                      ! the bundle B
               REAL(KIND=c_double), DIMENSION(:), INTENT(IN)  :: x       ! the current iteration point
               REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: d       ! the vector giving minimum norm
               
               REAL(KIND=c_double), INTENT(IN) :: t                      ! the proximity parameter t
               REAL(KIND=c_double), INTENT(OUT) :: obj                   ! the optimal value of the quadratic norm problem
               
               INTEGER(KIND=c_int), INTENT(IN)  :: NF    ! the number of variables
               INTEGER(KIND=c_int), INTENT(IN)  :: NA    ! the bundle size of B is 'give_size_b1(B) + 1' 
               
              
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=c_double), PARAMETER :: ETA0  = 1.0E-15_c_double  
               REAL(KIND=c_double), PARAMETER :: ETA2  = 1.0E-12_c_double
               REAL(KIND=c_double), PARAMETER :: ETA9  = 1.0E+60_c_double
               REAL(KIND=c_double), PARAMETER :: EPS7  = 1.0E-14_c_double
               REAL(KIND=c_double), PARAMETER :: EPS9  = 1.0E-12_c_double 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=c_double) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=c_double) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=c_double) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  

               REAL(KIND=c_double) :: u      ! help variable
               
               INTEGER(KIND=c_int) :: j, k, l      ! help variables
               
               INTEGER(KIND=c_int), PARAMETER :: NC = 0             ! number of constraints is zero
               INTEGER(KIND=c_int) :: IDECF=10, KBC=0, KBF=0, MFP=2 ! IDECF=10 diagonal matrix; KBC=0 no linear constraints; KBF=0 no simple bounds; MFP=2 optimum feasible point
               INTEGER(KIND=c_int) :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER(KIND=c_int), DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER(KIND=c_int), DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER(KIND=c_int), DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER(KIND=c_int), DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=c_double), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=c_double), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=c_double), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=c_double), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=c_double), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=c_double), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=c_double), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=c_double), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=c_double), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=c_double), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=c_double), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=c_double), DIMENSION(NF*NA) :: grad_m_b           ! subgradient matrix of B


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_c_double / t
               H = u              ! diagonal matrix
               
               grad_m_b = grad_matrix(B)  ! subgradient matrix of B 
               
               ! NO linearization errors
               AF = 0.0_c_double      
  
               ! matrix containing subgradients   
               DO j = 1, NA
                   k = (j-1)*NF
                   DO l = 1, NF
                      AG(k+l) = grad_m_b(k+l) 
                   END DO
               END DO   
                        
               !Calls PLQDF1 by Ladislav Luksan
               CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                           & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,MFP,KBF, &
                           & KBC,IDECF,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                       
                           & UMAX,GMAX,N,ITERQ)

               ! the vector giving minimum norm               
               DO j = 1, NF
                   d(j) = S(j) 
               END DO
              
              ! the optimal value of the quadratic norm problem
               obj = 0.5_c_double * DOT_PRODUCT(d, d)         
    
                   
               
           END SUBROUTINE quadratic_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        
        


        
          
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           SUBPROBLEM SOLVER                                    |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE subproblem_solver(x, NF, NA, B1, B2, t, subprob_counter) 
            ! Solves the subproblems in the original search direction problem (used in 'main iteration')
            ! Calls PLQDF1 by Ladislav Luksan
            !
            ! INPUT: * 'x'               : The current iteratioin point
            !        * 'NF' and 'NA'     : The number of variables and the size of the bundle 'B1', respectively    
            !        * 'B1' and 'B2'     : The bundle 'B1' and 'B2'
            !        * 't'               : The proximity parameter
            !
            ! OUTPUT: * 'subprob_counter'  : the number of subproblems solved
            !
            ! NOTICE: The dimensions of vectors 'x' has to be ' user_n' when SUBROUTINE subproblem_solver is used in SUBROUTINE main_iteration.
              
              IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B1                  ! bundle B_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                  ! bundle B_2
               REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: x        ! current iteration point
               
               REAL(KIND=c_double), INTENT(IN) :: t                      ! proximity parameter t
               
               INTEGER(KIND=c_int), INTENT(IN)  :: NF    ! number of variables
               INTEGER(KIND=c_int), INTENT(IN)  :: NA    ! bundle size of B1         IF ( NA = give_size_b1(B1) + 2 ) THEN aggregation is used
               
               INTEGER(KIND=c_int), INTENT(OUT) :: subprob_counter    ! number of subproblems solved    
                       
              
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=c_double), PARAMETER :: ETA0  = 1.0E-15_c_double  
               REAL(KIND=c_double), PARAMETER :: ETA2  = 1.0E-12_c_double
               REAL(KIND=c_double), PARAMETER :: ETA9  = 1.0E+60_c_double
               REAL(KIND=c_double), PARAMETER :: EPS7  = 1.0E-14_c_double
               REAL(KIND=c_double), PARAMETER :: EPS9  = 1.0E-12_c_double 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=c_double) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=c_double) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=c_double) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  
               
               REAL(KIND=c_double) :: alpha_b2     ! a linearization error of B_2
               REAL(KIND=c_double) :: u, obj, a    ! help variables
               
               INTEGER(KIND=c_int) :: i, j, k, l         ! help variables
               
               INTEGER(KIND=c_int), PARAMETER :: NC = 0             ! number of constraints is zero
               INTEGER(KIND=c_int) :: IDECF=10, KBC=0, KBF=0, MFP=2 ! IDECF=10 diagonal matrix; KBC=0 no linear constraints; KBF=0 no simple bounds; MFP=2 optimum feasible point
               INTEGER(KIND=c_int) :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER(KIND=c_int), DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER(KIND=c_int), DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER(KIND=c_int), DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER(KIND=c_int), DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=c_double), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=c_double), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=c_double), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=c_double), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=c_double), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=c_double), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=c_double), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=c_double), DIMENSION(NF) :: direction       ! Actual direction vector used (notice dimension)
               REAL(KIND=c_double), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=c_double), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=c_double), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=c_double), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=c_double), DIMENSION(NF*NA) :: grad_m_b1           ! subgradient matrix of B_1
               REAL(KIND=c_double), DIMENSION(NA) :: alpha_m_b1             ! linearization error matrix of B_1
               REAL(KIND=c_double), DIMENSION(NF) :: grad_b2                ! a subgradient of B_2


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_c_double / t
               H = u              ! diagonal matrix
               
               IF ( (give_size_b1(B1)+2) == NA) THEN     ! if this is TRUE then aggregation is in use
                   grad_m_b1 = grad_matrix_agg(B1)       ! subgradient matrix of B_1 with aggregation
                   alpha_m_b1 = lin_error_matrix_agg(B1) ! linearization error matrix of B_1 with aggregation
               ELSE
                   grad_m_b1 = grad_matrix(B1)       ! subgradient matrix of B_1
                   alpha_m_b1 = lin_error_matrix(B1) ! linearization error matrix of B_1
               END IF
               
               subprob_counter = give_size_b2(B2) 
                    
                    
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                    

                        grad_b2  = give_subgrad_b2(B2, i)  ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i)   ! linearization error of B_2 in the subproblem i

                        DO j = 1, NA
                           k = (j-1)*NF
                           DO l = 1, NF
                               AG(k+l) = grad_m_b1(k+l) - grad_b2(l)
                           END DO
                           AF(j) = - alpha_m_b1(j) + alpha_b2 
                        END DO
                        
                        !Calls PLQDF1 by Ladislav Luksan
                        CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,MFP,KBF, &
                              & KBC,IDECF,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
                              & UMAX,GMAX,N,ITERQ)

                              
                        DO j = 1, NF
                           direction(j) = S(j) 
                        END DO
              
                        a = DOT_PRODUCT(direction, direction)           

                        obj =   XNORM + (a * u) / 2
                        CALL add_solution(B2, i , direction, XNORM, obj )   
                    
                    END DO subproblems1
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      

               subprob_counter = subprob_counter + 1 
               
               
           END SUBROUTINE subproblem_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

        
      END MODULE dbdc

!----------------------------------------------------------------------------------------
!
!                           *************  LMBM  *************
!
!----------------------------------------------------------------------------------------
MODULE param  ! Parameters
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

! Intrinsic Functions
  INTRINSIC TINY,HUGE

! Parameters
  INTEGER(KIND=c_int), PARAMETER, PUBLIC :: maxeps = 20_c_int, maxnrs = 2000_c_int
  REAL(KIND=c_double), PARAMETER, PUBLIC :: &
       zero    = 0.0_c_double,    & ! 
       half    = 0.5_c_double,    & ! 
       one     = 1.0_c_double,    & ! 
       large   = 3.40282347*10.**38,  & !
       small   = 1.17549435*10.**(-38)     ! 

END MODULE param


MODULE exe_time  ! Execution time.
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

  PUBLIC :: getime

CONTAINS
  SUBROUTINE getime(tim)  ! Execution time.
  IMPLICIT NONE
      
! Scalar Arguments
    REAL(KIND=c_double), INTENT(OUT):: tim  ! Current time, REAL argument.

! Intrinsic Functions
    INTRINSIC CPU_TIME

    CALL CPU_TIME(tim)

  END SUBROUTINE getime

END MODULE exe_time


MODULE lmbm_sub  ! Subprograms for lmbm

  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

! MODULE lmbm_sub includes the following subroutines (S) and functions (F).
  PUBLIC :: &
       vdot, &   ! F Dot product of two vectors.
       vneg, &   ! S Change the signs of vector elements.
       copy, &   ! S Copying a vector.
       copy2, &  ! S Copying two vectors.
       xdiffy, & ! S Difference of two vectors z:= x - y.
       xsumy, &  ! S Sum of two vectors z:= x + y.
       scdiff, & ! S Difference of the scaled vector and a vector z:= a*x - y.
       scsum, &  ! S Sum of a vector and the scaled vector z:= y + a*x.
       vxdiag, & ! S Vector is multiplied by a diagonal matrix y:=d*x.
       symax, &  ! S Multiplication of a dense symmetric matrix A by a vector x.
       cwmaxv, & ! S Multiplication of a columnwise stored dense rectangular matrix by a vector.
       rwaxv2, & ! S Multiplication of two rowwise stored dense rectangular  
                 !   matrices A and B by vectors X and Y.
       trlieq, & ! S Solving x from linear equation u*x=y or u'*x=y, 
                 !   where u is an upper triangular matrix.
       lineq, &  ! S Solver from linear equation.
       calq      ! S Solving x from linear equation A*x=y. Contains:
                 !     S mxdpgf   Gill-Murray decomposition of a dense symmetric matrix.

CONTAINS

  FUNCTION vdot(n,x,y) RESULT(xty) ! Dot product of two vectors.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    REAL(KIND=c_double) xty
    INTEGER :: i

    xty = zero
    DO i = 1,n
       xty = xty + x(i)*y(i)
    END DO

  END FUNCTION vdot

  SUBROUTINE vneg(n,x,y)  ! Change the signs of vector elements.
    IMPLICIT NONE
      
! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= -x.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO i = 1,n
       y(i) = -x(i)
    END DO

  END SUBROUTINE vneg

  SUBROUTINE scalex(n,a,x,y)  ! Scaling a vector y:= a*x.
    IMPLICIT NONE
      
! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(IN) :: &
         a           ! Scaling parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO i = 1,n
       y(i) = a*x(i)
    END DO
      
  END SUBROUTINE scalex

  SUBROUTINE xdiffy(n,x,y,z)  ! Difference of two vectors z:= x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= x - y.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO  i = 1,n
       z(i) = x(i) - y(i)
    END DO
 
  END SUBROUTINE xdiffy

  SUBROUTINE xsumy(n,x,y,z)  ! Sum of two vectors z:= x + y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= x + y.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO  i = 1,n
       z(i) = x(i) + y(i)
    END DO

  END SUBROUTINE xsumy

  SUBROUTINE scdiff(n,a,x,y,z)  ! Difference of the scaled vector and a vector z:= a*x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x - y.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO  i = 1,n
       z(i) = a*x(i) - y(i)
    END DO
 
  END SUBROUTINE scdiff

  SUBROUTINE scsum(n,a,x,y,z)  ! Sum of a vector and the scaled vector z:= y + a*x.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x + y.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO  i = 1,n
       z(i) = a*x(i) + y(i)
    END DO

  END SUBROUTINE scsum

  SUBROUTINE copy(n,x,y)  ! Copying of vector Y:= X.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO i = 1,n
       y(i) = x(i)
    END DO
      
  END SUBROUTINE copy

  SUBROUTINE copy2(n,x,y,z,v)  ! Copying of two vectors: y:=x, v:=z.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,z         ! Input vectors.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y,v         ! Output vectors.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO i = 1,n
       y(i) = x(i)
       v(i) = z(i)
    END DO

  END SUBROUTINE copy2

  SUBROUTINE vxdiag(n,d,x,y)  ! Vector is multiplied by a diagonal matrix y:=d*x.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector.
         d           ! Diagonal matrix stored as a vector with n elements.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= d*x.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER(KIND=c_int) :: i

    DO  i = 1,n
       y(i) = x(i)*d(i)
    END DO
 
  END SUBROUTINE vxdiag

  SUBROUTINE symax(n,m,iold,a,x,y)  ! Multiplication of a dense symmetric matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector stored in a circular order.
         a           ! Dense symmetric matrix stored in the packed form: a(n*(n+1)/2).
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x. Vector y has the same circular order than x.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &         ! Order of matrix A.
         m, &         ! Length of vector x, m >= n, note that only n
                      ! components from vector x are used.
         iold         ! Index, which controlls the circular order of
                      ! the vector x.

! Local Scalars
    INTEGER(KIND=c_int) :: i,j,k,l

    DO j=1,n
       l=j+iold-1
       IF (l > m) l=l-m
       y(l) = zero
       k=l
       DO i=j,n
          y(l) = a((i-1)*i/2+j)*x(k)+y(l)
          k=k+1
          IF (k > m) k=k-m
       END DO
    END DO

    DO j=2,n
       l=j+iold-1
       IF (l > m) l=l-m
       k=iold
       DO i=1,j-1
          IF (k > m) k=k-m
          y(l) = a((j-1)*j/2+i)*x(k)+y(l)
          k=k+1
       END DO
    END DO
      
  END SUBROUTINE symax

  SUBROUTINE cwmaxv(n,m,a,x,y)  ! Multiplication of a columnwise stored dense 
                                ! rectangular matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector (dimension m).
         a           ! Rectangular matrix stored columnwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector equal to s*a*x. If m = 0 y is a zero vector. 

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &        ! Number of rows of the matrix A.
         m           ! Number of columns of the matrix A.

! Local Scalars
    INTEGER(KIND=c_int) :: i,j,k
      
    DO i = 1,n
       y(i) = zero
    END DO

    k = 1
    DO j = 1,m
       CALL scsum(n,x(j),a(k:),y,y)
       k = k + n
    END DO

  END SUBROUTINE cwmaxv

  SUBROUTINE rwaxv2(n,m,a,b,x,y,v,w)  ! Multiplication of two rowwise stored dense rectangular  
                                      ! matrices A and B by vectors X and Y.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x,y, &      ! Input vectors (dimension n).
         a,b         ! Rectangular matrices stored rowwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         v,w         ! Output vectors v=a*x and w=b*y. 

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &        ! Number of columns of the matrices A and B.
         m           ! Number of rows of the matrices A and B.

! Local Scalars
    REAL(KIND=c_double) :: tmp1,tmp2
    INTEGER(KIND=c_int) :: i,j,k
      
    k = 0
    DO i = 1,m
       tmp1 = zero
       tmp2 = zero
       DO j = 1,n
          tmp1 = tmp1 + a(k+j)*x(j)
          tmp2 = tmp2 + b(k+j)*y(j)
       END DO
       v(i) = tmp1
       w(i) = tmp2
       k = k + n
    END DO

  END SUBROUTINE rwaxv2


  SUBROUTINE trlieq(n,m,iold,u,x,y,job,ierr)  ! Solving x from linear equation u*x=y or u'*x=y, 
                                              ! where u is an upper triangular matrix.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order.
         u           ! Triangular matrix.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector y has the same circular order than x.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &         ! Order of matrix U.
         m, &         ! Length of vectors x and y, m >= n, note that only n
                      ! components from vectors are used.
         iold, &      ! Index, which controlls the circular order of
                      ! the vectors x and y.
         job          ! Option:
                      !   0  - x:=(u')**(-1)*y, u upper triangular.
                      !   1  - x:=u**(-1)*y, u upper triangular.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
         ierr         ! Error indicador: 
                      !   0   - Everything is ok.
                      !  -3   - Error; 0 at diagonal.

! Local Scalars
    INTEGER(KIND=c_int) :: i,ii,ij,j,k,l,ji

! Intrinsic Functions
    INTRINSIC ABS
      
    ierr = -3
      
    DO i=1,m
       x(i)=y(i)
    END DO
      
    IF (job == 0) THEN
     
! x=u'**(-1)*y, u' = [u1         ] is lower triangular.
!                    [u2 u3      ]
!                    [u4 u5 u6   ]
!                    [.  .  .  . ]
         
       ii = 0
       DO  i = 1,n
          ii=ii+i
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          x(l) = x(l)/u(ii)
          DO j = i+1,n
             ji = (j-1)*j/2+i
             k=j+iold-1
             IF (k > m) k=k-m
             x(k) = x(k) - u(ji)*x(l)
          END DO
       END DO
             
         
    ELSE IF (job == 1) THEN
     
! x=u**(-1)*y, u = [u1 u2 u4 . ] is upper triangular.
!                  [   u3 u5 . ]
!                  [      u6 . ]
!                  [         . ]
         
       ii = n* (n+1)/2
       DO i = n,1,-1
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          ij = ii
          DO j = i + 1,n
             k=j+iold-1
             IF (k > m) k=k-m
             ij = ij + j - 1
             x(l) = x(l) - u(ij)*x(k)
          END DO
          x(l)=x(l)/u(ii)
          ii = ii - i
       END DO
         
         
    ELSE
         
       RETURN
    END IF
      
    ierr = 0

  END SUBROUTINE trlieq
      
  SUBROUTINE lineq(n,m,iold,a,x,y,ierr)  ! Solving X from linear equation A*X=Y. 
                                         ! Positive definite matrix A+E is given using 
                                         ! the factorization A+E=L*D*L' obtained by the
                                         ! subroutine mxdpgf.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order (dimension m).
         a           ! Factorization a+e=l*d*l' obtained by the subroutine mxdpgf.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold        ! Index, which controlls the circular order of
                     ! the vectors x and y.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
         ierr        ! Error indicador: 
                     !   0   - Everything is ok.
                     !  -2   - Error; indefinite matrix.

! Local Scalars
    INTEGER(KIND=c_int) :: i,ii,ij,j,k,l


    ierr = -2
      
! Phase 1: x=l**(-1)*x

    ij = 0
    DO i = 1,n
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = y(l)
         
       DO j = 1,i - 1
          ij = ij + 1
          k=j+iold-1
          IF (k > m) k=k-m
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ij = ij + 1
    END DO

! Phase 2 : x:=d**(-1)*x

    ii = 0
    DO i = 1,n
       ii = ii + i
       IF (a(ii) <= small) RETURN
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = x(l)/a(ii)
    END DO

! Phase 3 : x:=trans(l)**(-1)*x

    ii = n* (n-1)/2
    DO i = n - 1,1,-1
       ij = ii
       l=i+iold-1
       IF (l > m) l=l-m
       DO j = i + 1,n
          k=j+iold-1
          IF (k > m) k=k-m
          ij = ij + j - 1
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ii = ii - i
    END DO

    ierr = 0

  END SUBROUTINE lineq

  SUBROUTINE calq(n,m,iold,a,x,y,iprint)  ! Solving x from linear equation A*x=y.
    USE param, ONLY : zero,small,one
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         y           ! Input vector stored in a circular order (dimension m).
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         a           ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold, &     ! Index, which controlls the circular order of
                     ! the vectors x and y.
         iprint      ! Printout specification.
!    INTEGER, INTENT(OUT) :: &
!         ierr        ! Error indicador: 
!                     !   0   - Everything is ok.
!                     !  -2   - Error; indefinite matrix.

! Local Scalars
    REAL(KIND=c_double) :: eta,bet
    INTEGER(KIND=c_int) :: inf,ierr

      
    eta = small+small
      
    CALL mxdpgf(n,a,inf,eta,bet)


    IF (iprint == 2) THEN
       IF (inf < 0) THEN
          WRITE (6,FMT='(1X,''Warning: Insufficiently positive'' &
               '' definite matrix detected. '')')
          WRITE (6,FMT='(1X,''Correction added.'')')
         
       ELSE IF (inf > 0) THEN
          WRITE (6,FMT='(1X,''Warning: Indefinite'' &
            '' matrix detected. Correction added.'')')
       END IF
    END IF
      
    CALL lineq(n,m,iold,a,x,y,ierr)
    IF (ierr /= 0) THEN
       PRINT*,'hihu'
       IF (iprint == 2) THEN
          WRITE (6,FMT='(1X,''Warning: Indefinite matrix detected. '')')
       END IF
    END IF

  CONTAINS

    SUBROUTINE mxdpgf(n,a,inf,alf,tau)  ! Factorization A+E=L*D*trans(L) of a dense symmetric positive
                                        ! definite matrix A+E, where D and E are diagonal positive 
                                        ! definite matrices and L is a lower triangular matrix. 
                                        ! If A is sufficiently positive definite then E=0.
      
! Array Arguments
      REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
           a         ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).

! Scalar Arguments
      REAL(KIND=c_double), INTENT(INOUT) :: &
           alf       ! On input a desired tolerance for positive definiteness. 
                     ! On output the most negative diagonal element used in the factorization
                     ! process (if inf>0).
      REAL(KIND=c_double), INTENT(OUT) :: &
           tau       ! Maximum diagonal element of matrix E.

      INTEGER(KIND=c_int), INTENT(IN) :: &
           n         ! Order of matrix a.
      INTEGER(KIND=c_int), INTENT(OUT) :: &
           inf       ! An information obtained in the factorization process:
                     !    inf=0  - A is sufficiently positive definite and E=0. 
                     !    inf<0  - A is not sufficiently positive definite and E>0.
                     !    inf>0  - A is indefinite and inf is an index of the most negative 
                     !             diagonal element used in the factorization process.

! Local Scalars
      REAL(KIND=c_double) :: bet,del,gam,rho,sig,tol
      INTEGER(KIND=c_int) :: i,ij,ik,j,k,kj,kk,l

! Intrinsic Functions
      INTRINSIC ABS,MAX

      l = 0
      inf = 0
      tol = alf
      

! Estimation of the matrix norm

      alf = zero
      bet = zero
      gam = zero
      tau = zero
      kk = 0

      DO k = 1,n
         kk = kk + k
         bet = MAX(bet,ABS(a(kk)))
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
      END DO
      bet = MAX(tol,bet,gam/n)

      del = tol*MAX(bet,one)
      kk = 0
      DO k = 1,n
         kk = kk + k

!     Determination of a diagonal correction

         sig = a(kk)
         IF (alf > sig) THEN
            alf = sig
            l = k
         END IF

         gam = zero
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
         gam = gam*gam
         rho = MAX(ABS(sig),gam/bet,del)
         IF (tau < rho-sig) THEN
            tau = rho - sig
            inf = -1
         END IF
         
! Gaussian elimination

         a(kk) = rho
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = a(kj)
            a(kj) = gam/rho
            ik = kk
            ij = kj
            DO i = k + 1,j
               ik = ik + i - 1
               ij = ij + 1
               a(ij) = a(ij) - a(ik)*gam
            END DO
         END DO
      END DO
      IF (l > 0 .AND. ABS(alf) > del) inf = l
      
    END SUBROUTINE mxdpgf
  END SUBROUTINE calq

END MODULE lmbm_sub


MODULE initializat  ! Initialization of parameters and x for LDGBM and LMBM

  USE, INTRINSIC :: iso_c_binding
  USE functions                  ! Contains INFORMATION from the USER

  IMPLICIT NONE

    ! Parameters 
    INTEGER(KIND=c_int), SAVE :: &
    
        na, &                         ! Size of the bundle na >= 2.
        !na     =    5_c_int, &       ! Size of the bundle na >= 2.
        mcu, &                        ! Upper limit for maximum number of stored corrections, mcu >= 3.
        !mcu    =    7_c_int, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
        mcinit, &                     ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
        !mcinit =    7_c_int, &       ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                      ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                      ! However, the value mcinit = 7 is recommented.
        inma    = 3_c_int, &          ! Selection of line search method:
                                      !   inma = 0, Armijo line search, Not here
                                      !   inma = 1, nonmonotone Armijo line search. Not here
                                      !   inma = 2, weak Wolfe line search.
                                      !   inma = 3, nonmonotone  weak Wolfe line search.
        mnma    = 10_c_int, &         ! Maximum number of function values used in nonmonotone line search.
        maxnin  = 20_c_int            ! Maximum number of interpolations, maxnin >= 0.
                                      ! The value maxnin = 2-20 is recommented with inma=0,
                                      ! maxnin >= 20 with inma=1 and 3, and maxnin =200 with inma=2.
                                      ! For example:
                                      !   inma = 0, maxin = 20.
                                      !   inma = 1, mnma = 20, maxin = 30.
                                      !   inma = 2, maxnin = 200.
                                      !   inma = 3, mnma=10, maxnin = 20.

    REAL(KIND=c_double), PARAMETER :: &
        time =  1000.0_c_double        !  Maximum CPU-time in seconds. If time <= 0.0 the maximum time 
                                       !  is ignored. REAL argument.

    ! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
    REAL(KIND=c_double), SAVE :: &
        tolf, &                          ! Tolerance for change of function values (default = 1.0E-8).
        !tolf  = 1.0E-5_c_double, &      ! Tolerance for change of function values (default = 1.0E-8).
        tolf2, &                         ! Second tolerance for change of function values.
        !tolf2 = 1.0E+4_c_double, &      ! Second tolerance for change of function values.
                                         !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                         !   criterion will be ignored. 
                                         !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
        tolb  =  -3.40282347*10.**38, &  ! + 1.17549435*10.**(-38), &  ! Tolerance for the function value (default = -large).
        
		tolg, &                          ! Tolerance for the first termination criterion (default = 1.0E-6).
		!tolg  = 1.0E-5_c_double, &      ! Tolerance for the first termination criterion (default = 1.0E-6).
        tolg2, &                         ! Tolerance for the second termination criterion (default = tolg). clustering code small data
        !tolg2 = 1.0E-3_c_double, &      ! Tolerance for the second termination criterion (default = tolg). clustering code small data
        !tolg2 = 1.0E-5_c_double, &      ! Tolerance for the second termination criterion (default = tolg).
    
        ! eta=0.5 on toistaiseksi antanut parhaita tuloksia
        eta, &                           ! Distance measure parameter, eta >= 0. 
        !eta   =    0.5_c_double, &      ! Distance measure parameter, eta >= 0. 
        !eta   =    1.0E-4_c_double, &   ! Distance measure parameter, eta >= 0. 
                                         !   - If eta < 0  the default value 0.5 will be used. 
          
        epsL, &    		                  ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.
        !epsl  =    1.0E-4_c_double, &    ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4.
        !epsl  =    0.24E+00_c_double, &  ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4.
        !epsl  =    0.125E+00_c_double, & ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4.
    
        !xmax  =    1.5_c_double      ! Maximum stepsize, 1 < XMAX (default = 1.5).
        !xmax  =    10.0_c_double     ! Maximum stepsize, 1 < XMAX (default = 1.5).
        xmax  =    100.0_c_double     ! Maximum stepsize, 1 < XMAX (default = 1.5).

    INTEGER(KIND=c_int), SAVE :: n        ! The direction

    INTEGER(KIND=c_int), SAVE :: nproblem ! The solved problem 


    ! Integer parameters (if value <= 0 the default value of the parameter will be used).
    INTEGER, SAVE :: &
    
        mittt    = 5000, &        ! Maximun number of iterations (default = 10000).
        !mittt    = 5000, &       ! Maximun number of iterations (default = 10000). This is in original.
        mfe    = 50000, &         ! Maximun number of function evaluations (default = n*mittt).
        mtesf  =     10, &        ! Maximum number of iterations with changes of
                                  ! function values smaller than tolf (default = 10).
        iiprint =      0, &       ! Printout specification:
                                  !     0  - Only the error messages.
                                  !     1  - The final values of the objective function 
                                  !          (default used if iiprint < -1).
                                  !     2  - The final values of the objective function and the 
                                  !          most serious warning messages.
                                  !     3  - The whole final solution. 
                                  !     4  - At each iteration values of the objective function.
                                  !     5  - At each iteration the whole solution
        method =      0, &        ! Selection of the method:
                                  !     0  - Limited memory bundle method (default).
                                  !     1  - L-BFGS bundle method.
                                  !     2  - Limited memory discrete gradient bundle method.
        iscale =      0           ! Selection of the scaling:
                                  !     0  - Scaling at every iteration with STU/UTU (default).
                                  !     1  - Scaling at every iteration with STS/STU.
                                  !     2  - Interval scaling with STU/UTU.
                                  !     3  - Interval scaling with STS/STU.
                                  !     4  - Preliminary scaling with STU/UTU.
                                  !     5  - Preliminary scaling with STS/STU.
                                  !     6  - No scaling.      


    REAL(KIND=c_double), DIMENSION(:), ALLOCATABLE, SAVE :: x   ! Vector of variables (initialize x in subroutine init_x)
    
    TYPE(set_info), SAVE :: LMBM_set          ! The set of information               
    


CONTAINS

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!-------------------------------------------------------------------------------------------
    SUBROUTINE allocate_xn(user_N) 
        !
        ! 
        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        INTEGER(KIND=c_int), INTENT(IN) :: user_N                 ! the of features/inputs in a data point   
        !************************************************************************************
        n=user_N
   
        ALLOCATE(x(n))
  
    END SUBROUTINE allocate_xn
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  


!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE deallocate_x() 
        !
        ! 
        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        !***********************************************************************************
  
        DEALLOCATE(x)
  
    END SUBROUTINE deallocate_x  
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
  
  
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE init_x(x_apu)  ! User supplied subroutine for initialization of vector x(n) and the solved problem

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: x_apu 
                   
        !***********************************************************************************
        INTEGER(KIND=c_int) :: i
    
        DO i = 1, n
            x(i) = x_apu(i)
        END DO   
    
    END SUBROUTINE init_x
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
  
  
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE init_LMBMinfo(problem, set)  ! User supplied subroutine for initialization of vector x(n) and the solved problem

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        INTEGER(KIND=c_int), INTENT(IN) :: problem
        TYPE(set_info), INTENT(INOUT) :: set                         ! The set of information                       
        !***********************************************************************************
        INTEGER(KIND=c_int) :: i
    
        nproblem = problem
        LMBM_set = set
    
    END SUBROUTINE init_LMBMinfo
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  


!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
!------------------------------------------------------------------------------------------- 
    SUBROUTINE copy_x(x_apu)  ! User supplied subroutine to copy the vector x(n).

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: x_apu 
        !***********************************************************************************
        INTEGER :: i
    
        DO i = 1, n
            x_apu(i) = x(i)
        END DO   
    
    END SUBROUTINE copy_x  
!-------------------------------------------------------------------------------------------    
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/


!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
!------------------------------------------------------------------------------------------- 
    SUBROUTINE set_rho_LMBM(value_rho)  ! User supplied subroutine to pass rho value

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=c_double), INTENT(IN) :: value_rho 
        !***********************************************************************************
        INTEGER :: i
    
       ! ***** penalization parameter 'rho' for L0-norm *****
       IF ( (value_rho<0.0_c_double) ) THEN
         LMBM_set%user_rho = 0.0_c_double
       ELSE       
         LMBM_set%user_rho = value_rho 
       END IF  
       
    END SUBROUTINE set_rho_LMBM
!-------------------------------------------------------------------------------------------    
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

  

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
!------------------------------------------------------------------------------------------- 
    SUBROUTINE set_lambda_LMBM(value_lambda)  ! User supplied subroutine to pass rho value

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=c_double), INTENT(IN) :: value_lambda 
        !***********************************************************************************
        INTEGER :: i
    
       ! ***** penalization parameter 'rho' for L1-norm *****
       IF ( (value_lambda<0.0_c_double) ) THEN
         LMBM_set%user_lambda = 0.0_c_double
       ELSE       
         LMBM_set%user_lambda = value_lambda 
       END IF  
       
    END SUBROUTINE set_lambda_LMBM
!-------------------------------------------------------------------------------------------    
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
     
    
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!-------------------------------------------------------------------------------------------  
    SUBROUTINE init_par(in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, in_tolg2, in_eta, in_epsL)  
                             ! User supplied subroutine for further initialization of parameters (when needed).
                             ! This subroutine may be usefull if you want to run some test problems in loop.
                             ! May be left empty.
 
        IMPLICIT NONE
        INTEGER(KIND=c_int) ::  in_na, in_mcu, in_mcinit
        REAL(KIND=c_double) ::  in_tolf, in_tolf2
        REAL(KIND=c_double) ::  in_tolg, in_tolg2
        REAL(KIND=c_double) ::  in_eta, in_epsL
        !na     =    5_c_int, &        ! Size of the bundle na >= 2.
        !mcu    =    7_c_int, &        ! Upper limit for maximum number of stored corrections, mcu >= 3.
        !mcinit =    7_c_int, &        ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                       ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                       ! However, the value mcinit = 7 is recommented.
        !tolf  = 1.0E-5_c_double, &       ! Tolerance for change of function values (default = 1.0E-8).
        !tolf2 = 1.0E+4_c_double, &       ! Second tolerance for change of function values.
                                         !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                         !   criterion will be ignored. 
                                         !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
        !tolg  = 1.0E-5_c_double, &        ! Tolerance for the first termination criterion (default = 1.0E-6).
        !tolg2 = 1.0E-3_c_double, &        ! Tolerance for the second termination criterion (default = tolg). clustering code small data
        !eta   =    0.5_c_double, &     s  ! Distance measure parameter, eta >= 0. 
        !epsL  =    0.24E+00_c_double, &   ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.

        IF ( in_na >= 2) THEN 
            na = in_na
        ELSE
            na = 3
        END IF      
        
        IF ( in_mcu >= 3) THEN 
            mcu = in_mcu
        ELSE
            mcu = 3
        END IF    
        
        IF ( (in_mcinit >= 3) .AND. (in_mcinit <= mcu)) THEN 
            mcinit = in_mcinit
        ELSE
            mcinit = mcu
        END IF  
        
        IF ( in_tolf > 0.0_c_double) THEN 
            tolf = in_tolf
        ELSE
            tolf = 1.0E-5_c_double
        END IF          

        IF ( in_tolf2 == 0.0_c_double) THEN 
            tolf2 = 1.0E+4_c_double
        ELSE
            tolf2 = in_tolf2
        END IF      

        IF ( in_tolg > 0.0_c_double) THEN 
            tolg = in_tolg
        ELSE
            tolg = 1.0E-5_c_double
        END IF          

        IF ( in_tolg2 > 0.0_c_double) THEN 
            tolg2 = in_tolg2
        ELSE
            tolg2 = tolg
        END IF          
        
        IF ( in_eta >= 0.0_c_double) THEN 
            eta = in_eta
        ELSE
            eta = 0.5_c_double
        END IF              
        
        IF ( (in_epsL > 0.0_c_double) .AND. (in_epsL < 0.25_c_double)) THEN 
            epsL = in_epsL
        ELSE
            epsL = 1.0E-4_c_double
        END IF  

    END SUBROUTINE init_par
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

END MODULE initializat


MODULE obj_fun  ! Computation of the value and the subgradient of the 
                ! objective function (supplied by the user).

  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

  PUBLIC :: &
       myf, &    ! Computation of the value of the objective. User spesified function. 
       myg       ! Computation of the subgradient of the objective. User spesified function. 

CONTAINS
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE myf *                                               *
!*                                                                      *
!*     Computation of the value of the objective. User spesified        *
!*     function.                                                        *
!*                                                                      *
!************************************************************************
     
  SUBROUTINE myf(n,x,f,iterm,nproblem)

    USE functions, ONLY : func       ! Computation of the value for problem next.
                                     
    USE initializat, ONLY : LMBM_set ! Data info

    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(n), INTENT(IN) :: &
         x  ! Vector of variables.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(OUT) :: f          ! Value of the objective function. 
    INTEGER(KIND=c_int), INTENT(IN) :: n           ! Number of variables.
    INTEGER(KIND=c_int), INTENT(IN) :: nproblem    ! The number of the problem solved   
    INTEGER(KIND=c_int), INTENT(OUT) :: iterm      ! Cause of termination:
                                                   !   0  - Everything is ok.
                                                   !  -3  - Failure in function calculations 
                                                   !        (assigned by the user).

    iterm = 0

! Function evaluation (give your function here).

    CALL func(n,x,f,nproblem,LMBM_set)         ! Only the value of the objective
    !IF (next < 1) iterm = -3 ! Error checking. 
      
  END SUBROUTINE myf

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE MYG *                                               *
!*                                                                      *
!*     Computation of the subgradient of the objective function. User   *
!*     spesified function. If discrete gradients are used, this         *
!*     subroutine may be left empty.                                    *
!*                                                                      *
!************************************************************************
     
  SUBROUTINE myg(n,x,g,iterm,nproblem)

    USE functions, ONLY : subgra       ! Computation of the subgradient for problem next.
    USE initializat, ONLY : LMBM_set ! Data info

    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
    REAL(KIND=c_double), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: n            ! Number of variables.
    INTEGER(KIND=c_int), INTENT(IN) :: nproblem     ! The number of the problem solved
    INTEGER(KIND=c_int), INTENT(OUT) :: iterm       ! Cause of termination:
                                                    !   0  - Everything is ok.
                                                    !  -3  - Failure in subgradient calculations 
                                                    !        (assigned by the user).

    iterm = 0

! Subgradient evaluation (give your function here).

    CALL subgra(n,x,g,nproblem,LMBM_set)         ! Only the subgradient
    !IF (next < 1) iterm = -3  ! Error checking. 
      
  END SUBROUTINE myg

END MODULE obj_fun


MODULE lmbm_mod  ! Limited memory bundle method

  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

! MODULE lmbm_mod includes the following subroutines (S) and functions (F).
  PUBLIC :: &
       init_lmbm, &  ! S Initialization for limited memory bundle subroutine.
       lmbm          ! S Limited memory bundle subroutine for nonsmooth 
                     !   large-scale optimization. Contains:
                     !     S restar  Initialization and reinitialization.
                     !     S dobun   Bundle construction.
  PRIVATE :: &
       indic1, &     ! S Initialization of indices.
       tinit, &      ! S Calculation of initial step size. Contains:
                     !     S destep  Stepsize selection using polyhedral 
                     !                 approximation for descent steps.
                     !     S nulstep Stepsize selection using polyhedral 
                     !                 approximation for null steps.
       lls, &        ! S Line search. Contains:
                     !     F qint    Quadratic interpolation.
       nmlls, &      ! S Nonmonotonic Weak Wolfe line search. Contains:
                      !     F qint    Quadratic interpolation.
       dlbfgs, &     ! S Computing the search direction by limited memory 
                     !   BFGS update. Contains:
                     !     F sclpar  Selection of scaling parameter.
       dlskip, &     ! S Skipping the updates and computing the search
                     !   direction by limited memory BFGS update.    
       dlsr1, &      ! S Computing the search direction by limited
                     !   memory SR1 update.    
       agbfgs, &     ! S Simplified subgradient aggregation.
       aggsr1, &     ! S Subgradient aggregation.
       agskip, &     ! S Subgradient aggregation using BFGS update.
       wprint, &     ! S Printout the error and warning messages.
       rprint        ! S Printout the results.

CONTAINS
!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE init_lmbm *                                        *
!*                                                                     *
!*     Initialization for limited memory bundle subroutine for         *
!*     large-scale unconstrained nonsmooth optimization.               *
!*                                                                     *
!***********************************************************************
      
  SUBROUTINE init_lmbm(mc,iterm)
    USE param, ONLY : zero,half,small,large
    USE initializat, ONLY : &
         n, &         ! Number of variables.
         na, &        ! Maximum bundle dimension, na >= 2.
         mcu, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
         iiprint, &    ! Printout specification (see initializat for more details).
         method, &    ! Selection of the method (see initializat for more details).
         iscale, &    ! Selection of the scaling (see initializat for more details).
         tolf, &      ! Tolerance for change of function values.
         tolf2, &     ! Second tolerance for change of function values.
         tolb, &      ! Tolerance for the function value.
         tolg, &      ! Tolerance for the first termination criterion. 
         tolg2, &     ! Tolerance for the second termination criterion.
         xmax, &      ! Maximum stepsize, 1 < XMAX.
         eta, &       ! Distance measure parameter, eta >= 0. 
         epsl, &      ! Line search parameter, 0 < epsl < 0.25.
         mtesf, &     ! Maximum number of iterations with changes of
                      ! function values smaller than tolf.
         mittt, &       ! Maximun number of iterations.
         mfe          ! Maximun number of function evaluations.

    IMPLICIT NONE

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(INOUT) :: mc     ! Initial maximum number of stored corrections.
    INTEGER(KIND=c_int), INTENT(OUT) :: iterm    ! Cause of termination:
                                                 !   0  - Everything is ok.
                                                 !  -5  - Invalid input parameters.
 
! Initialization and error checking

    IF (iiprint < -1) iiprint  = 1     ! Printout specification. 
    iterm = 0
      
    IF (n <= 0) THEN
       iterm = -5
       CALL wprint(iterm,iiprint,1)
       RETURN 
    END IF

    IF (na < 2) THEN
       iterm = -5
       CALL wprint(iterm,iiprint,3)
       RETURN
    END IF

    IF (epsl >= 0.25_c_double) THEN
       iterm = -5
       CALL wprint(iterm,iiprint,4)
       RETURN
    END IF

    IF (mcu <= 3) THEN
       iterm = -5
       CALL wprint(iterm,iiprint,2)
       RETURN
    END IF

    IF (mc > mcu) THEN
       mc = mcu
       CALL wprint(iterm,iiprint,-1)
    END IF

! Default values

    IF (mc    <= 0) mc       = 3               ! Initial maximum number of corrections.
    IF (mittt   <= 0) mittt      = 10000           ! Maximum number of iterations.
    IF (mfe   <= 0) mfe      = n*mittt           ! Maximum number of function evaluations.
    IF (tolf  <= zero) tolf  = 1.0E-08_c_double    ! Tolerance for change of function values.
    IF (tolf2 == zero) tolf2 = 1.0E+04_c_double    ! Second tolerance for change of function values.
    IF (tolb  == zero) tolb  = -large + small  ! Tolerance for the function value.
    IF (tolg  <= zero) tolg  = 1.0E-06_c_double    ! Tolerance for the first termination criterion.
    IF (tolg2 <= zero) tolg2 = tolg            ! Tolerance for the second termination criterion.
    IF (xmax  <= zero) xmax  = 1.5_c_double        ! Maximum stepsize.
    IF (eta   <  zero) eta   = half            ! Distance measure parameter
    IF (epsl  <= zero) epsl  = 1.0E-04_c_double    ! Line search parameter,
    IF (mtesf <= 0) mtesf    = 10              ! Maximum number of iterations with changes 
                                               ! of function values smaller than tolf.
    IF (method > 1 .OR. method < 0) method = 0 ! Selection of the method.
    IF (iscale > 6 .OR. iscale < 0) iscale = 0 ! Selection of the scaling.

     
  END SUBROUTINE init_lmbm

!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE lmbm *                                             *
!*                                                                     *
!*     Limited memory bundle subroutine for nonsmooth optimization.    *
!*                                                                     *
!***********************************************************************
      
  SUBROUTINE lmbm(mc,f,nit,nfe,nge,iterm,strtim)
    USE param, ONLY : small,large,zero,half,one
    USE exe_time, ONLY : getime  ! Execution time.
    USE initializat, ONLY : &
         n, &         ! Number of variables.
         x, &         ! Vector of variables
         nproblem, &  ! The solved problem
         na, &        ! Maximum bundle dimension, na >= 2.
         mcu, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
         inma, &      ! Selection of line search method.
         mnma, &      ! Maximum number of function values used in nonmonotone line search.
         iiprint, &    ! Printout specification (see initializat for more details).
         method, &    ! Selection of the method (see initializat for more details).
         iscale, &    ! Selection of the scaling (see initializat for more details).
         tolf, &      ! Tolerance for change of function values.
         tolf2, &     ! Second tolerance for change of function values.
         tolb, &      ! Tolerance for the function value.
         tolg, &      ! Tolerance for the first termination criterion. 
         tolg2, &     ! Tolerance for the second termination criterion.
         xmax, &      ! Maximum stepsize, 1 < XMAX.
         eta, &       ! Distance measure parameter, eta >= 0. 
         epsl, &      ! Line search parameter, 0 < epsl < 0.25.
         mtesf, &     ! Maximum number of iterations with changes of
                      ! function values smaller than tolf.
         mittt, &     ! Maximun number of iterations.
         mfe, &       ! Maximun number of function evaluations.
         time         ! Maximum time
    USE obj_fun, ONLY : &
         myf, &       ! Computation of the value of the objective function. 
         myg          ! Computation of the subgradient of the objective function. 
    USE lmbm_sub, ONLY : &
         vdot, &      ! Dot product.
         vneg, &      ! Copying of a vector with change of the sign.
         xdiffy, &    ! Difference of two vectors.
         copy, &      ! Copying of a vector.
         copy2        ! Copying of two vectors.


    IMPLICIT NONE

! Scalar Arguments
    REAL(KIND=c_double), INTENT(OUT) :: &
         f            ! Value of the objective function.
    INTEGER(KIND=c_int), INTENT(OUT) :: & 
         nit,nfe,nge  ! Number of iterations, and function and subgradient evaluations.
    INTEGER(KIND=c_int), INTENT(INOUT) :: &
         mc           ! Maximum number of stored corrections.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.
                      !  -6  - Unspecified error.


! Local Arrays
    REAL(KIND=c_double), DIMENSION(n) :: &
         xo, &        ! Previous vector of variables.
         g, &         ! Subgradient of the objective function.
         gp, &        ! Previous subgradient of the ohjective function.
         ga, &        ! Aggregate subgradient.
         s, &         ! Difference of current and previous variables.
         u, &         ! Difference of current and previous subgradients.
         d, &         ! Direction vector.
         tmpn1        ! Auxiliary array.
    REAL(KIND=c_double), DIMENSION(n*(mcu+1)) :: &
         sm,um        ! Matrises whose columns are stored differences of
                      ! variables (sm) and subgradients (um).
    REAL(KIND=c_double), DIMENSION((mcu+2)*(mcu+1)/2) :: &
         rm, &        ! pper triangular matrix stored columnwise in the 
                      ! one-dimensional array.
         umtum        ! Matrix whose columns are stored subgradient differences.
    REAL(KIND=c_double), DIMENSION(mcu+1) :: &
         cm, &        ! Diagonal matrix.
         smtgp, &     ! smtgp = trans(sm)*gp.
         umtgp, &     ! umtgp = trans(um)*gp.
         tmpmc1, &    ! Auxiliary array.
         tmpmc2       ! Auxiliary array.
    REAL(KIND=c_double), DIMENSION(n*na) :: &
         ax,ag        ! Matrix whose columns are bundle points and subgradients.
    REAL(KIND=c_double), DIMENSION(na) :: &
         af           ! Vector of bundle values.
    REAL(KIND=c_double), DIMENSION(mnma) :: &
         fold         ! Old function values.
         
! Local Scalars
    REAL(KIND=c_double) :: &
         alfn, &      ! Locality measure.
         alfv, &      ! Aggregate locality measure.
         epsr, &      ! Line search parameter.
         dnorm, &     ! Euclidean norm of the direction vector.
         gnorm, &     ! Euclidean norm of the aggregate subgradient.
         xnorm, &     ! Stopping criterion.
         pxnorm, &    ! Previous stopping criterion.
         p, &         ! Directional derivative.
         tmax, &      ! Maximum stepsize.
         t, &         ! Stepsize.
         theta, &     ! Correction parameter for stepsize.
         fo, &        ! Previous value of the objective.
         gamma, &     ! Scaling parameter.
         rmse_best    ! The best RMSE for the validation set.
    INTEGER(KIND=c_int) :: i, &
         mcinit, &    ! Initial maximum number of stored corrections.
         mcc, &       ! Current number of stored corrections.
         inew, &      ! Index for the circular arrays.
         inewnma, &   ! Index for the circular arrays containing function values.
         ibfgs, &     ! Index of the type of BFGS update.
         isr1, &      ! Index of the type of SR1 update.
         iters, &     ! Null step indicator.
                      !   0  - Null step.
                      !   1  - Serious step.
         nnk, &       ! Consecutive null steps counter.
         ibun, &      ! Index for the circular arrays in bundle updating.
         mal, &       ! Current size of the bundle.
         ic, &        ! Correction indicator.
         icn, &       ! Correction indicator for null steps.
         iflag, &     ! Index for adaptive version.
         neps, &      ! Number of consecutive equal stopping criterions.
         ntesf, &     ! Number of tests on function decrease.
         ncres, &     ! Number of restarts.
         nres, &      ! Number of consecutive restarts.
         nress, &     ! Number of consecutive restarts in case of tmax < tmin.
         nout         ! Auxilary printout specification.

! Intrinsic Functions
    INTRINSIC ABS,MAX,SQRT

! Computational Time
    REAL(KIND=c_double) strtim,elapstim
       
! Parameters
    REAL(KIND=c_double), PARAMETER :: &
         fmin    = -large, &        ! Smallest acceptable value of the function.
         tmin    = 1.0E-12_c_double, &  ! Minimum stepsize.
         lengthd = 1.0E+20_c_double, &  ! Direction vector length.
         rho     = 1.0E-12_c_double     ! Correction parameter.
    INTEGER(KIND=c_int), PARAMETER :: &
         maxeps = 20, &             ! Maximum number of consecutive equal stopping criterions.
         maxnrs = 2000              ! Maximum number of restarts.
      
      
    IF (iiprint > 3) THEN
       IF (method == 0) WRITE (6,FMT='(1X,''Entry to LMBM:'')')
       IF (method == 1) WRITE (6,FMT='(1X,''Entry to LBB:'')')
    END IF
      
     
! Initialization

    inewnma = 1
    nout   = 0
    nit    = 0
    nfe    = 0
    nge    = 0
    ntesf  = 0
    nres   = 1
    ncres  = -1
    nress  = 0
    neps   = 0
    iterm  = 0
    iters  = 1
    nnk    = 0
    isr1   = 0
    alfn   = zero
    alfv   = zero
    mcinit = mc
      
    tmax   = xmax
    xnorm  = large
    fold   = -large
    

    !epsr   = 0.25_c_double+small
    epsr   = 1.0E-05_c_double+small
    IF (epsl+epsl >= epsr) THEN
       epsr = epsl+epsl + small
       IF (epsr >= half) THEN
          CALL wprint(iterm,iiprint,-2)
       END IF
    END IF
        
! Computation of the value and the subgradient of the objective
! function and the search direction for the first iteration

    CALL myf(n,x,f,iterm,nproblem)
    CALL myg(n,x,g,iterm,nproblem)

       !WRITE(*,*)
       !WRITE(*,*) 'x:'
       !DO i =1, 31        
       !WRITE(*,*) 'x(',i,')=', x(i)      
       !END DO     
       !WRITE(*,*)  
       !WRITE(*,*) 'g='            
       !DO i =1, 31        
       !WRITE(*,*) 'g(',i,')=', g(i)      
       !END DO  
       !WRITE(*,*) 


    nfe = nfe + 1
    nge = nge + 1

    IF (iterm /= 0) GO TO 900
    
    CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
         alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
      
    CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)

     
! Start of the iteration
            
    iteration: DO
    
! Computational time

       IF (time > 0.0E+00) THEN
          CALL getime(elapstim)
          IF (elapstim-strtim > time) THEN
             iterm = 6
             EXIT iteration
          END IF
       END IF


! Computation of norms

       IF (iters > 0) THEN
          gnorm = vdot(n,g,g)
          dnorm = SQRT(vdot(n,d,d))

          p = vdot(n,g,d)
          
       ELSE
          gnorm = vdot(n,ga,ga)
          dnorm = SQRT(vdot(n,d,d))

          p = vdot(n,ga,d)
       END IF

   
! Test on descent direction

       IF (p+small*SQRT(gnorm)*dnorm <= zero) THEN
          nres = 0
          
       ELSE
          nres = nres + 1
          IF (nres == 1) THEN
             CALL wprint(iterm,iiprint,-3)
             
             CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                  alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
             IF (ncres > maxnrs) THEN
                nout = maxnrs
                iterm = -2
                EXIT iteration
             END IF
             
             CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
             
             CYCLE iteration
          END IF
          nout = -1
          iterm = -1
          EXIT iteration
       END IF
       
       
! Stopping criterion

       nit = nit + 1
       pxnorm = xnorm
       xnorm = -p + 2.0_c_double*alfv

     
! Tests for termination

       IF (xnorm <= 1.0E+03_c_double*tolg .AND. &
            (mcc > 0 .OR. ibfgs == 2)) THEN
          
          IF(half*gnorm + alfv <= tolg2 .AND. &
               xnorm <= tolg) THEN
          
             iterm = 1
             EXIT iteration
          END IF
       
          IF (mc < mcu .AND. iflag == 0) THEN
             mc = mc+1
             iflag = 1
          END IF
       END IF


       IF (nfe >= mfe) THEN
          nout = mfe
          iterm = 4
          EXIT iteration
       END IF

      
       IF (nit >= mittt) THEN
          nout = mittt
          iterm = 5
          EXIT iteration
       END IF

      
       IF (f <= tolb) THEN
          iterm = 7
          EXIT iteration
       END IF
    
      
       IF (iters == 0) THEN
          IF (ABS(xnorm - pxnorm) <= small) THEN
             neps = neps + 1
          
             IF (neps > maxeps) THEN
                iterm = -4
                EXIT iteration
             END IF

          ELSE
             neps = 0
          END IF

       ELSE
          neps = 0
       END IF


! Correction
    
       IF (-p < rho*gnorm .OR. icn == 1) THEN

          xnorm = xnorm + rho*gnorm
          dnorm = SQRT(dnorm*dnorm-2.0_c_double*rho*p+rho*rho*gnorm)
         
          IF (iters > 0) THEN
             DO i=1,n
                d(i) = d(i)-rho*g(i)
             END DO

          ELSE
             DO i=1,n
                d(i) = d(i)-rho*ga(i)
             END DO
             icn = 1
          END IF

          ic = 1
            
       ELSE
          ic = 0
       END IF


       IF (pxnorm < xnorm .AND. nnk > 2) THEN
          CALL wprint(iterm,iiprint,-4)
       END IF
      

       CALL rprint(n,nit,nfe,nge,x,f,xnorm,half*gnorm+alfv,iterm,iiprint)
      
     
!     Preparation of line search

       fo = f
       
       IF (iters > 0) THEN
          CALL copy2(n,x,xo,g,gp)
       END IF

       IF (dnorm > zero) tmax = xmax/dnorm
       
       IF (tmax > tmin) THEN
          nress = 0

       ELSE
          nress = nress + 1
          IF (nress == 1) THEN
             CALL wprint(iterm,iiprint,-5)
             
             CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                  alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
             
             IF (ncres > maxnrs) THEN
                nout = maxnrs
                iterm = -2
                EXIT iteration
             END IF
             
             CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
             CYCLE iteration
          END IF
          iterm = -1
          EXIT iteration
       END IF
      

! Initial step size

       CALL tinit(n,na,mal,x,af,ag,ax,ibun,d,f,p,t,tmax,tmin, &
            eta,iters,iterm)

       IF (iterm /= 0) EXIT iteration
     

! Line search with directional derivatives which allows null steps

       theta = one
       IF (dnorm > lengthd) THEN
          theta=lengthd/dnorm
       END IF
    
       IF (inma == 2) THEN! Line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.
      
                CALL lls(n,x,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,xnorm, &
                    theta,epsl,epsr,eta,iters,nfe,nge,nnk,iterm)
            ELSE ! Nonmonotone line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.

                CALL nmlls(n,x,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,xnorm, &
                    theta,epsl,epsr,eta,iters,inewnma,nfe,nge,nnk,iterm)

            END IF
          
       IF (iterm /= 0) EXIT iteration
       
       IF (tolf2 >= 0) THEN
          IF (ABS(fo-f) <= tolf2*small*MAX(ABS(f),ABS(fo),one) &
               .AND. iters == 1) THEN
             
             iterm = 3
             EXIT iteration
          END IF
       END IF

       IF (ABS(fo-f) <= tolf) THEN
          ntesf = ntesf + 1
          
          if (ntesf >= mtesf .AND. iters == 1) THEN
             iterm = 2
             EXIT iteration
          END IF
          
       ELSE
          ntesf = 0
       END IF
      

! Bundle updating

       CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
  

! Computation of variables difference 

       CALL xdiffy(n,x,xo,s)
  

! Computation of aggregate values and gradients difference

       IF (iters == 0) THEN
          nnk = nnk + 1

          IF (nnk == 1) THEN
             CALL copy(n,gp,tmpn1)
             CALL xdiffy(n,g,gp,u)
             
             CALL agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
                  rm,cm,umtum,alfn,alfv,gamma,ic,rho)
             
          ELSE
             IF (method == 0) THEN
                CALL copy(n,ga,tmpn1)
                CALL aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv &
                     ,umtum,rm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm &
                     ,um,icn,rho)
                CALL xdiffy(n,g,gp,u)
             ELSE
                CALL copy(n,ga,tmpn1)
                CALL xdiffy(n,g,gp,u)
                CALL agskip(n,mc,mcc,inew,iflag,g,gp,ga,d,u,alfn,alfv &
                     ,umtum,rm,cm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm &
                     ,um,icn,rho)
             END IF
          END IF
          
          CALL copy(n,xo,x)
          f = fo
          
       ELSE
          IF (nnk /= 0) THEN
             CALL copy(n,ga,tmpn1)
          ELSE
             CALL copy(n,gp,tmpn1)
          END IF
          nnk = 0
          CALL xdiffy(n,g,gp,u)
       END IF

     
! Serious step initialization

       IF (iters > 0) THEN
          icn = 0
          alfn = zero
          alfv = zero
       END IF

      
! Direction finding
    
       IF (iters > 0) THEN
         
     
! BFGS update and direction determination

          CALL dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
               umtum,cm,smtgp,umtgp,gamma,tmpn1,method,iscale)

       ELSE
          IF (method == 0) THEN
            
     
! SR1 update and direction determination
             
             CALL dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm, &
                  umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iiprint)
             ibfgs=0
          ELSE
             
     
! BFGS skipping and direction determination

             CALL dlskip(n,mc,mcc,inew,ibfgs,iflag,d,ga,sm,um,rm, &
                  umtum,cm,tmpmc1,tmpmc2,gamma,iscale)
          END IF
       END IF
    
    END DO iteration
      

900 CONTINUE
        

! Printout the final results

    IF (iiprint > 3) THEN
       IF (method == 0) WRITE (6,FMT='(1X,''Exit from LMBM:'')')
       IF (method == 1) WRITE (6,FMT='(1X,''Exit from LBB:'')')
    END IF
    CALL wprint(iterm,iiprint,nout)
    CALL rprint(n,nit,nfe,nge,x,f,xnorm,half*gnorm+alfv,iterm,iiprint)

  CONTAINS

    SUBROUTINE restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
         alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)  ! Initialization
      
!      USE param, ONLY : zero,one  ! given in host
!      USE lmbm_sub, ONLY : vneg,copy  ! given in host
      IMPLICIT NONE

! Array Arguments
      REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
           gp        ! Basic subgradient of the objective function.
      REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
           g         ! Current (auxiliary) subgradient of the objective function.
      REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
           d         ! Search direction.

! Scalar Arguments
      INTEGER, INTENT(IN) :: & 
           n, &      ! Number of variables
           mcinit    ! Initial maximum number of stored corrections.
      INTEGER, INTENT(OUT) :: &
           mc, &     ! Current maximum number of stored corrections.
           mcc, &    ! Current number of stored corrections.
           inew, &   ! Index for the circular arrays.
           ibun, &   ! Index for the circular arrays in bundle updating.
           ibfgs, &  ! Index of the type of BFGS update.
           nnk, &    ! Consecutive null steps counter.
           ic, &     ! Correction indicator.
           icn, &    ! Correction indicator for null steps.
           mal, &    ! Current size of the bundle.
           iflag     ! Index for adaptive version.
      INTEGER, INTENT(INOUT) :: &
           iters, &  ! Null step indicator.
                     !   0  - Null step.
                     !   1  - Serious step.
           ncres     ! Number of restarts.
      REAL(KIND=c_double), INTENT(OUT) :: & 
           alfn, &   ! Locality measure.
           alfv, &   ! Aggregate locality measure.
           gamma     ! Scaling parameter.


! Restart
      mc    = mcinit
      mcc   = 0
      inew  = 1
      ibun  = 1
      ibfgs = 0
      ic    = 0
      icn   = 0
      mal   = 0
      ncres = ncres + 1
      iflag = 0
        
      IF (iters == 0) THEN
         CALL copy(n,gp,g)
         iters = 1
         nnk = 0
         alfv=zero
         alfn=zero
      END IF

      gamma = one
      CALL vneg(n,g,d)
      
    END SUBROUTINE restar
    
      
    SUBROUTINE dobun(n,ma,mal,x,g,f,ay,ag,af,iters,ibun)  ! Bundle construction

!      USE lmbm_sub, ONLY : copy2 ! given in host
      IMPLICIT NONE

! Array Arguments
      REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables
           g         ! Subgradient of the objective function.
      REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
           ay, &     ! Matrix whose columns are bundle points 
                     ! (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients 
                     ! (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions 
                     ! (stored in one-dimensional ma -array).

! Scalar Arguments
      INTEGER(KIND=c_int), INTENT(IN) :: & 
           n, &      ! Number of variables
           iters, &  ! Null step indicator.
                     !   0  - Null step.
                     !   1  - Serious step.
           ma        ! Maximum size of the bundle.
      INTEGER(KIND=c_int), INTENT(INOUT) :: &
           ibun, &   ! Index for the circular arrays in bundle updating.
           mal       ! Current size of the bundle.
      REAL(KIND=c_double), INTENT(IN) :: & 
           f         ! Value of the objective function.

! Local Scalars
      INTEGER(KIND=c_int) :: i,j
      
      IF (iters == 1) THEN
     
! Serious step

         af(ibun) = f
         i = (ibun-1)*n + 1
         CALL copy2(n,g,ag(i:),x,ay(i:))
         
      ELSE

! Null step

         IF (mal < ma) THEN
            
            af(ibun) = af(mal)
            af(mal) = f
            
            i = mal*n + 1
            CALL copy2(n,ag(i-n:),ag(i:),ay(i-n:),ay(i:))
            CALL copy2(n,g,ag(i-n:),x,ay(i-n:))
            
         ELSE
            i = ibun-1
            IF (i < 1) i = mal
            af(ibun) = af(i)
            af(i) = f
            
            i = (ibun-2)*n + 1
            IF (i < 1) i = (mal-1)*n + 1
            j = (ibun-1)*n + 1
            CALL copy2(n,ag(i:),ag(j:),ay(i:),ay(j:))
            CALL copy2(n,g,ag(i:),x,ay(i:))
         END IF
         
      END IF
      
      mal = mal + 1
      IF (mal > ma) mal = ma
      
      ibun = ibun + 1
      IF (ibun > ma) ibun = 1
      
    END SUBROUTINE dobun
    
  END SUBROUTINE lmbm

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE tinit *                                             *
!*                                                                      *
!*     Initial stepsize selection for limited memory bundle method      *
!*                                                                      *
!************************************************************************
  SUBROUTINE tinit(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,tmax,tmin,eta,iters,iterm)

    USE param, ONLY : zero,half,one,large
!    USE param, ONLY : zero,one  ! these are the one needed in tinit itself 
                                 ! half and large are used in destep and nulstep
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x, &      ! Vector of variables (n array).
         d, &      ! Direction vector (n array).
         ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
         ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Scalar Arguments
    REAL(KIND=c_double), INTENT(OUT) :: & 
         t         ! Initial stepsize
    REAL(KIND=c_double), INTENT(IN) :: & 
         p, &      ! Directional derivative.
         eta, &    ! Distance measure parameter.
         f, &      ! Value of the objective function.
         tmax, &   ! Upper limit for stepsize parameter.
         tmin      ! Lower limit for stepsize parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         na, &     ! Maximum size of the bundle.
         mal, &    ! Current size of the bundle.
         ibun, &   ! Index for the circular arrays in bundle updating.
         iters     ! Null step indicator.
                   !   0  - Null step.
                   !   1  - Serious step.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
         iterm     ! Cause of termination:
                   !   0  - Everything is ok.
                   !  -6  - Error.


! Intrinsic Functions
    INTRINSIC MAX,MIN

    t = MIN(one,tmax)

    IF (p == zero) RETURN

    IF (iters == 1) THEN
       CALL destep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    ELSE
       CALL nulstep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    END IF

    t = MIN(MAX(t,tmin),tmax)

  CONTAINS

    SUBROUTINE destep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

!      USE param, ONLY : zero,half,one,large  ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=c_double), INTENT(INOUT) :: & 
           t         ! Initial stepsize
      REAL(KIND=c_double), INTENT(IN) :: & 
           df, &     ! Directional derivative.
           eta, &    ! Distance measure parameter.
           f         ! Value of the objective function.
      INTEGER(KIND=c_int), INTENT(IN) :: & 
           n, &      ! Number of variables
           ma, &     ! Maximum size of the bundle.
           mal, &    ! Current size of the bundle.
           ibun      ! Index for the circular arrays in bundle updating.
      INTEGER(KIND=c_int), INTENT(OUT) :: &
           iterm     ! Cause of termination:
                     !   0  - Everything is ok.
                     !  -6  - Error.

! Array Arguments
      REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables (n array).
           d, &      ! Direction vector (n array).
           ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Local Arrays
      REAL(KIND=c_double), DIMENSION(2*ma) :: &
           tmparray  ! Auxiliary array.

! Local Scalars
      REAL(KIND=c_double) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
      INTEGER(KIND=c_int) :: i,j,jn,k,l,lq,ib

! Intrinsic Functions
      INTRINSIC ABS,REAL,MAX,MIN,SQRT

      iterm = 0
      alfl = zero
      betl = zero
      
      w = df*t* (one - t*half)
      
     
! Initial choice of possibly active lines
      
      k = 0
      l = -1
      jn = (ibun-1)*n
      betr = - large
      DO j=1,mal-1
         ib = ibun - 1 + j
         IF (ib > mal) ib = ib - mal
         IF (jn >= mal*n) jn = jn - mal*n
         r = zero
         bet = zero
         alfl = af(ib) - f
         DO i=1,n
            dx = x(i) - ay(jn+i)
            q = ag(jn+i)
            r = r + dx*dx
            alfl = alfl + dx*q
            bet = bet + d(i)*q
         END DO
         alf = MAX(ABS(alfl),eta*r)
         r = one - bet/df
         IF (r*r + (alf+alf)/df > 1.0E-6_c_double) THEN
            k = k + 1
            tmparray(k) = alf
            tmparray(ma+k) = bet
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = k
            END IF
         END IF
         
         betr = MAX(betr,bet-alf)
         jn = jn + n
      END DO
      lq = -1
      IF (betr <= df*half) RETURN
      lq = 1
      IF (l <= 0) RETURN
      betr = tmparray(ma+l)
      IF (betr <= zero) THEN
         IF (t < one .OR. betr == zero) RETURN
         lq = 2
      END IF
      
      alfr = tmparray(l)


! Iteration loop

      ds_iteration: DO 
         
         IF (lq >= 1) THEN
            q = one - betr/df
            r = q + SQRT(q*q + (alfr+alfr)/df)
            IF (betr >= zero) r = - (alfr+alfr)/ (df*r)
            r = MIN(1.95_c_double,MAX(zero,r))
         ELSE
            IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_c_double*df) RETURN
            IF (betr-betl  == zero) THEN
               iterm = -6
               RETURN
            END IF
            r = (alfr-alfl)/ (betr-betl)
         END IF
       
         IF (ABS(t-r) < 1.0E-4_c_double) RETURN
         t = r
         tmparray(l) = - one
         w = t*betr - alfr
         l = -1
         DO j = 1,k
            alf = tmparray(j)
            IF (alf < zero) EXIT
            bet = tmparray(ma+j)
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = j
            END IF
         END DO
         IF (l < 0) RETURN
       
         bet = tmparray(ma+l)
         IF (bet == zero) RETURN
     

!     New interval selection

         alf = tmparray(l)
         IF (bet < zero) THEN
            IF (lq == 2) THEN
               alfr = alf
               betr = bet
               
            ELSE
               alfl = alf
               betl = bet
               lq = 0
            END IF
       
         ELSE
            IF (lq == 2) THEN
               alfl = alfr
               betl = betr
               lq = 0
            END IF
       
            alfr = alf
            betr = bet
         END IF
       
      END DO ds_iteration
      
    END SUBROUTINE destep

     
    SUBROUTINE nulstep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

!      USE param, ONLY : zero,half,one,large  ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=c_double), INTENT(INOUT) :: & 
           t         ! Initial stepsize
      REAL(KIND=c_double), INTENT(IN) :: & 
           df, &     ! Directional derivative.
           eta, &    ! Distance measure parameter.
           f         ! Value of the objective function.
      INTEGER(KIND=c_int), INTENT(IN) :: & 
           n, &      ! Number of variables
           ma, &     ! Maximum size of the bundle.
           mal, &    ! Current size of the bundle.
           ibun      ! Index for the circular arrays in bundle updating.
      INTEGER(KIND=c_int), INTENT(OUT) :: &
           iterm     ! Cause of termination:
                     !   0  - Everything is ok.
                     !  -6  - Error.

! Array Arguments
      REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables (n array).
           d, &      ! Direction vector (n array).
           ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions (stored in one-dimensional 4*ma -array).

! Local Arrays
      REAL(KIND=c_double), DIMENSION(2*ma) :: &
           tmparray  ! Auxiliary array.

! Local Scalars
      REAL(KIND=c_double) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
      INTEGER(KIND=c_int) :: i,j,jn,k,l,ib

! Intrinsic Functions
      INTRINSIC ABS,REAL,MAX,MIN,SQRT

      iterm = 0
      w = df*t
     
! Initial choice of possibly active parabolas

      k = 0
      l = -1
      jn = (ibun-1)*n
      betr = - large
      DO j = 1,mal - 1
         ib = ibun - 1 + j
         IF (ib > mal) ib = ib - mal
         IF (jn >= mal*n) jn = jn - mal*n
         bet = zero
         r = zero
         alfl = af(ib) - f
         DO i = 1,n
            dx = x(i) - ay(jn+i)
            r = r + dx*dx
            q = ag(jn+i)
            alfl = alfl + dx*q
            bet = bet + d(i)*q
         END DO
         alf = MAX(ABS(alfl),eta*r)
         betr = MAX(betr,bet-alf)
         IF (alf < bet-df) THEN
            k = k + 1
            r = t*bet - alf
            tmparray(k) = alf
            tmparray(ma+k) = bet
            IF (r > w) THEN
               w = r
               l = k
            END IF
         END IF
         jn = jn + n
      END DO
      IF (l <= 0) RETURN
      betr = tmparray(ma+l)
      alfr = tmparray(l)
      alf = alfr
      bet = betr
      alfl = zero
      betl = df
    
     
!     Iteration loop
    
      ns_iteration: DO

         w = bet/df
         IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_c_double*df) RETURN
         IF (betr-betl  == zero) THEN
            iterm = -6
            RETURN
         END IF
         r = (alfr-alfl)/ (betr-betl)
         IF (ABS(t-w) < ABS(t-r)) r = w
         q = t
         t = r
         IF (ABS(t-q) < 1.0E-3_c_double) RETURN
         tmparray(l) = - one
         w = t*bet - alf
         l = -1
         DO j=1,k
            alf = tmparray(j)
            IF (alf < zero) EXIT
            bet = tmparray(ma+j)
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = j
            END IF
         END DO

         IF (l <= 0) RETURN
         bet = tmparray(ma+l)
         q = bet - t*df
         IF (Q == zero) RETURN

     
!     New interval selection

         alf = tmparray(l)
         IF (q < zero) THEN
            alfl = alf
            betl = bet
         ELSE
            alfr = alf
            betr = bet
         END IF

      END DO ns_iteration

    END SUBROUTINE nulstep

  END SUBROUTINE tinit
      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE lls *                                               *
!*                                                                      *
!*     Special line search for limited memory bundle method             *
!*                                                                      *
!************************************************************************

  SUBROUTINE lls(n,x,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
       eta,iters,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
         myf,myg   ! Computation of the value and the subgradient of
                   ! the objective function. 
    USE lmbm_sub, ONLY : &
         scsum, &  ! Sum of a vector and the scaled vector.
         vdot      ! Dot product of vectors.
    USE initializat, ONLY : &
        nproblem, &      ! the solved problem
        maxint => maxnin ! Maximum number of interpolations.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         xo        ! Previous vector of variables.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         x         ! Vector of variables.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         g         ! Subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         epsl,&    ! Linesearch parameter.
         epsr,&    ! Linesearch parameter.
         t, &      ! Stepsize
         p         ! Directional derivative.
    REAL(KIND=c_double), INTENT(IN) :: & 
         theta, &  ! Scaling parameter.
         eta, &    ! Distance measure parameter.
         fo, &     ! Previous value of the objective function.
         wk, &     ! Stopping parameter.
         dnorm, &  ! Euclidean norm of the direction vector.
         tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=c_double), INTENT(OUT) :: & 
         f, &      ! Value of the objective function.
         alfn      ! Locality measure.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         nnk       ! Number of consequtive null steps.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         nfe, &    ! Number of function evaluations.
         nge       ! Number of subgradient evaluations.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
         iters, &  ! Null step indicator.
                   !   0  - Null step.
                   !   1  - Serious step.
         iterm     ! Cause of termination:
                   !   0  - Everything is ok.
                   !  -3  - Failure in function or subgradient calculations.

! Local Scalars
    INTEGER(KIND=c_int) :: nin ! Number of interpolations.
    REAL(KIND=c_double) :: &
         tl,tu, &  ! Lower and upper limits for t used in interpolation.
         fl,fu, &  ! Values of the objective function for t=tl and t=tu.
         epsa,epst,epslk,epsrk, & ! Line search parameters.
         thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

! Parameters
!    INTEGER, PARAMETER :: maxint = 200  ! Maximum number of interpolations.

! Intrinsic Functions
    INTRINSIC ABS,MAX
    INTEGER(KIND=c_int) :: i
      

! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo

    IF (theta < one) THEN
       epst  = theta*epst
       epsa  = theta*epsa
       epslk = epsl
       epsl  = theta*epsl
       epsrk = epsr
       epsr  = theta*epsr
    END IF
 
    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk


! Function evaluation at a new point

    lls_iteration: DO
       
       CALL scsum(n,theta*t,d,xo,x)

       CALL myf(n,x,f,iterm,nproblem)
       nfe = nfe + 1
       IF (iterm /= 0) RETURN

    
! Null/descent step test (ITERS=0/1)

       iters = 1
       IF (f <= fo - t*epstwk) THEN
          tl = t
          fl = f
       ELSE
          tu = t
          fu = f
       END IF
      

! Additional interpolation

       IF (f > fo .AND. tu-tl >= tmin*0.1_c_double &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

          nin=nin+1
          IF (tl == zero .AND. wk > zero) THEN
             t = qint(tu,fl,fu,wk,one-half/(one-epst))
          ELSE
             t = half*(tu+tl)
          END IF
          CYCLE lls_iteration
       END IF

       CALL myg(n,x,g,iterm,nproblem)
       
       !WRITE(*,*)
       !WRITE(*,*) 'x:'
       !DO i =1, 31        
       !WRITE(*,*) 'x(',i,')=', x(i)      
       !WRITE(*,*) 'x(',i,')=', x(i)      
       !END DO     
       !WRITE(*,*)  
       !WRITE(*,*) 'g='            
       !DO i =1, 31        
       !WRITE(*,*) 'g(',i,')=', g(i)      
       !END DO  
       !WRITE(*,*) 
       nge = nge + 1
       IF (iterm /= 0) RETURN

       p = theta*vdot(n,g,d)
       alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


! Serious step

       IF (f <= fo - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


! Null step

       IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_c_double .OR. &
            nin >= maxint) THEN
          ITERS = 0
          EXIT lls_iteration
       END IF


! Interpolation

       nin=nin+1
       IF (tl == zero .AND. wk > zero) THEN
          t = qint(tu,fl,fu,wk,one-half/(one-epst))
       ELSE
          t = half*(tu+tl)
       END IF

    END DO lls_iteration

    IF (theta /= one) THEN
       epsl = epslk
       epsr = epsrk
    END IF

  CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

!      USE param, ONLY : half,one  ! given in host
      IMPLICIT NONE
      
! Scalar Arguments
      REAL(KIND=c_double), INTENT(IN) :: & 
           fl, &  ! Value of the objective function.
           fu, &  ! Value of the objective function for t=tu.
           wk, &  ! Directional derivative.
           tu, &  ! Upper value of the stepsize parameter.
           kappa  ! Interpolation parameter.
      REAL(KIND=c_double) :: & 
           t      ! Stepsize.

! Local Scalars
      REAL(KIND=c_double) :: tmp1,tmp2

! Intrinsic Functions
      INTRINSIC MAX


      tmp1 = (fu-fl)/ (-wk*tu)

     
! Quadratic interpolation with one directional derivative

      tmp2 = 2.0_c_double * (one - tmp1)

      IF (tmp2 > one) THEN
   
! Interpolation accepted
       
         t = MAX(kappa*tu,tu/tmp2)
         RETURN
      END IF
      
     
! Bisection
    
      t = half*tu
      
    END FUNCTION qint

  END SUBROUTINE lls

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE nmlls *                                             *
!*                                                                      *
!*     Nonmonotonic weak Wolfe-type line search                         *
!*                                                                      *
!************************************************************************

SUBROUTINE nmlls(n,x,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
    eta,iters,inewnma,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
        myf,myg    ! Computation of the value and the subgradient of
                   ! the objective function.
    USE lmbm_sub, ONLY : &
        scsum      ! Sum of a vector and the scaled vector.
    USE initializat, ONLY : &
        nproblem, &  ! The solved problem
        mnma, &      ! Maximum number of function values used in line search.
        maxint => maxnin   ! Maximum number of interpolations.

    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        xo        ! Previous vector of variables.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
        x         ! Vector of variables.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
        g         ! Subgradient of the objective function.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
        fold           ! Old function values.

    ! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: &
        epsl,&    ! Linesearch parameter.
        epsr,&    ! Linesearch parameter.
        t, &      ! Stepsize
        p         ! Directional derivative.
    REAL(KIND=c_double), INTENT(IN) :: &
        theta, &  ! Scaling parameter.
        eta, &    ! Distance measure parameter.
        fo, &     ! Previous value of the objective function.
        wk, &     ! Stopping parameter.
        dnorm, &  ! Euclidean norm of the direction vector.
        tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=c_double), INTENT(OUT) :: &
        f, &      ! Value of the objective function.
        alfn      ! Locality measure.
    INTEGER(KIND=c_int), INTENT(IN) :: &
        n, &      ! Number of variables
        nnk       ! Number of consequtive null steps.
    INTEGER(KIND=c_int), INTENT(INOUT) :: &
        inewnma, &! index for array.
        nfe, &    ! Number of function evaluations.
        nge       ! Number of subgradient evaluations.
    INTEGER(KIND=c_int), INTENT(OUT) :: &
        iters, &  ! Null step indicator.
                  !   0  - Null step.
                  !   1  - Serious step.
        iterm     ! Cause of termination:
                  !   0  - Everything is ok.
                  !  -3  - Failure in function or subgradient calculations.

    ! Local Scalars
    INTEGER(KIND=c_int) :: nin ! Number of interpolations.
    REAL(KIND=c_double) :: &
        tl,tu, &  ! Lower and upper limits for t used in interpolation.
        fl,fu, &  ! Values of the objective function for t=tl and t=tu.
        ftmp, &   ! Maximum function value from mnma last iterations.
        epsa,epst,epslk,epsrk, & ! Line search parameters.
        thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

    ! Intrinsic Functions
    INTRINSIC ABS,MAX,MAXVAL,DOT_PRODUCT



    ! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo

    IF (theta < one) THEN
        epst  = theta*epst
        epsa  = theta*epsa
        epslk = epsl
        epsl  = theta*epsl
        epsrk = epsr
        epsr  = theta*epsr
    END IF

    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk

    ! Updating circular array fold

    IF (iters==1) THEN
        fold(inewnma) = fo
        inewnma = inewnma + 1
        IF (inewnma > mnma) inewnma = 1
    END IF
    ftmp = MAXVAL(fold)


    ! Function evaluation at a new point

    lls_iteration: DO

        CALL scsum(n,theta*t,d,xo,x)

        CALL myf(n,x,f,iterm,nproblem)
        nfe = nfe + 1
        IF (iterm /= 0) RETURN


        ! Null/descent step test (ITERS=0/1)

        iters = 1
        IF (f <= ftmp - t*epstwk) THEN
            tl = t
            fl = f
        ELSE
            tu = t
            fu = f
        END IF


        ! Additional interpolation

        IF (f > ftmp .AND. tu-tl >= tmin*0.1_c_double &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

            nin=nin+1
            IF (tl == zero .AND. wk > zero) THEN
                t = qint(tu,fl,fu,wk,one-half/(one-epst))
            ELSE
                t = half*(tu+tl)
            END IF
            CYCLE lls_iteration
        END IF

        CALL myg(n,x,g,iterm,nproblem)
        nge = nge + 1
        IF (iterm /= 0) RETURN

        p = theta*DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


        ! Serious step

        IF (f <= ftmp - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


        ! Null step

        IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_c_double .OR. &
            nin >= maxint) THEN
            ITERS = 0
            EXIT lls_iteration
        END IF


        ! Interpolation

        nin=nin+1
        IF (tl == zero .AND. wk > zero) THEN
            t = qint(tu,fl,fu,wk,one-half/(one-epst))
        ELSE
            t = half*(tu+tl)
        END IF

    END DO lls_iteration

    IF (theta /= one) THEN
        epsl = epslk
        epsr = epsrk
    END IF

CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

        !      USE param, ONLY : half,one  ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=c_double), INTENT(IN) :: &
            fl, &  ! Value of the objective function.
            fu, &  ! Value of the objective function for t=tu.
            wk, &  ! Directional derivative.
            tu, &  ! Upper value of the stepsize parameter.
            kappa  ! Interpolation parameter.
        REAL(KIND=c_double) :: &
            t      ! Stepsize.

        ! Local Scalars
        REAL(KIND=c_double) :: tmp1,tmp2

        ! Intrinsic Functions
        INTRINSIC MAX


        tmp1 = (fu-fl)/ (-wk*tu)


        ! Quadratic interpolation with one directional derivative

        tmp2 = 2.0_c_double * (one - tmp1)

        IF (tmp2 > one) THEN

            ! Interpolation accepted

            t = MAX(kappa*tu,tu/tmp2)
            RETURN
        END IF


        ! Bisection

        t = half*tu

    END FUNCTION qint

END SUBROUTINE nmlls



!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlbfgs *                                            *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*g  *
!*     by the limited memory BFGS update.                               *
!*                                                                      *
!************************************************************************
    
  SUBROUTINE dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
       umtum,cm,smtgp,umtgp,gamma,tmpn1,method,iscale)
      
    USE param, ONLY : zero,small,one,half ! half is used at sclpar
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.
         vneg, &   ! Copying of a vector with change of the sign.
         xdiffy, & ! Difference of two vectors.
         xsumy, &  ! Sum of two vectors.
         scdiff, & ! Difference of the scaled vector and a vector.
         scsum, &  ! Sum of a vector and the scaled vector.
         vxdiag, & ! Multiplication of a vector and a diagonal matrix.
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         copy2     ! Copying of two vectors.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         g, &      ! Current subgradient of the objective function.
         gp        ! Previous subgradient of the objective function.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         s, &      ! Difference of current and previous variables.
         u, &      ! Difference of current and previous subgradients.
         tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         method, & ! Selection of the method:
                   !   0  - Limited memory bundle method.
                   !   1  - L-BFGS bundle method.
         iscale    ! Selection of the scaling:
                   !   0  - Scaling at every iteration with STU/UTU.
                   !   1  - Scaling at every iteration with STS/STU.
                   !   2  - Interval scaling with STU/UTU.
                   !   3  - Interval scaling with STS/STU.
                   !   4  - Preliminary scaling with STU/UTU.
                   !   5  - Preliminary scaling with STS/STU.
                   !   6  - No scaling.      
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER(KIND=c_int), INTENT(OUT) :: & 
         ibfgs     ! Index of the type of BFGS update:
                   !   1  - BFGS update: the corrections are stored.
                   !   2  - BFGS update: the corrections are not stored.
                   !   3  - BFGS update is skipped.
      
! Local Arrays
    REAL(KIND=c_double), DIMENSION(mcc+1) :: &
         tmpmc1,tmpmc2,tmpmc3,tmpmc4

! Local Scalars
    REAL(KIND=c_double) :: &
         stu, &    ! stu = trans(s)*u. 
         sts       ! sts = trans(s)*s. 
    INTEGER(KIND=c_int) :: i,j,k, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         iflag2, & ! Index for adaptive version.
         ierr      ! Error indicator

! Intrinsic Functions
    INTRINSIC SQRT,MIN,MAX


    ierr = 0
    ibfgs = 0
    iflag2 = 0
    stu = vdot(n,s,u)
    sts = vdot(n,s,s)


! Positive definiteness

    IF (stu > zero) THEN
       IF (-vdot(n,d,u)-vdot(n,tmpn1,s) < -small .OR. method == 1) THEN
          
     
! Update matrices
         
          ibfgs = 1

! Initialization of indices.

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
            
     
! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
           
     
! Computation of trans(sm)*g and trans(um)*g

          IF (inew >= mcnew) THEN
             CALL rwaxv2(n,mcnew,sm((inew-mcnew)*n+1:),&
                  um((inew-mcnew)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
          ELSE
             CALL rwaxv2(n,inew,sm,um,g,g,tmpmc1,tmpmc2)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),&
                  um((iold-1)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
          END IF
            
            
! Computation of trans(sm)*u and trans(um)*u

          IF (inew >= mcnew) THEN
             DO i=iold,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          ELSE
             DO i=1,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
             DO i=iold,mcnew+1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          END IF
          tmpmc3(inew) = tmpmc1(inew) - vdot(n,s,gp)
          smtgp(inew)  = tmpmc1(inew)
          tmpmc4(inew) = tmpmc2(inew) - vdot(n,u,gp)
          umtgp(inew)  = tmpmc2(inew)
            
         
! Update rm and umtum

          IF (mcc >= mc .AND. iflag2 /= 1) THEN
             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO
          END IF
                      
          IF (inew >= mcnew) THEN
             CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                  tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:)&
                  ,tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
             CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),&
                  tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF
            

! Update CM

          cm(inew) = stu
            
! Computation of gamma

          gamma = sclpar(mcc,iscale,method,sts,stu,tmpmc4(inew))
            
          inew = inew + 1
          IF (inew > mc + 1) inew = 1
          IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
          
       ELSE

            
! BFGS update, corrections are not saved.
     
          ibfgs = 2

! Initialization of indices.

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
          
! Update sm and um
          
          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
            
! Computation of trans(sm)*g and trans(um)*g

          CALL rwaxv2(n,mcnew,sm,um,g,g,tmpmc1,tmpmc2)

! Computation of trans(sm)*u and trans(um)*u
          
          IF (iold /= 1) THEN
             DO i=1,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
             DO i=iold,mcnew
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          ELSE
             DO i=1,mcnew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          END IF
          tmpmc3(inew) = tmpmc1(inew) - vdot(n,s,gp)
          smtgp(inew)  = tmpmc1(inew)
          tmpmc4(inew) = tmpmc2(inew) - vdot(n,u,gp)
          umtgp(inew)  = tmpmc2(inew)


! Update rm and umtum

          IF (iold /= 1) THEN
             CALL copy2(mcnew-inew,tmpmc3(iold:),&
                  rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                  umtum((mcnew-1)*mcnew/2+1:))
             CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                  umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          ELSE
             CALL copy2(mcnew,tmpmc3,rm((mcnew-1)*mcnew/2+1:)&
                  ,tmpmc4,umtum((mcnew-1)*mcnew/2+1:))
          END IF
            

! Update cm

          cm(inew) = stu
            

! Computation of gamma

          gamma = sclpar(mcc,iscale,method,sts,stu,tmpmc4(inew))
               
       END IF
    ELSE
         
     
! BFGS update is skipped

       ibfgs = 3

       IF (mcc == 0) THEN
          iflag = 0
          CALL vneg(n,g,d)
          RETURN
       END IF
         

!     Initialization of indices.

       CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)


!     Computation of gamma

       IF (iscale >= 4) gamma = one
               
         
!     Computation of trans(sm)*g and trans(um)*g and the two intermediate values

       IF (iold <= 2) THEN
          CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),g,g,&
               smtgp(iold:),umtgp(iold:))
       ELSE
          CALL rwaxv2(n,inew-1,sm,um,g,g,smtgp,umtgp)
          CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
               um((iold-1)*n+1:),g,g,smtgp(iold:),umtgp(iold:))
       END IF
    END IF


! Computation of two intermediate values tmpmc1 and tmpmc2

    IF (iold == 1 .OR. ibfgs == 2) THEN
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mcnew,cm,tmpmc1,tmpmc2)
       CALL scsum(mcnew,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum(mcnew,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE IF (iflag == 0) THEN
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mc+1,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mc+1,cm,tmpmc1,tmpmc2)
       CALL scsum(mc+1,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum(mc+1,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE
       CALL trlieq(mcnew,mc,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mc,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mc,cm,tmpmc1,tmpmc2)
       CALL scsum(mc,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum(mc,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mc,iold,rm,tmpmc2,tmpmc3,0,ierr)
    END IF

      
! Computation of the search direction d

    IF (iold == 1 .OR. ibfgs == 2) THEN
       CALL cwmaxv(n,mcnew,um,tmpmc1,d)
       CALL xdiffy(n,d,g,d)
       CALL cwmaxv(n,mcnew,sm,tmpmc2,tmpn1)
       CALL scdiff(n,gamma,d,tmpn1,d)
    ELSE 
       CALL cwmaxv(n,inew-1,um,tmpmc1,d)
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc1(iold:),tmpn1)
       CALL xsumy(n,d,tmpn1,d)
       CALL xdiffy(n,d,g,d)
       CALL cwmaxv(n,inew-1,sm,tmpmc2,tmpn1)
       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc2(iold:),tmpn1)
       CALL xdiffy(n,d,tmpn1,d)
    END IF

  CONTAINS

    FUNCTION sclpar(mcc,iscale,method,sts,stu,utu) RESULT(spar) ! Calculation of the scaling parameter.
      
!      USE param, ONLY : small,one,half ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=c_double), INTENT(IN) :: & 
           sts, &    ! sts = trans(s)*s. 
           stu, &    ! stu = trans(s)*u. 
           utu       ! utu = trans(u)*u. 
      REAL(KIND=c_double) :: & 
           spar      ! Scaling parameter.
      INTEGER, INTENT(IN) :: & 
           mcc, &    ! Current number of stored corrections.
           method, & ! Selection of the method:
                     !   0  - Limited memory bundle method.
                     !   1  - L-BFGS bundle method.
           iscale    ! Selection of the scaling:
                     !   0  - Scaling at every iteration with STU/UTU.
                     !   1  - Scaling at every iteration with STS/STU.
                     !   2  - Interval scaling with STU/UTU.
                     !   3  - Interval scaling with STS/STU.
                     !   4  - Preliminary scaling with STU/UTU.
                     !   5  - Preliminary scaling with STS/STU.
                     !   6  - No scaling.      

! Intrinsic Functions
      INTRINSIC SQRT


! Computation of scaling parameter.

      SELECT CASE(iscale)
     
! Scaling parameter = STU/UTU

      CASE(0,2,4)
         IF (utu < SQRT(small)) THEN
            spar = one
            RETURN
         ELSE
            spar = stu/utu
         END IF
    
! Scaling parameter = STS/STU

      CASE(1,3,5)
         IF (stu < SQRT(small)) THEN
            spar = one
            RETURN
         ELSE
            spar = sts/stu
         END IF

! No scaling

      CASE DEFAULT
         spar = one
         RETURN
      END SELECT

               
!     Scaling
            
      IF (MCC == 0) THEN
         IF (spar < 0.01_c_double) spar=0.01_c_double
         IF (spar > 100.0_c_double) spar=100.0_c_double

      ELSE 

         SELECT CASE(iscale)

! Interval scaling
         CASE(2)
            IF (method == 0) THEN
               IF (spar < 0.6_c_double .OR. spar > 6.0_c_double) spar = one
            ELSE
               IF (spar < 0.01_c_double .OR. spar > 100.0_c_double) spar = one
            END IF

         CASE(3)
            IF (spar < half .OR. spar > 5.0_c_double) spar = one
               
! Preliminary scaling
         CASE(4,5)
            spar = one

! Scaling at every iteration
         CASE DEFAULT
            CONTINUE
         END SELECT

      END IF

      IF (spar < 1.0E+03_c_double*small) spar = 1.0E+03_c_double*small
         
    END FUNCTION sclpar

  END SUBROUTINE dlbfgs

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlskip *                                            *
!*                                                                      *
!*     Matrix skipping and computation of the search direction          *
!*     d = -dm*ga by the limited memory BFGS update.                    *
!*                                                                      *
!************************************************************************

  SUBROUTINE dlskip(n,mc,mcc,inew,ibfgs,iflag,d,ga,sm,um,rm,&
       umtum,cm,tmpmc1,tmpmc2,gamma,iscale)
   
    USE param, ONLY : zero,small,one
    USE lmbm_sub, ONLY : &
         vneg, &   ! Copying of a vector with change of the sign.
         xsumy, &  ! Sum of two vectors.
         xdiffy, & ! Difference of two vectors.
         scsum, &  ! Sum of a vector and the scaled vector.
         scdiff, & ! Difference of the scaled vector and a vector.
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         vxdiag    ! Multiplication of a vector and a diagonal matrix.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         ga, &     ! Current aggregate subgradient of the objective function.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum     ! Matrix umtum = trans(um) * um.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         tmpmc1, & ! Auxiliary array. On output: trans(sm)*ga.
         tmpmc2, & ! Auxiliary array. On output: trans(um)*ga.
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc        ! Declared number of stored corrections.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag, &  ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
         iscale    ! Selection of the scaling:
                   !   0   - Scaling at every iteration with STU/UTU.
                   !   1  - Scaling at every iteration with STS/STU.
                   !   2  - Interval scaling with STU/UTU.
                   !   3  - Interval scaling with STS/STU.
                   !   4  - Preliminary scaling with STU/UTU.
                   !   5  - Preliminary scaling with STS/STU.
                   !   6  - No scaling.      
    INTEGER(KIND=c_int), INTENT(OUT) :: & 
         ibfgs     ! Index of the type of L-SR1 update:
                   !   3  - BFGS update is skipped.
      
! Local Arrays
    REAL(KIND=c_double), DIMENSION(n) :: tmpn1
    REAL(KIND=c_double), DIMENSION(mcc+1) :: tmpmc3,tmpmc4,tmpmc5
 
! Local Scalars
    INTEGER(KIND=c_int) :: &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr,&    ! Error indicator.
         iflag2    ! Index for adaptive version.


! BFGS update is skipped

    ierr = 0
    ibfgs = 3

    IF (mcc == 0) THEN
       iflag = 0
       CALL vneg(n,ga,d)
       RETURN
    END IF
      
   
! Initialization of indices.
 
    CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
      
         
! Computation of trans(sm)*ga and trans(um)*ga
    
    IF (iold <= 2) THEN
       CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),ga,ga,&
            tmpmc1(iold:),tmpmc2(iold:))
    ELSE
       CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
       CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
            ga,ga,tmpmc1(iold:),tmpmc2(iold:))
    END IF


! Computation of GAMMA

    IF (iscale >= 4) gamma = one


! Computation of two intermediate values tmpmc3 and tmpmc4

    IF (iold == 1) THEN
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mcnew,cm,tmpmc3,tmpmc4)
       CALL scsum(mcnew,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum(mcnew,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc4,tmpmc5,0,ierr)

    ELSE IF (iflag == 0) THEN
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mc+1,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mc+1,cm,tmpmc3,tmpmc4)
       CALL scsum(mc+1,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum(mc+1,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc4,tmpmc5,0,ierr)
    ELSE
       CALL trlieq(mcnew,mc,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mc,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mc,cm,tmpmc3,tmpmc4)
       CALL scsum(mc,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum(mc,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mc,iold,rm,tmpmc4,tmpmc5,0,ierr)
    END IF
      

! Computation of the search direction d
    
    IF (iold == 1) THEN
       CALL cwmaxv(n,mcnew,um,tmpmc3,d)
    ELSE 
       CALL cwmaxv(n,inew-1,um,tmpmc3,d)
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
       CALL xsumy(n,d,tmpn1,d)
    END IF

    CALL xdiffy(n,d,ga,d)
      
    IF (iold == 1) THEN
       CALL cwmaxv(n,mcnew,sm,tmpmc4,tmpn1)
       CALL scdiff(n,gamma,d,tmpn1,d)
    ELSE
       CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
       CALL xdiffy(n,d,tmpn1,d)
    END IF

  END SUBROUTINE dlskip
      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlsr1 *                                             *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*ga *
!*     by the limited memory SR1 update.                                *
!*                                                                      *
!************************************************************************

  SUBROUTINE dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm,&
       umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iiprint)
      
    USE param, ONLY : zero,small,one
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.   
         vneg, &   ! Copying of a vector with change of the sign.
         scalex, & ! Scaling a vector.
         xdiffy, & ! Difference of two vectors.
         scdiff, & ! Difference of the scaled vector and a vector.
         xsumy, &  ! Sum of two vectors.
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.
         calq, &   ! Solving x from linear equation A*x=y.
         copy, &   ! Copying of a vector.
         copy2     ! Copying of two vectors.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         ga, &     ! Current aggregate subgradient of the objective function.
         gp        ! Basic subgradient of the objective function.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         s, &      ! Difference of current and previous variables.
         u, &      ! Difference of current and previous subgradients.
         tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         tmpmc1, & ! Auxiliary array. On output: trans(sm)*ga.
         tmpmc2, & ! Auxiliary array. On output: trans(um)*ga.
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         nnk, &    ! Consecutive null steps counter.
         iiprint    ! Printout specification.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER(KIND=c_int), INTENT(OUT) :: & 
         isr1      ! Index of the type of L-SR1 update:
                   !   1  - SR1 update: the corrections are stored.
                   !   3  - SR1 update is skipped.
      
! Local Arrays
    REAL(KIND=c_double), DIMENSION(n) :: tmpn2
    REAL(KIND=c_double), DIMENSION((mcc+1)*(mcc+2)/2) :: tmpmat
    REAL(KIND=c_double), DIMENSION(mcc+1) :: tmpmc3,tmpmc4,tmpmc5,tmpmc6

! Local Scalars
    REAL(KIND=c_double) :: &
         stu, &    ! stu = trans(s)*u. 
         a, &      ! a = trans(ga) dm_(k-1) ga.
         b         ! b = trans(ga) dm_k ga.
    INTEGER(KIND=c_int) :: i,j,k, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         iflag2    ! Index for adaptive version.


    iflag2 = 0
    isr1 = 0 
      

! Initialization of indices
      
    CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,3)
      
         
! Computation of gamma 

    gamma = one


! Computation of trans(sm)*ga and trans(um)*ga

    IF (iold <= 2) THEN
       CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),ga,ga,&
            tmpmc1(iold:),tmpmc2(iold:))
    ELSE
       CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
       CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
            ga,ga,tmpmc1(iold:),tmpmc2(iold:))
    END IF


! Positive definiteness

    IF (-vdot(n,d,u) - vdot(n,tmpn1,s) >= -small) THEN
      
     
! SR1 update is skipped

       isr1 = 3
       
       IF (mcc == 0) THEN
          iflag = 0
          CALL vneg(n,ga,d)
          RETURN
       END IF

    ELSE
      
       stu = vdot(n,s,u)
        
       tmpmc1(inew) = vdot(n,s,ga)
       tmpmc2(inew) = vdot(n,u,ga)
      

!     Convergence conditions

       IF ((nnk == 1 .OR. mcc < mc) .OR. &
            (iflag == 1 .AND. (inew == 1 .OR. inew == mc))) THEN

! SR1 Update

          isr1 = 1
 

! Initialization of indices

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,1)        

          if (iflag2 == 1 .and. iold == 2) then
             tmpmc1(inew) = tmpmc1(1)
             tmpmc2(inew) = tmpmc2(1)
          end if


! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

         
! Update trans(sm)*gp and trans(um)*gp

          smtgp(inew) = vdot(n,s,gp)
          umtgp(inew) = vdot(n,u,gp)

     
! Computation of trans(sm)*u and trans(um)*u

          IF (iold <= 2) THEN
             CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                  u,u,tmpmc3(iold:),tmpmc4(iold:))
          ELSE
             CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                  u,u,tmpmc3(iold:),tmpmc4(iold:))
          END IF

          tmpmc3(inew) = stu
          tmpmc4(inew) = vdot(n,u,u)

         
! Update rm and umtum

          IF (mcc >= mc .AND. iflag2 /= 1) THEN
             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO
          END IF


          IF (inew >= mcnew) THEN
             CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                  tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy2(mcnew-inew,tmpmc3(iold:),&
                  rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                  umtum((mcnew-1)*mcnew/2+1:))
             CALL COPY2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                  umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF

      
! Update CM

          cm(inew) = stu

          inew = inew + 1
          IF (inew > mc + 1) inew = 1
          IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1

       ELSE


! Calculation of matrix (umtum-rm-trans(rm)+cm) from previous iteration
         
          DO  i=1,mcnew*(mcnew+1)/2
             tmpmat(i)= gamma * umtum(i) - rm(i)
          END DO

     
! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

          IF (iold == 1) THEN
             CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc5,0)

          ELSE IF (iflag == 0) THEN
             CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc5,0)

          ELSE
             CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc5,0)
          END IF


! Computation of a = -trans(ga)*dm_(k-1)*ga
      
          IF (iold <= 2) THEN
             CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
             CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
             CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             
          ELSE
             CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
             CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
             CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
                  tmpn1)
             CALL xdiffy(n,tmpn2,tmpn1,tmpn2)
             CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
                  tmpn1)
             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
          END IF
          
          a = vdot(n,ga,tmpn2)
          
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 2
             IF (iold > mc+1) iold = iold - mc - 1
          ELSE
             mcnew = mc - 1
             iold = inew + 2
             IF (iold > mc) iold = iold - mc
          END IF
      

! Calculation of the new canditate for search direction
! Updates are not necessarily saved

! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

     
! Computation of trans(sm)*u and trans(um)*u

          IF (iold == 1 .OR. iold == 2) THEN
             CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                  tmpmc3(iold:),tmpmc4(iold:))
          ELSE
             CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                  tmpmc3(iold:),tmpmc4(iold:))
          END IF
       
          tmpmc3(inew) = stu
          tmpmc4(inew) = vdot(n,u,u)


! Calculation of matrix (umtum-rm-trans(rm)+cm) without updating
! matrices rm, umtum and cm
      
          DO i=1,mcnew*(mcnew+1)/2
             tmpmat(i)= gamma * umtum(i) - rm(i)
          END DO

          DO i=1,mcnew-1
             j=(i-1)*i/2+1
             k=i*(i+1)/2+2
             CALL copy(i,tmpmat(k:),tmpmat(j:))
          END DO
         
          CALL scdiff(mcnew+1,gamma,tmpmc4,tmpmc3,tmpmc5)
         
          IF (inew >= mcnew) THEN
             CALL copy(mcnew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy(mcnew-inew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
             CALL copy(inew,tmpmc5,tmpmat((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF
      
          IF (iflag == 0) THEN
             CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc5,tmpmc5,iiprint)

          ELSE
             CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc,iold,tmpmat,tmpmc5,tmpmc5,iiprint)
          END IF


! Calculation of the new canditate for search direction d = -dm_k*ga
! and computation of b = -trans(ga)*dm_k*ga
      
          IF (iold <= 2) THEN
             CALL scalex(mcnew,gamma,tmpmc5(iold:),tmpmc6(iold:))
             CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc5(iold:),tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,d)
             CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc6(iold:),tmpn1)
             CALL xsumy(n,d,tmpn1,d)
         
          ELSE
             CALL scalex(mcnew+1,gamma,tmpmc5,tmpmc6)
             CALL cwmaxv(n,inew,sm,tmpmc5,tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,d)
             CALL cwmaxv(n,mcnew-inew,sm((iold-1)*n+1:),tmpmc5(iold:),&
                  tmpn1)
             CALL xdiffy(n,d,tmpn1,d)
             CALL cwmaxv(n,inew,um,tmpmc6,tmpn1)
             CALL xsumy(n,d,tmpn1,d)
             CALL cwmaxv(n,mcnew-inew,um((iold-1)*n+1:),tmpmc6(iold:),&
                  tmpn1)
             CALL xsumy(n,d,tmpn1,d)
          END IF

          b = vdot(n,ga,d)


! Checking the convergence conditions

          IF (b - a < zero) THEN
             isr1 = 3
             CALL copy(n,tmpn2,d)
            
          ELSE

             isr1 = 1
         
     
! Update trans(sm)*gp and trans(um)*gp
     
             smtgp(inew) = vdot(n,s,gp)
             umtgp(inew) = vdot(n,u,gp)

                     
! Update rm and umtum

             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO

             IF (inew >= mcnew) THEN
                CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                     umtum((mcnew-1)*mcnew/2+1:))
             ELSE
                CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                     umtum((mcnew-1)*mcnew/2+1:))
                CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                     umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
             END IF
            

! Update cm

             cm(inew) = stu
                     
             inew = inew + 1
             IF (inew > mc + 1) inew = 1
             IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
            
          END IF
         
          RETURN
         
       END IF
    END IF
      
    DO i=1,mcnew*(mcnew+1)/2
       tmpmat(i)= gamma * umtum(i) - rm(I)
    END DO
      
     
! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

    IF (iold == 1) THEN
       CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc4)
       CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
    ELSE IF (iflag == 0) THEN
       CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc4)
       CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
    ELSE
       CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc4)
       CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
    END IF
      

! Computation of the search direction d
      
    IF (iold <= 2) THEN
       CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
       CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
       CALL scdiff(n,-gamma,ga,tmpn1,d)
       CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
       CALL xsumy(n,d,tmpn1,d)
    ELSE
       CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
       CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
       CALL scdiff(n,-gamma,ga,tmpn1,d)
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
            tmpn1)
       CALL xdiffy(n,d,tmpn1,d)
       CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
       CALL xsumy(n,d,tmpn1,d)
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
            tmpn1)
       CALL xsumy(n,d,tmpn1,d)
    END IF
      
  END SUBROUTINE dlsr1

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE indic1 *                                            *
!*                                                                      *
!*     Initialization of indices.                                       *
!*                                                                      *
!************************************************************************

  SUBROUTINE indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,itype)

    IMPLICIT NONE

! Scalar Arguments

    INTEGER(KIND=c_int), INTENT(IN) :: & 
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         itype     ! Type of Initialization:
                   !   1  - corrections are stored,
                   !   2  - corrections are not stored,
                   !   3  - update is skipped.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         inew, &   ! Index for circular arrays.
         iflag, &  ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
         iflag2    ! Index for adaptive version.
                   !   0  - iflag has not been changed.
                   !   1  - iflag has been changed.
    INTEGER(KIND=c_int), INTENT(OUT) :: & 
         mcnew, &  ! Current size of vectors.
         iold      ! Index of the oldest corrections.
      
    IF (itype == 1) THEN
       IF (mcc < mc) THEN
          mcnew = mcc + 1
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 2
             IF (iold > mc+1) iold = iold - mc - 1
          ELSE
             IF (inew == 1) THEN
                inew = mc + 1
                mcnew = mc
                iold = 2
                iflag = 0
                iflag2 = 1
             ELSE IF (inew == mc) THEN
                mcnew = mc
                iold = 1
                iflag = 0
                iflag2 = 1
             ELSE
                mcnew = mc - 1
                iold = inew + 2
                IF (iold > mc) iold = iold - mc
             END IF
          END IF
       END IF
      
    ELSE IF (itype == 2) THEN
       IF (mcc < mc) THEN
          mcnew = mcc + 1
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc + 1
             iold = inew + 1
             IF (iold > mc + 1) iold = 1
          ELSE
             mcnew = mc
             iold = inew + 1
             IF (iold > mc) iold = 1
          END IF
       END IF

    ELSE 
       IF (mcc < mc) THEN
          mcnew = mcc
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 1
             IF (iold > mc + 1) iold = 1
          ELSE
             mcnew = mc - 1
             iold = inew + 1
             IF (iold > mc) iold = 1
          END IF
       END IF
    END IF
  END SUBROUTINE indic1
  
!************************************************************************
!*
!*     * SUBROUTINE agbfgs *
!*
!*     Computation of aggregate values by the limited memory BFGS update.
!*
!************************************************************************

  SUBROUTINE agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
       rm,cm,umtum,alfn,alfv,gamma,ic,rho)

    USE param, ONLY : zero,half,one
    USE lmbm_sub, ONLY : &
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         vdot      ! Dot product.

    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         u, &      ! Difference of trial and aggregate gradients.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         cm        ! Diagonal matrix.
    REAL(KIND=c_double), DIMENSION(:), INTENT(OUT) :: &
         ga        ! Next aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(OUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=c_double), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         ibfgs, &  ! Index of the type of BFGS update.
         ic        ! Correction indicator.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.

! Local arrays
    REAL(KIND=c_double), DIMENSION(mcc+1) :: tmpmc1, tmpmc2

! Local Scalars
    REAL(KIND=c_double) :: &
         p, &      ! p = trans(d)*u - alfn.
         q, &      ! q = trans(u)*dm*u, where dm is the inverse approximation of 
                   ! the Hessian calculated by using the L-BFGS formula.
         lam, &    ! Multiplier used to calculate aggregate values.
         w         ! Correction.
    INTEGER(KIND=c_int) :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.

! Intrinsic Functions
    INTRINSIC MAX,MIN,SIGN

    ierr = 0

    IF (mcc < mc) THEN
       IF (ibfgs == 2) THEN
          mcnew = mcc + 1
       ELSE
          mcnew = mcc
       END IF
       iold = 1

    ELSE
       IF (iflag == 0) THEN
          IF (ibfgs == 2) THEN
             mcnew = mc + 1
          ELSE
             mcnew = mc
          END IF
          iold = inew + 1
          IF (iold > mc+1) iold = 1

       ELSE
          IF (ibfgs == 2) THEN
             mcnew = mc
          ELSE
             mcnew = mc - 1
          END IF
          iold = inew + 1
          IF (iold > mc) iold = 1
       END IF
    END IF
      
      
! Computation of trans(d) * u - alfn

    p = vdot(n,d,u) - alfn
    q = vdot(n,u,u)

    IF (ic == 1) THEN
       w = rho * q
    ELSE
       w = zero
    END IF
         
     
! Computation of the product trans(u)*dm*u

    IF (mcc > 0 .OR. ibfgs == 2) THEN

       IF (iold == 1 .OR. ibfgs == 2) THEN
          CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc1,tmpmc2)
          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,tmpmc1,1,ierr)

          q = q - 2.0_c_double*vdot(mcnew,tmpmc2,tmpmc1)
          q = gamma*q
            
          DO i=1,mcnew
             tmpmc2(i) = cm(i)*tmpmc1(i)
          END DO

          q = q + vdot(mcnew,tmpmc1,tmpmc2)

          CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc2)

          q = q + gamma*vdot(mcnew,tmpmc1,tmpmc2)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc1,tmpmc2)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:), &
               u,u,tmpmc1(iold:),tmpmc2(iold:))
          CALL trlieq(mcnew,mcc,iold,rm,tmpmc1,tmpmc1,1,ierr)

          q = q - 2.0_c_double*(vdot(mcc-inew,tmpmc2(iold:),tmpmc1(iold:)) + &
               vdot(inew-1,tmpmc2,tmpmc1))
          q = gamma*q

          DO i=1,mcc
             tmpmc2(i) = cm(i)*tmpmc1(i)
          END DO

          q = q + vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
               vdot(inew-1,tmpmc1,tmpmc2)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc1,tmpmc2)

          q = q + gamma*(vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
               vdot(inew-1,tmpmc1,tmpmc2))
       END IF

    END IF
    
    q = q + w
    
    lam = half + SIGN(half,p)

    IF (q > zero) lam = MIN(one,MAX(zero,p/q))
      

! Computation of the aggregate values

    p = one - lam
    DO i=1,n
       ga(i)=lam*g(i) + p*gp(i)
    END DO
      
    alfv = lam*alfn
      
  END SUBROUTINE agbfgs

!************************************************************************
!*
!*     * SUBROUTINE aggsr1 *
!*
!*     Computation of aggregate values by the limited memory SR1 update.
!*
!************************************************************************
      
  SUBROUTINE aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv, &
       umtum,rm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,one,small
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.
         scalex, & ! Scaling a vector.
         xsumy, &  ! Sum of two vectors.
         xdiffy, & ! Difference of two vectors.
         scsum, &  ! Sum of a vector and the scaled vector.
         scdiff, & ! Difference of the scaled vector and a vector.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.   
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         lineq, &  ! Solver for linear equation.
         calq      ! Solving x from linear equation A*x=y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         smtga, &  ! vector smtga = trans(sm)*ga.
         umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         ga        ! Aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=c_double), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         icn       ! Correction indicator.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
! Local arrays
    REAL(KIND=c_double), DIMENSION(n) :: tmpn2,tmpn3,tmpn4
    REAL(KIND=c_double), DIMENSION((mcc+1)*(mcc)/2) :: tmpmat
    REAL(KIND=c_double), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

! Local Scalars
    REAL(KIND=c_double) :: &
         pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                   ! presents the L-SR1- approximation of Hessian.
         rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
         prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
         rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
         qr, &     ! qr = trans(g-ga) dm (g-ga).
         pq, &     ! pq = trans(g-gp) dm (g-gp).
         qqp, &    ! qqp = trans(g-gp) dm g + alfn.
         lam1, &   ! Multiplier used to calculate aggregate values.
         lam2, &   ! Multiplier used to calculate aggregate values.
         w, &      ! Correction.
         tmp1, &   ! Auxiliary scalar.
         tmp2      ! Auxiliary scalar.
    INTEGER(KIND=c_int) :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.

      
! Intrinsic Functions
    INTRINSIC MIN,MAX


    ierr = 0
      
    IF (mcc < mc) THEN
       iold = 1
       mcnew = mcc
    ELSE IF (iflag == 0) THEN
       mcnew = mc
       iold = inew + 1
       IF (iold > mc+1) iold = 1
    ELSE
       mcnew = mc - 1
       iold = inew + 1
       IF (iold > mc) iold = 1
    END IF
      
    CALL xdiffy(n,gp,ga,tmpn2)
     
      
! Calculation of tmpn3 = trans(gp - ga)dm

    IF (mcc > 0) THEN

       DO i=1,mcnew*(mcnew+1)/2
          tmpmat(i)= gamma * umtum(i) - rm(i)
       END DO

       IF (iold == 1) THEN
          CALL xdiffy(mcnew,umtgp,umtga,tmpmc4)
          CALL scdiff(mcnew,gamma,tmpmc4,smtgp,tmpmc4)
          CALL xsumy(mcnew,tmpmc4,smtga,tmpmc4)
          
          CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,0)
          CALL scalex(mcnew,gamma,tmpmc3,tmpmc4)
          
          CALL cwmaxv(n,mcnew,sm,tmpmc3,tmpn4)
          CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
          CALL cwmaxv(n,mcnew,um,tmpmc4,tmpn4)
          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)

       ELSE
          CALL xdiffy(mcc,umtgp,umtga,tmpmc4)
          CALL scdiff(mcc,gamma,tmpmc4,smtgp,tmpmc4)
          CALL xsumy(mcc,tmpmc4,smtga,tmpmc4)

          CALL calq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,0)
          CALL scalex(mcc,gamma,tmpmc3,tmpmc4)
          
          CALL cwmaxv(n,inew-1,sm,tmpmc3,tmpn4)
          CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
          CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc3(iold:)&
               ,tmpn4)
          CALL xsumy(n,tmpn3,tmpn4,tmpn3)
          CALL cwmaxv(n,inew-1,um,tmpmc4,tmpn4)
          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
          CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc4(iold:)&
               ,tmpn4)
          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
       END IF
       
       IF (icn == 1) THEN
          CALL scsum(n,rho,tmpn2,tmpn3,tmpn3)
       END IF
         
       pr = vdot(n,tmpn3,tmpn2)
       rrp = vdot(n,tmpn3,ga) 
       CALL xdiffy(n,g,ga,tmpn4)
       prqr = vdot(n,tmpn3,tmpn4)
       rrq = -vdot(n,tmpn4,d)

    ELSE

       pr = vdot(n,tmpn2,tmpn2)
       rrp = vdot(n,tmpn2,ga) 
       CALL xdiffy(n,g,ga,tmpn4)
       prqr = vdot(n,tmpn2,tmpn4)
       rrq = -vdot(n,tmpn4,d)
    END IF

! calculation of qr = trans(g - ga) dm (g - ga)

    qr = vdot(n,tmpn4,tmpn4)
    IF (icn == 1) THEN
       w = rho*qr
    ELSE
       w = zero
    END IF
      
    IF (mcc > 0) THEN
       qr = gamma*qr

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
          CALL scsum(mcnew,-gamma,tmpmc3,tmpmc4,tmpmc4)
          CALL lineq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,ierr)
            
          qr = qr - vdot(mcnew,tmpmc4,tmpmc3) + w

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
          CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
               um((iold-1)*n+1:),tmpn4,tmpn4,tmpmc4(iold:),tmpmc3(iold:))
          CALL scsum(mcc,-gamma,tmpmc3,tmpmc4,tmpmc4)
          CALL lineq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,ierr)
          
          qr = qr - vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) -&
               vdot(inew-1,tmpmc4,tmpmc3) + w
       END IF

    END IF
      
    pq = qr - prqr - prqr + pr
    qqp = pq + prqr + rrq - pr - rrp + alfn
    rrp = rrp - alfv
    rrq = rrq + alfn - alfv

! computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
       tmp1 = rrq/qr
       tmp2 = prqr/qr
       w = pr - prqr*tmp2
       IF (w /= zero) THEN
          lam1 = (tmp1*prqr - rrp)/w
          lam2 = -tmp1 - lam1*tmp2
          IF (lam1*(lam1 - one) < zero .AND. &
               lam2*(lam1 + lam2 - one) < zero) GO TO 200
       END IF
    END IF

! Minimum on the boundary

100 continue
    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq+rrq)*lam2
!    w = (lam2*qr + 2.0_c_double*rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
!    tmp2 = (tmp1*pr + 2.0_c_double*rrp)*tmp1
    tmp2 = (tmp1*pr + rrp+rrp)*tmp1
    IF (tmp2 < w) THEN
       w = tmp2
       lam1 = tmp1
       lam2 = zero
    END IF
    
    IF (qqp*(qqp - pq) < zero) THEN
       IF (qr + rrq + rrq - qqp*qqp/pq < W) THEN
          lam1 = qqp/pq
          lam2 = one - lam1
       END IF
    END IF
    
200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
         .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
         lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)

! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
       ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
    
  END SUBROUTINE aggsr1
      
!************************************************************************
!*
!*     * SUBROUTINE agskip *
!*
!*     Computation of aggregate values after consecutive null steps
!*     by the limited memory BFGS update.
!*
!************************************************************************
      
  SUBROUTINE agskip(n,mc,mcc,inew,iflag,g,gp,ga,d,u,alfn,alfv, &
       umtum,rm,cm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,half,one,small
    USE lmbm_sub, ONLY : &
         xdiffy, & ! Difference of two vectors.
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation l*x=y or trans(l)*x=y.
         vdot      ! Dot product
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         u, &      ! Difference of trial and aggregate gradients.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         smtga, &  ! vector smtga = trans(sm)*ga.
         umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=c_double), DIMENSION(:), INTENT(INOUT) :: &
         ga        ! Aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=c_double), INTENT(INOUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=c_double), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         icn       ! Correction indicator.
    INTEGER(KIND=c_int), INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
! Local arrays
    REAL(KIND=c_double), DIMENSION(n) :: tmpn2
    REAL(KIND=c_double), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

! Local Scalars
    REAL(KIND=c_double) :: &
         pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                   ! presents the L-SR1- approximation of Hessian.
         rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
         prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
         rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
         qr, &     ! qr = trans(g-ga) dm (g-ga).
         pq, &     ! pq = trans(g-gp) dm (g-gp).
         qqp, &    ! qqp = trans(g-gp) dm g + alfn.
         lam1, &   ! Multiplier used to calculate aggregate values.
         lam2, &   ! Multiplier used to calculate aggregate values.
         w, &      ! Correction.
         tmp1, &   ! Auxiliary scalar.
         tmp2      ! Auxiliary scalar.
    INTEGER(KIND=c_int) :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.
    
! Intrinsic Functions
    INTRINSIC MIN,MAX

    ierr = 0
           
    IF (mcc < mc) THEN
       iold = 1
       mcnew = mcc
    ELSE
       IF (iflag == 0) THEN
          mcnew = mc
          iold = inew + 1
          IF (iold > mc+1) iold = 1
       ELSE
          mcnew = mc - 1
          iold = inew + 1
          IF (iold > mc) iold = 1
       END IF
    END IF

      
! Calculation of pq = trans(g-gp) dm (g-gp) = trans(u) dm u.

    pq = vdot(n,u,u)

    IF (icn == 1) THEN
       w = rho * pq
    ELSE
       w = zero
    END IF
    
    IF (mcc > 0) THEN

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc3,tmpmc4)
          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
          
          pq = pq - 2.0_c_double*vdot(mcnew,tmpmc4,tmpmc3)
          pq = gamma*pq
          
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pq = pq + vdot(mcnew,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
          pq = pq + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
               u,u,tmpmc3(iold:),tmpmc4(iold:))
          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

          pq = pq - 2.0_c_double*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          pq = gamma*pq
          
          DO i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pq = pq + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
          pq = pq + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               + vdot(inew-1,tmpmc3,tmpmc4))
       END IF

    END IF

    pq = pq + w
      

! Calculation of pr = trans(gp-ga) dm (gp-ga).
      
    CALL xdiffy(n,gp,ga,tmpn2)
    pr = vdot(n,tmpn2,tmpn2)

    IF (icn == 1) THEN
       w = rho * pr
    ELSE
       w = zero
    END IF

    IF (mcc > 0) THEN
       
       IF (iold == 1) THEN 
          DO i=1, mcnew
             tmpmc3(i)=smtgp(i)-smtga(i)
             tmpmc4(i)=umtgp(i)-umtga(i)
          END DO
          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
             
          pr = pr - 2.0_c_double*vdot(mcnew,tmpmc4,tmpmc3)
          pr = gamma*pr
            
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pr = pr + vdot(mcnew,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)

          pr = pr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          DO i=1, mcc
             tmpmc3(i)=smtgp(i)-smtga(i)
             tmpmc4(i)=umtgp(i)-umtga(i)
          END DO
          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

          pr = pr - 2.0_c_double*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          pr = gamma*pr

          DO  i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pr = pr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)

          pr = pr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               + vdot(inew-1,tmpmc3,tmpmc4))
       END IF

    END IF

    pr = pr + w

      
! Calculation of rrp = trans(gp-ga) dm ga - alfv.
      
    rrp = - vdot(n,tmpn2,d) - alfv
      

! Calculation of qr = trans(g-ga) dm (g-ga).

    CALL xdiffy(n,g,ga,tmpn2)
    qr = vdot(n,tmpn2,tmpn2)

    IF (icn == 1) THEN
       w = rho * qr
    ELSE
       w = zero
    END IF

    IF (mcc > 0) THEN

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)

          qr = qr - 2.0_c_double*vdot(mcnew,tmpmc4,tmpmc3)
          qr = gamma*qr
            
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          qr = qr + vdot(mcnew,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
          qr = qr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
               tmpn2,tmpn2,tmpmc3(iold:),tmpmc4(iold:))
          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

          qr = qr - 2.0_c_double*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          qr = gamma*qr
          
          DO i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO
          
          qr = qr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
          qr = qr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               +vdot(inew-1,tmpmc3,tmpmc4))
       END IF
       
    END IF
    
    qr = qr + w
      

! Calculation of rrq = trans(g-ga) dm ga - alfv + alfn.

    rrq = - vdot(n,tmpn2,d) - alfv + alfn

     
! Calculation of prqr = trans(gp-ga) dm (g-ga).
      
    prqr = half*(qr - pq + pr)

     
! Calculation of qqp = trans(g-gp) dm g + alfn.

    qqp = pq + prqr + rrq - pr - rrp

     
! Computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
       tmp1 = rrq/qr
       tmp2 = prqr/qr
       w = pr - prqr*tmp2
       IF (w /= zero) THEN

          lam1 = (tmp1*prqr - rrp)/w
          lam2 = -tmp1 - lam1*tmp2

          IF (lam1*(lam1 - one) < zero .AND. &
               lam2*(lam1 + lam2 - one) < zero) GO TO 200
       END IF
    END IF


! Minimum on the boundary

    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq + rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
    tmp2 = (tmp1*pr + rrp + rrp)*tmp1
    IF (tmp2 < w) THEN
       w = tmp2
       lam1 = tmp1
       lam2 = zero
    END IF
      
    IF (qqp*(qqp - pq) < zero) THEN
       IF (qr + rrq + rrq - qqp*qqp/pq < w) THEN
          lam1 = qqp/pq
          lam2 = one - lam1
       END IF
    END IF

200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
         .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
         lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)
      

! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
       ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
      
  END SUBROUTINE agskip

!************************************************************************
!*
!*     * SUBROUTINE wprint *
!*
!*     Printout the warning and error messages.
!*
!************************************************************************
      
  SUBROUTINE wprint(iterm,iiprint,nout)
    IMPLICIT NONE

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: &
         iiprint, &    ! Printout specification:
                      !  -1  - No printout.
                      !   0  - Only the error messages.
                      !   1  - The final values of the objective
                      !        function.
                      !   2  - The final values of the objective
                      !        function and the most serious
                      !        warning messages.
                      !   3  - The whole final solution. 
                      !   4  - At each iteration values of the
                      !        objective function.
                      !   5  - At each iteration the whole
                      !        solution
         nout, &      ! Auxilary printout specification.
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    IF (iiprint >= 0) THEN

! Initial error messages

       IF (iterm == -5) THEN
          IF (nout == 1) WRITE (6,FMT='(1X,''Error: '' &
               ''Number of variables (n) is too small, iterm='',I3)') iterm
          IF (nout == 2) WRITE (6,FMT='(1X,''Error: '' &
               ''The maximum number of stored corrections (mcu) '' &
               ''is too small, iterm='',I3)') iterm
          IF (nout == 3) WRITE (6,FMT='(1X,''Error: '' &
               ''The size of the bundle (na) is too small, iterm='' &
               ,I3)') iterm
          IF (nout == 4) WRITE (6,FMT='(1X,''Error: '' &
               ''Line search parameter epsl >= 0.25, iterm='',I3)') iterm
          RETURN
       END IF

        
! Warning messages

       IF (iiprint >= 2) THEN
          IF (iterm == 0) THEN
             IF (nout == -1) WRITE (6,FMT='(1X,''Warning: '' &
                  ''mc > mcu. Assigned mc = mcu.'')')
             IF (nout == -2) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A line search parameter epsr >= 0.5.'')')
             IF (nout == -3) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A nondescent search direction occured. Restart.'')')
             IF (nout == -4) WRITE (6,FMT='(1X,''Warning: '' &
                  ''Does not converge.'')')
             IF (nout == -5) WRITE (6,FMT='(1X,''Warning: '' &
                  ''tmax < tmin. Restart.'')')
             RETURN
          END IF
         

! Printout the final results
            
          IF (iterm == 6) WRITE (6,FMT='(1X,''Abnormal exit: Time is up.'')')
          IF (iterm == 7) WRITE (6,FMT='(1X,''Abnormal exit: f < tolb.'')')
          IF (iterm == 2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Too many steps without significant progress.'')')
          IF (iterm == 3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''The value of the function does not change.'')')
          IF (iterm == 5) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of iterations > '',I5)') nout
          IF (iterm == 4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of function evaluations > '',I5)') nout
          IF (iterm == -1) THEN
             IF (nout == -1) THEN
                WRITE (6,FMT='(1X,''Abnormal exit: Two consecutive restarts.'')')
             ELSE
                WRITE (6,FMT='(1X,''Abnormal exit: tmax < tmin in two subsequent iterations.'')')
             END IF
          END IF
          IF (iterm == -2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of restarts > '',I5''.'')') nout
          IF (iterm == -3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in function or subgradient calculations.'')')
          IF (iterm == -4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in attaining the demanded accuracy.'')')
       END IF
    END IF
      
  END SUBROUTINE wprint

            
!************************************************************************
!*
!*     * SUBROUTINE rprint *
!*      
!*     Printout the (final) results.
!*
!************************************************************************
      
  SUBROUTINE rprint(n,nit,nfe,nge,x,f,wk,qk,iterm,iiprint)
    IMPLICIT NONE

! Scalar Arguments
    INTEGER(KIND=c_int), INTENT(IN) :: & 
         n, &         ! Number of variables 
         nit, &       ! Number of used iterations.
         nfe, &       ! Number of used function evaluations.
         nge, &       ! Number of used subgradient evaluations.
         iiprint, &    ! Printout specification:
                      !  -1  - No printout.
                      !   0  - Only the error messages.
                      !   1  - The final values of the objective
                      !        function.
                      !   2  - The final values of the objective
                      !        function and the most serious
                      !        warning messages.
                      !   3  - The whole final solution. 
                      !   4  - At each iteration values of the
                      !        objective function.
                      !   5  - At each iteration the whole
                      !        solution
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    REAL(KIND=c_double), INTENT(IN) :: &
         f, &         ! Value of the objective function.
         wk, &        ! Value of the first stopping criterion.
         qk           ! Value of the second stopping criterion.

! Array Arguments
    REAL(KIND=c_double), DIMENSION(:), INTENT(IN) :: &
         x            ! Vector of variables
         
! Local Scalars
    INTEGER(KIND=c_int) :: i

! Intermediate results
    
    IF (iterm == 0) THEN
       IF (iiprint > 3) WRITE (6,FMT='(1X,''nit='',I5,2X, &
            ''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
            ''qk='',D11.4,2X)') nit,nfe,nge,f,wk,qk
       IF (iiprint == 5) WRITE (6,FMT='(1X,''x='', &
            5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
       RETURN
    END IF
         

! Final results

    IF (iiprint > 0) WRITE (6,FMT='(1X,''nit='',I5,2X, &
         ''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
         ''qk='',D11.4,2X,''iterm='',I3)') nit,nfe,nge,f,wk,qk,iterm
    IF (iiprint .EQ. 3 .OR. iiprint .EQ. 5) &
         WRITE (6,FMT='(1X,''x='',5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
      
  END SUBROUTINE rprint
      
END MODULE lmbm_mod
!LMBM ENDS


  MODULE oscar
    
     USE, INTRINSIC :: iso_c_binding
     USE omp_lib

     USE functions                 ! INFORMATION from the USER
     USE bundle1                   ! Bundle 1
     USE bundle2                   ! Bundle 2
     USE dbdc                      ! DBDC method
	 USE initializat		   	   ! Initialization of parameters for LMBM
     USE lmbm_mod		   		   ! LMBM method implementation
         
     IMPLICIT NONE   
        
     PRIVATE    
     PUBLIC :: oscar_cox 
     PUBLIC :: oscar_mse 
     PUBLIC :: oscar_logistic 
        
        CONTAINS 
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      oscar_cox        : The double bundle method for Cox's model                 | |
        !| |      oscar_mse        : The double bundle method for mean squared error          | |
        !| |      oscar_logistic   : The double bundle method for logistic regression         | |
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
        ! Subroutine for oscar with Cox's proportional hazard model
        
    !--------------------------------------------------------------------------
    ! ^^^^ START: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------        
          SUBROUTINE oscar_cox(x, y, kits, costs, nrow, ncol, nkits, beta, fperk, &
       & in_print, in_start, in_k_max, &
       & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
       & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, &
       & nKitOnes, betakits, &
       & solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
       & in_tolg, in_tolg2, in_eta, in_epsL ) & 
        & BIND(C, name = "oscar_cox_f_")               
    !--------------------------------------------------------------------------
    ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------


          !--------------------------------------------------------------------------            
          ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------    
          ! SUBROUTINE oscar_cox(infileX, infileY, infileK, infileC,  &
          !                     & nrow, ncol, nkits, &
          !                     & beta, fperk, in_print, in_start, in_k_max, &
          !                     & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
          !                     & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, & 
          !                     & nKitOnes, betakits, solver_id, & 
          !                     & in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
          !                     & in_tolg, in_tolg2, in_eta, in_epsL)    
          !--------------------------------------------------------------------------            
          ! ^^^^ END: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------   

            !_____________________________________________________________________________________
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            !
            !  !! NOTICE !!  The DATA is SCALED during this SUBROUTINE !
            !
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            ! 
            !           
            ! * 'problem = 3'        : the used objective function is Cox's proportional hazards model with L0-norm for kits
            !                          (parametric approach since the number of nonzero elements is fixed) 
            !
            !   Solves in a loop all L0-norm problems with the fixed number of k=1,...,in_k_max of kits    (if in_start > 0)
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'ncol'        : The dimension of the problem = the number of features in a predictor, INTEGER
            !         * 'nrow'        : The number of records (data points), INTEGER
            !         * 'nkits'       : the (maximum) number of kits, INTEGER
            !
            !         * 'in_print'    : specifies the print, INTEGER
            !         * 'in_start'    : specifies how starting points are selected when the L0-norm problem is solved for  
            !                           the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'in_k_max'    : specifies how many kits are used in the last problem of the loop  (NEEDS to be 1 <= 'in_k_max' <= nkits !)
            !
            !         * 'nKitOnes'    : the number of value 1 in kit matrix
			!
            !         * 'solver_id'   : the used solver (1=DBDC and 2=LMBM) INTEGER
            !         
            !         PARAMETERS in DBDC method:
            !
            !         * in_mrounds     : the maximum number of rounds in one main iteratioin
            !         * in_mit         : the maximum number of main iterations
            !         * in_mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance          
            !
            !         PARAMETERS in LMBM method:
            !
            !         * in_na        : Size of the bundle na >= 2.
            !         * in_mcu       : Upper limit for maximum number of stored corrections, mcu >= 3.
            !         * in_mcinit    : Initial maximum number of stored corrections, mcu >= mcinit >= 3.
            !                          If mcinit <= 0, the default value mcinit = 3 will be used. 
            !                          However, the value mcinit = 7 is recommented.             
            !         * in_tolf      : Tolerance for change of function values (default = 1.0E-5).
            !         * in_tolf2     : Second tolerance for change of function values.
            !                          - If tolf2 < 0 the the parameter and the corresponding termination 
            !                            criterion will be ignored. 
            !                          - If tolf2 = 0 the default value 1.0E+4 will be used
            !         * in_tolg      : Tolerance for the first termination criterion (default = 1.0E-5).   
            !         * in_tolg2     : Tolerance for the second termination criterion (default = tolg). clustering code small data  
            !         * in_eta       : Distance measure parameter, eta >= 0.   
            !         * in_epsL      : Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.)  
            !
            !         NOTICE: DATA IS GIVEN IN VECTOR FORMAT ! Due to this values for observations/kits are given in these vectors one after another
            !
            !         * 'x'          : the matrix X of input variables in vector form --> will be stored to in_mX (each row is one observation)
            !                               REAL, DIMENSION(nrow*ncol)  
            !         * 'y'          : the matrix Y of y_time and y_event in vector form --> will be stored to in_mY (each row cosists of time and event)
            !                               INTEGER, DIMENSION(nrow*2) 
            !         * 'costs'      : the cost vector C for kits --> will be stored to in_mC
            !                               REAL, DIMENSION(knits)
            !         * 'kits'       : the matrix K of kit structures in vector form (values are either 0 or 1) --> will be stored to in_mK (each row is one kit)
            !                               INTEGER, DIMENSION(knits*ncol)
            !
            !
            ! OUTPUT: * 'fperk'      : The vector containing objective function values at for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(nkits) 
            !         * 'beta'       : The vector containing solution vectors beta obtained for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(ncol*nkits) 
            !
            ! NOTICE: * The dimension of vectors 'fperk' has to be 'nkits' ('nkits' is the maximum number of kits)
            !         * The dimension of vectors 'beta' has to be 'nft*nkits' ('nft' is the dimension of the problem and 'nkits' is the number of kits)
            !
            ! OPTIONAL INPUT (CAN BE INCLUDED AS INPUT PARAMETERS IF NEEDED. AT THE MOMENT DEFAULT VALUES ARE USED FOR THEM.):
            !         * 'mit'            : The maximum number of 'main iterations', INTEGER                  
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration', INTEGER
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Escape procedure', INTEGER
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm, LOGICAL           
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm, LOGICAL
            !         * 'scale_in_use'   : If .TRUE. then the data is scaled
            !         * 'CPUtime'       : the CPU time (in seconds) REAL
            !
            ! NOTICE: * 'in_print' has to be 0, 1, 2 or 3. If it is NOT then DEFAULT value 1 is used.             
            !         * 'in_start' has to be  1, 2, 3 or 4. If it is NOT then DEFAULT value 2 is used.             
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! INPUTs
               ! INTEGER(KIND = c_int), INTENT(IN) :: nrow                  ! Number of rows in x (i.e. records)
               ! INTEGER(KIND = c_int), INTENT(IN) :: ncol                  ! Number of cols in x (i.e. features)
               ! INTEGER(KIND = c_int), INTENT(IN) :: nkits                 ! Number of kits for features
               
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_start              ! Starting point procedure used
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_print              ! Print used
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_k_max              ! The maximum number of kits in the loop
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds              ! the maximum number of rounds in one main iteration
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mit                  ! the maximum number of main iterations
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               ! REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               ! REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               ! REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: nKitOnes                ! the number of value 1 in kit matrix
			   ! INTEGER(KIND=c_int), INTENT(IN) :: solver_id               ! the used solver (1=DBDC and 2=LMBM)

!**
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_na                  ! Size of the bundle na >= 2.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcu                 ! Upper limit for maximum number of stored corrections, mcu >= 3.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcinit              ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! ! However, the value mcinit = 7 is recommented.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf                ! Tolerance for change of function values (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf2               ! Second tolerance for change of function values.
                                                                          ! !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          ! !   criterion will be ignored. 
                                                                          ! !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg                ! Tolerance for the first termination criterion (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg2               ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               ! REAL(KIND=c_double), INTENT(IN) ::  in_eta                 ! Distance measure parameter, eta >= 0.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_epsL                ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**

            ! INPUTs
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: nrow                  ! Number of rows in x (i.e. records)
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: ncol                  ! Number of cols in x (i.e. features)
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: nkits                 ! Number of kits for features
               
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: in_start              ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: in_print              ! Print used
               INTEGER(KIND = c_int), INTENT(IN), VAlUE :: in_k_max              ! The maximum number of kits in the loop
               
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_mrounds              ! the maximum number of rounds in one main iteration
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_mit                  ! the maximum number of main iterations
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN), VAlUE :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN), VAlUE :: in_crit_tol             ! the stopping tolerance: criticality tolerance
             
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: nKitOnes                ! the number of value 1 in kit matrix
             
			   INTEGER(KIND=c_int), INTENT(IN), VALUE :: solver_id               ! the used solver (1=DBDC and 2=LMBM)

  !**
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_na           ! Size of the bundle na >= 2.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcu          ! Upper limit for maximum number of stored corrections, mcu >= 3.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcinit       ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! However, the value mcinit = 7 is recommented.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf         ! Tolerance for change of function values (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf2        ! Second tolerance for change of function values.
                                                                          !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          !   criterion will be ignored. 
                                                                          !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg         ! Tolerance for the first termination criterion (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg2        ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_eta          ! Distance measure parameter, eta >= 0.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_epsL         ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**
            !--------------------------------------------------------------------------
            ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nrow*ncol)  :: x        !Vector of data values
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nrow*2)     :: y        !Vector of response values (2-column survival, time + event)
               INTEGER(KIND = c_int), INTENT(IN), DIMENSION(nkits*ncol) :: kits     !Vector of kit indicator values (binary indicators)
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nkits)      :: costs    !Costs associated to each kit
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------


            !--------------------------------------------------------------------------            
            ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
            !--------------------------------------------------------------------------
            !   CHARACTER(LEN=80), INTENT(IN) :: infileX         ! The name of "predictor matrix" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileY         ! The name of "observed time and label" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileC         ! The name of "kit costs" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileK         ! The name of "predictor matrix" file   
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^          
            !--------------------------------------------------------------------------
            
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(ncol*nkits)   :: beta       !Output variable for beta coefficients per k
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nkits)        :: fperk      !Output variable target function value per k     
               INTEGER(KIND = c_int), INTENT(OUT), DIMENSION(nkits*nkits)  :: betakits   !!Output variable telling kits in beta coefficients per k     
           
           !***************************** LOCAL VARIABLES ************************************      
           
               TYPE(set_info) :: info                         ! The set of information                     
 
               REAL(KIND=c_double) :: CPUtime                 ! the CPU time (in seconds)
			   
               INTEGER(KIND=c_int) :: solver                  ! defines the solver used (1=DBDC, 2=LMBM)			   

               INTEGER(KIND=c_int) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int) :: nk                      ! the number of kits 
               INTEGER(KIND=c_int) :: nk_max                  ! the maximum number of kits in the loop 
   
               REAL(KIND=c_double), DIMENSION(nKitOnes) :: beta_solution    ! the solution vector beta obtained for the problem

               REAL(KIND=c_double), DIMENSION(nKitOnes,nkits) :: points     ! the beta_solutions for problem 3 for fixed k ('nkits' different starting points) 
               REAL(KIND=c_double), DIMENSION(nkits) ::      f_points       ! the objective function values for problem 3 for fixed k ('nkits' different starting points)
               REAL(KIND=c_double), DIMENSION(nKitOnes) ::       x_koe      ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nKitOnes) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller
               
               REAL(KIND=c_double), DIMENSION(nKitOnes,nkits) :: x_solution    ! the beta solutions for the problem with fixed starting points and number of kits 
               REAL(KIND=c_double), DIMENSION(nkits) :: f_solution             ! the objective function valuse at the solution 'x_solution'
               REAL(KIND=c_double) :: f_solution_DBDC                          ! the f_solution obtained from DBDC method
               
               REAL(KIND=c_double), DIMENSION(nrow,ncol) :: in_mX      ! predictor matrix (row is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrow,2) :: in_mY             ! observed times and labels matrix (row is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nkits,ncol) :: in_mK     ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nkits) :: in_mC              ! kit costs          
       
               REAL(KIND=c_double), DIMENSION(nKitOnes,nrow) :: mXt        ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(2,nrow) :: mYt               ! observed times and labels matrix (column is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nkits,nKitOnes) :: mK        ! Modified kit matrix (There is only one value one in each column)  
               REAL(KIND = c_double), DIMENSION(nKitOnes*nkits)  :: beta_nft   !Output variable for beta coefficients per k

               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta          ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta_ed       ! indices of kits in the previous solution 'x_ed'

               REAL(KIND=c_double), DIMENSION(nKitOnes) :: x_0             ! the starting point
               
               REAL(KIND=c_double), DIMENSION(nrow) :: mTimes          ! Times for the observations   
               INTEGER(KIND=c_int), DIMENSION(nrow) :: mTimesInd       ! Labels of times for the observations (0=alive, 1=death)  
               
               INTEGER(KIND=c_int), DIMENSION(ncol) :: KitOnes        ! The number of kits where each variables is located. If kits do not intersect then each variable is only in one kit!
              
               REAL(KIND=c_double), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=c_int) :: nstart              ! the number of starting point
               INTEGER(KIND=c_int) :: start_max           ! the number of starting point when 'start = 5'

               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                       
               INTEGER(KIND=c_int) :: mit                ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds            ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds_clarke     ! the maximum number of rounds during one 'Clarke stationary' algorithm
                                                         ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used

               ! Possible USER PARAMETER
               INTEGER(KIND=c_int) :: iprint   ! specifies the print
                                               !   iprint = 0 : print is suppressed
                                               !   iprint = -1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Not table form)
                                               !   iprint = 1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Table form)
                                               !   iprint = 2 : prints for each L0-norm problem the final value of Cox's proportional hazard model obtained from each starting point
                                               !   iprint = 3 : prints for each starting point and L0-norm problem all intermediate results together the final results
 
               INTEGER(KIND=c_int) :: start    !   start = 1  : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                               !                (starting point is the solution obtained for L0-norm problem with k-1 nonzero components)
                                               !
                                               !   start = -1 : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1
                                               !                (starting point is the solution obtained for L0-norm problem with k+1 nonzero components)  
                                               !
                                               !   start = 2  : for L0-norm problem with k nonzero components uses 'nkits-k+1' starting points
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
                               
 
               LOGICAL :: agg_used              ! .TRUE. if aggregation is used in the algorithm. Otherwise .FALSE.                                                           
               LOGICAL :: stepsize_used         ! .TRUE. if simple stepsize determination is used in the algorithm. Otherwise .FALSE.
               
               LOGICAL :: scale_in_use          ! If .TRUE. data is scaled
               
               LOGICAL :: kit_in_use            ! If .TRUE. kit is used in the solution
               LOGICAL :: run_stop              ! If .TRUE. run is stopped for selected k
               LOGICAL :: mukana                ! If .TRUE. specific kit is in the solution
               
               LOGICAL :: new_start             ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
 
               LOGICAL :: ed_sol_in_pen         ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
 
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
               REAL(KIND=c_double) :: small           ! The cost of solution beta
               
               REAL(KIND=c_double) :: s_time, f_time  ! The start and finish times
               REAL(KIND=c_double) :: cpu             ! The cpu time
               
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2

               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
                              
               REAL(KIND=c_double) :: random_num                     ! Random number
               REAL(KIND=c_double), DIMENSION(nkits) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=c_int) :: problem1              ! The DC component f1 
               INTEGER(KIND=c_int) :: problem2              ! The DC component f2 

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
                       
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: nremoved
               INTEGER(KIND=c_int) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=c_int) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii
               INTEGER(KIND=c_int) :: l, l2
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization

               
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mPrNum   ! The number of problems for threads
               INTEGER(KIND=c_int) :: tread_num                  ! The number of thread
               INTEGER(KIND=c_int) :: jaannos                    ! The residue in division
               INTEGER(KIND=c_int) :: kokoosa                    ! The integer obtained in division
               INTEGER(KIND=c_int) :: startind, j3               ! Help variables

               REAL(KIND=c_double) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER(KIND=c_int) :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time 
			   
!**			   
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.               
 !**
 
               CALL cpu_time(s_time)   ! Start CPU timing
               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time     

!**
              ! The solver is defined/selected
               IF (solver_id == 1) THEN
                   solver = 1 !DBDC
               ELSE IF (solver_id == 2) THEN
                   solver = 2 !LMBM
               ELSE 
                   solver = 1 !DBDC              
               END IF
!**
              ! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
               
               ! Maximum number of possible treads in parallellization
               max_threads = omp_get_max_threads()
               
!**
			   IF (solver == 2) THEN  !No parallellization with LMBM can be used
                   max_threads = 1
			   END If	   
!**
               
               CALL omp_set_num_threads(max_threads)
               tread_num = max_threads

               WRITE(*,*) 'The number of threads', max_threads             
               ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----   

               ! The initialization of parametrs used in DBDC methods
               CALL allocate_parameters(info, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)

!**                                          
                 ! The initialization of parameters used in LMBM method
                 CALL init_par(in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, &
                                & in_tolg2, in_eta, in_epsL)
!**
               
               ! Set the number of rows and columns inside Fortran  + kits           
               nrecord = nrow
               nft = nKitOnes
               nk = nkits
            
               ! The maximum number of kits in the loop
               nk_max = min(nk,in_k_max)         ! Cannot be greater than nk
               nk_max = max(1,nk_max)            ! Cannot be smaller than 1               
              
               start = in_start       ! Starting point generation procedure
               iprint = in_print      ! Print option
              
               ! The default print is used if user specifieed value is not acceptable
               IF (iprint < -1 .OR. iprint > 4) THEN
                  iprint = 1_c_int
               END IF   

               ! The default start is used if user specifieed value is not acceptable
               IF ((start> 4) .OR. (start <= 0)) THEN
                  start = 2_c_int
               END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5_c_int
               
               mrounds = in_mrounds            ! maximum number of rounds during one 'main iterations'
               mit = in_mit                    ! maximum number of 'main iteration'
               mrounds_clarke = in_mrounds_esc ! maximum number of rounds during one 'Clarke stationary' algorithm
                    
               agg_used = .TRUE.         ! Aggregation is used
               stepsize_used = .FALSE.   ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.     ! The scaling of data is used
            
               !Values for parameter rho used in penalization problem 
               !mrho = (/0.1_c_double, 0.2_c_double, 0.5_c_double, 1.0_c_double &
               !      & 2.0_c_double, 5.0_c_double, 10.0_c_double, 20.0_c_double /) 
               mrho = (/0.2_c_double, 0.5_c_double, 1.0_c_double, 2.0_c_double, &
                     & 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double /) 

!**
               IF (solver==2) THEN		!LMBM		
			      mrho = (/0.5_c_double, 1.0_c_double, 2.0_c_double, &
                         & 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double, 100.0_c_double /) 
               END IF
!**

                ! Problem 
                problem1 = 3_c_int
                problem2 = 3_c_int                
               
                user_n = nft   ! The dimension of the problem

                ed_sol_in_pen = .FALSE.
               !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem      

                tol_zero = (10.0_c_double)**(-6)               
                
       
               IF ( (iprint > 0) .OR. (iprint == -1) ) THEN 
                 IF (ed_sol_in_pen) THEN 
                  WRITE(*,*) 'When the penalized problems is solved we utilize&
                              & the previous solution as a starting point.'
                 ELSE
                  WRITE(*,*) 'When the penalized problems is solved we do NOT utilize&
                              & the previous solution as a starting point.'
                 END IF
                
                 IF ( start > 0 ) THEN 
                  WRITE(*,*) 'Problems are solved from smallest to largest.'
                 ELSE   
                  WRITE(*,*) 'Problems are solved from largest to smallest.'
                 END IF
               
                 WRITE(*,*) 'Starting points are generated with the procedure', start 
                 WRITE(*,*) 'The rho values in the penalized problem:', mrho
               
               END IF                 
               
              !--------------------------------------------------------------------------------------         
          
              ! The starting point
        
                  x_0 = 0.0_c_double
                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------

              !---------------------------------------------------------------------------             
              ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       READING DATA MATRICES
              !---------------------------------------------------------------------------
          
                ! OPEN(78,file=infileX,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) (in_mX(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileY,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) (in_mY(i,j),j=1,2)
                ! END DO
                ! CLOSE(78)        

                ! OPEN(78,file=infileK,status='old',form='formatted')      
                ! DO i=1,nkits
                   ! READ(78,*) (in_mK(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileC,status='old',form='formatted')
                ! DO i=1,nkits
                   ! READ(78,*) in_mC(i)
                ! END DO
                ! CLOSE(78)  
                
              !--------------------------------------------------------------------------
              ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------

 
              !---------------------------------------------------------------------------            
              ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       POPULATE MATRICES
              !--------------------------------------------------------------------------- 
              
               ! Populate the input matrix 'in_mX' of dim {nrecord,nft}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nrecord
                      ind = ind +1
                      in_mX(i,j) = x(ind)
                    END DO
                 END DO
                 
               ! Populate response matrix 'in_mY' of dim {nrecord,2}
                 ind = 0
                 DO j = 1, 2
                    DO i = 1, nrecord
                      ind = ind + 1
                      in_mY(i,j) = y(ind)
                    END DO
                 END DO
                 
              ! Populate kit structure matrix 'in_mK' of dim {nkits,nft}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nk
                      ind = ind + 1
                      in_mK(i,j) = kits(ind)
                    END DO
                 END DO
                 
              ! Populate the cost vector for kits 'in_mC' of dim {nkits}
                 ind = 0
                 DO i = 1, nk
                    ind = ind + 1
                    in_mC(i) = costs(ind)
                 END DO 
                 
              !--------------------------------------------------------------------------
              !^^^^ END: IF Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------
  
                ! The number of value ones in columns of kit matrix
                KitOnes = 0_c_int
                DO j = 1, ncol
                   DO i = 1, nk
                       KitOnes(j) = KitOnes(j) + in_mK(i,j)
                   END DO
                END DO
          
              !---------------------------------------------------------------------------
              !               ORDERING DATA IN INCREASING ORDER BASED ON time 
              !---------------------------------------------------------------------------
           
               DO i = 1, nrecord
                   mTimes(i) = in_mY(i,1)*1.0_c_double
                   mTimesInd(i) = i
               END DO        
            
               CALL heapsort_ind(mTimes,mTimesInd)      

               
               !Transpose of matrix mX
               DO i = 1, nrecord
                 ind = mTimesInd(i)
                 l2 = 1
                 DO j = 1, ncol
                   DO l = 1, KitOnes(j) 
                       mXt(l2,i) = in_mX(ind,j)
                       l2 = l2 + 1
                   END DO   
                 END DO
               END DO
            
               !Transpose of matrix mY          
               DO i = 1, nrecord
                 ind = mTimesInd(i)        
                 DO j = 1, 2
                    mYt(j,i) = in_mY(ind,j)
                 END DO
               END DO    

               ! Alters the kit matrix
               mK = 0_c_int
               l2 = 1_c_int
               DO j = 1, ncol 
                  DO i = 1, nk
                     IF (in_mK(i,j)==1_c_int) THEN  
                        mK(i,l2) = in_mK(i,j)
                        l2 = l2 + 1_c_int
                     END IF
                  END DO               
               END DO 
               
      
               ! Allocation of sizes of matrices in function.f95
               CALL allocate_data_cox(info,nft,nrecord,nk,user_n)   
                           
              !---------------------------------------------------------------------------
              !                     STORING DATA MATRICES 
              !---------------------------------------------------------------------------
              ! Notice: Matrices are transposed! => Each column presents either an observation or a kit!
          
              CALL allocate_matrices_cox(info, mXt, mYt, mK, in_mC,  & 
                                      & nrecord, nft, nk)

              ! Scaling of data 
              IF (scale_in_use) THEN
                  CALL scaling_cox(info)
                  
                 ! Original data matrices mXt and mYt are changed to the scaled ones 
                  CALL scaling_matrix_cox(info, mXt)
                  
              END IF  

              
              ! The initialization of beta vector
              beta_nft = 0.0_c_double
              beta = 0.0_c_double
              
              ! The print in DBDC method is supressed
              iprint_DBDC = 0
              
              !-------------------------------------------------------------------------------------------------
              ! The best beta_solution for Cox's proportional hazard model without regularization/penalization                    
                      
              CALL set_k(info, nkits)           ! All kits can be used

!**                      
			  IF (solver == 1) THEN	  
                  CALL DBDC_algorithm( f_solution_DBDC, x_koe, x_0, 0.0_c_double, 0.0_c_double, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, 3_c_int, 3_c_int, user_n, &
                            & info)            
			  ELSE IF (solver == 2) THEN 
			  
                  CALL allocate_xn(user_n) 
                  CALL init_LMBMinfo(problem1, info) 
                  CALL init_x(x_0)          
                  CALL set_rho_LMBM(0.0_c_double)                 
                  CALL set_lambda_LMBM(0.0_c_double)                  
                  CALL lmbm(mc,f_solution_DBDC,iout(1),iout(2),iout(3),iout(4),LMBMstart)     
                  CALL copy_x(x_koe)
                  			  
              END IF			  
!**
                            
              ! Notice: * solution x_koe is obtained by fitting Cox's model to data without regularization
              !         * x_koe is utilized in formation of starting points 

              !-------------------------------------------------------------------------------------------------              
              
               x_ed = 0.01_c_double    ! We initialize the previous solution, since we do not have a value for it 
               
               
               IF ((iprint > 0) .OR. (iprint == -1)) THEN
                 WRITE(*,*) 
                 WRITE(*,*) 'Value of Cox proportional hazard model without regularization:', f_solution_DBDC
                 WRITE(*,*) 
               END IF 
               
               IF (iprint==1) THEN
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
                 WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
               END IF              
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
               
                 kit_num_ed = 0           ! The number of kits in the previous solution
                 kits_beta_ed = 0         ! The kits in the previous solution
                 small = (10_c_double)**9 ! There is no previous solution known
                 
                 IF (start == 1) THEN
                    nstart = 1
                 ELSE IF (start == 2 .OR. start == 3) THEN
                    nstart = nk
                 ELSE IF (start == 4) THEN 
                    nstart = start_max        ! The number of random starting points
                 END IF 
                 
                 !---------------------------------------------------------------------------------------------
                 ! Solves in a loop all L0-norm problems with the fixed number of k=1,..., nk_max of kits    
                 !---------------------------------------------------------------------------------------------
                 DO k = 1, nk_max              ! In this loop we calculate the solution for problem 3 wiht L0-norm such that the number of kits varies from 1 to nk_max
                  
                  IF (iprint>=2 .OR. iprint == -1) THEN
                       WRITE(*,*) '-------------------------------------' 
                       WRITE(*,*) 'PROBLEM with', k, 'kits' 
                       
                  END IF
                  
                  CALL set_k(info, k)            ! The number of nonzero kits is fixed 
                  
                  jaannos = MOD(nk,tread_num)    ! The residue in the division
                  kokoosa = nk/tread_num         ! The integer part obtained from the division
                  
                  mPrNum = kokoosa               ! The initialization of the number of problems for threads
                  DO i = 1, jaannos              ! We increase the number of problems with one for specific threads
                     mPrNum(i) = mPrNum(i)+1                  
                  END DO

                  !----------------------------------------------------------
                  ! Starting points for problem with fixed number of k kits 
                  !----------------------------------------------------------
                  
                  !$OMP PARALLEL DO PRIVATE(x_solution, f_solution, i) & 
                  !$OMP FIRSTPRIVATE(x_ed, x_koe, kits_beta_ed, kit_num_ed)  &             
                  !$OMP FIRSTPRIVATE(small, k, iprint, mrho)                 &  
                  !$OMP FIRSTPRIVATE(mit, mrounds, mrounds_clarke)           &    
                  !$OMP FIRSTPRIVATE(agg_used, stepsize_used)                &
                  !$OMP FIRSTPRIVATE(problem1, problem2)                     &  
                  !$OMP FIRSTPRIVATE(in_b1, in_b2, in_m, in_c, in_r_dec)     &  
                  !$OMP FIRSTPRIVATE(in_r_inc, in_eps1, in_b, in_m_clarke)   &  
                  !$OMP FIRSTPRIVATE(in_eps, in_crit_tol)                    &  
                  !$OMP FIRSTPRIVATE(nft, nrecord, nk, user_n)               &  
                  !$OMP FIRSTPRIVATE(mXt, mYt, mK, in_mC)                    &  
                  !$OMP SHARED(points, f_points, mPrNum)               
                  
                   DO i = 1, tread_num              ! Different starting points in threads are looked through to solve the problem 3 with fixed number of nonzero kits                   
                    
                      CALL solution_with_k_kits_cox( x_solution, f_solution, k, i, small, &
                                               & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                               & nk, start, iprint, mrho, mit, mrounds, mrounds_clarke, &
                                               & agg_used, stepsize_used, user_n, problem1, problem2, &
                                               & mXt, mYt, mK, in_mC, nrecord,  & 
                                               & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                               & in_b, in_m_clarke, in_eps, in_crit_tol, mPrNum, solver)
                    
                    ! Storing of the obtained solutions and the corresponding objective values
                     !$OMP CRITICAL 
                     startind = 1
                     IF (i > 1) THEN    
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                     END IF                        
                        
                     DO j3 = 1, mPrNum(i)   
                        DO j = 1, user_n
                            points(j,startind+(j3-1)) = x_solution(j,startind+(j3-1))
                        END DO                            
                         f_points(startind+j3-1) = f_solution(startind+j3-1)
                     END DO  
                     !$OMP END CRITICAL                  
                   
                   END DO          
                  !$OMP END PARALLEL DO   
                   
                   IF (start == 1) THEN 
                      small = f_points(1)
                      min_ind = 1
                   END IF 
                   
                   IF (start == 2 .OR. start == 3) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0                     ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nk                    ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF

                   IF (start == 4) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0                     ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nstart                ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF        
                   
                   ind = (k-1)*user_n  
                   DO i = 1, user_n
                     x_ed(i) = points(i,min_ind)            ! The beta solution yielding the smallest objective function value is set as the previous solution
                     beta_nft(ind+i) = points(i,min_ind)        ! The beta solution yielding the smallest objective function value is stored to solution vector 
                   END DO
              
                   help_counter = 0           
                   DO j = 1, nft         ! The number of zero components in the previous beta vector is computed
                      IF ( ABS(x_ed(j)) < tol_zero ) THEN
                        help_counter = help_counter+1
                      END IF
                   END DO

                   cost = 0.0_c_double         ! The initialization of the previous solution
                   kit_num_ed = 0              ! The initialization of the number of kits in the previous solution
                   kits_beta_ed = 0            ! The initialization of kits in the previous solution                    

                   !Calculation of kits in the previous solution
                   DO j1 = 1, nk    ! Each kit is looked through
                     kit_in_use = .FALSE.
                     DO j2 = 1, nft
                       IF (info%mk(j2,j1)==1) THEN 
                         IF ( ABS(x_ed(j2)) >= tol_zero ) THEN 
                           kit_in_use = .TRUE.           
                         END IF
                       END IF
                     END DO
                     IF (kit_in_use) THEN              ! Executed if kit j1 is in the previous solution
                       cost = cost + info%mC(j1)            ! The cost of kit j1 is taken into account
                       kit_num_ed = kit_num_ed + 1     ! The number of kits is updated
                       kits_beta_ed(kit_num_ed) = j1   ! The index of kit j1 is updated to table kits_beta                 
                     END IF 
                   END DO       

                   ! The storing of kits to "betakits" vector
                   ind = (k-1)*nk
                   DO i = 1, nk
                      betakits(ind+i) = kits_beta_ed(i)                
                   END DO

 
                   IF ((iprint>=2) .OR. (iprint == -1)) THEN 
                     WRITE(*,*)
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*) 'Information about the best solution:'
                     WRITE(*,*)                  
                     WRITE(*,*)  '  f=',small 
                     WRITE(*,*)  '  number of zero elements=',help_counter
                     WRITE(*,*)  '  number of nonzero elements=',nft-help_counter 
                     WRITE(*,*)  '  cost=',cost
                     WRITE(*,*)  '  number of kits=',kit_num_ed 
                     WRITE(*,*)  '  kits=',kits_beta_ed(1:kit_num_ed)    
                     WRITE(*,*)                  
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*)
                     
                   END IF
 
                   IF (iprint==1) THEN
                     WRITE(*,*) small, help_counter, nft-help_counter, &
                     & cost, kit_num_ed, kits_beta_ed(1:kit_num_ed)                
                   END IF
                   
                 END DO
              
   
              END IF              
              
              !--------------------------------------------------------------------------
              !                 SOLVING OF L0-NORM PROBLEM ENDS 
              !---------------------------------------------------------------------------  

            
              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk 
              !---------------------------------------------------------------------------                

              IF (scale_in_use) THEN  ! Rescaling
                 
                 CALL set_k(info, nk)               
                 info%user_lambda = 0.0_c_double
                 CALL rescaling_cox(info)        ! the rescaling of data 
                 
                 DO k = 1, nk         ! Each solution is rescaled
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   CALL rescaling_beta_cox(info, beta_solution)         ! Rescaling of solution
                   f1_current = f1(info, beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info, beta_solution,problem1,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current               ! The objective function value for problem 3 with k nonzero kits 
                    WRITE(*,*) 'f',k,':', fperk(k)
                    
                   l2 = 1   
                   ind = (k-1)*ncol
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO       
                   
                 END DO 
                 
              ELSE   
                  CALL set_k(info, nk)               
                  info%user_lambda = 0.0_c_double
                 
                 DO k = 1, nk         
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   f1_current = f1(info,beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem1,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current               ! The objective function value for problem 3 with k nonzero kits 
                    WRITE(*,*) 'f',k,':', fperk(k)
                    
                   l2 = 1   
                   ind = (k-1)*ncol
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO     
                   
                 END DO            
              END IF            

              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk COMPLETED
              !---------------------------------------------------------------------------  
               
             CALL cpu_time(f_time)   ! Finish CPU timing    
             cpu = f_time-s_time             
             
             CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
             ! Calculate the elapsed 'clock' time in seconds:
             elapsed_time=(1.0_c_double*clock_end-clock_start)/clock_rate                
             
             
             IF ((iprint >= 1) .OR. (iprint == -1)) THEN               
                 WRITE(*,*) 'Used CPU:', cpu   , 'Elapsed time:', elapsed_time             
             END IF
       
             CALL deallocate_data_cox(info) 
             CALL deallocate_x()

         END SUBROUTINE oscar_cox       
        !--------------------------------------------------------------------------------------  



        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |               THE DOUBLE BUNDLE ALGORITHM FOR MEAN SQUARE ERROR                |  |
        !  |                                                                                |  |
        !  |                     SOLUTION FOR EVERY NUMBER OF KITS                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        ! Subroutine for oscar with the mean square error model
        
    !--------------------------------------------------------------------------
    ! ^^^^ START: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------        
          SUBROUTINE oscar_mse(x, y, kits, costs, nrow, ncol, nkits, beta, fperk, &
       & in_print, in_start, in_k_max, &
       & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
       & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, &
       & nKitOnes, betakits, &
       & solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
       & in_tolg, in_tolg2, in_eta, in_epsL ) &      
        & BIND(C, name = "oscar_mse_f_")               
    !--------------------------------------------------------------------------
    ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------


          !--------------------------------------------------------------------------            
          ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------    
          ! SUBROUTINE oscar_mse(infileX, infileY, infileK, infileC,  &
          !                     & nrow, ncol, nkits, &
          !                     & beta, fperk, in_print, in_start, in_k_max, &
          !                     & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
          !                     & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, &
          !                     & nKitOnes, betakits, solver_id, & 
          !                     & in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
          !                     & in_tolg, in_tolg2, in_eta, in_epsL)     
          !--------------------------------------------------------------------------            
          ! ^^^^ END: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------   
            
            !_____________________________________________________________________________________
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            !
            !  !! NOTICE !!  The DATA is SCALED during this SUBROUTINE !
            !
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            !
            !
            ! * 'problem = 4'        : the used objective function is mean square error with L0-norm for kits
            !                          (parametric approach since the number of nonzero elements is fixed) 
            !
            !
            !   Solves in a loop all L0-norm problems with the fixed number of k=1,...,in_k_max of kits    (if in_start > 0)
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'ncol'        : The dimension of the problem = the number of features in a predictor, INTEGER
            !         * 'nrow'        : The number of records (data points), INTEGER
            !         * 'nkits'       : the (maximum) number of kits, INTEGER
            !
            !         * 'in_print'    : specifies the print, INTEGER
            !         * 'in_start'    : specifies how starting points are selected when the L0-norm problem is solved for  
            !                           the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'in_k_max'    : specifies how many kits are used in the last problem of the loop  (NEEDS to be 1 <= 'in_k_max' <= nkits !)
            !
            !         * 'nKitOnes'    : the number of value 1 in kit matrix
			!
            !         * 'solver_id'   : the used solver (1=DBDC and 2=LMBM) INTEGER
            !
            !         PARAMETERS in DBDC method:
            !
            !         * in_mrounds     : the maximum number of rounds in one main iteratioin
            !         * in_mit         : the maximum number of main iterations
            !         * in_mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance          
            !
            !         PARAMETERS in LMBM method:
            !
            !         * in_na        : Size of the bundle na >= 2.
            !         * in_mcu       : Upper limit for maximum number of stored corrections, mcu >= 3.
            !         * in_mcinit    : Initial maximum number of stored corrections, mcu >= mcinit >= 3.
            !                          If mcinit <= 0, the default value mcinit = 3 will be used. 
            !                          However, the value mcinit = 7 is recommented.             
            !         * in_tolf      : Tolerance for change of function values (default = 1.0E-5).
            !         * in_tolf2     : Second tolerance for change of function values.
            !                          - If tolf2 < 0 the the parameter and the corresponding termination 
            !                            criterion will be ignored. 
            !                          - If tolf2 = 0 the default value 1.0E+4 will be used
            !         * in_tolg      : Tolerance for the first termination criterion (default = 1.0E-5).   
            !         * in_tolg2     : Tolerance for the second termination criterion (default = tolg). clustering code small data  
            !         * in_eta       : Distance measure parameter, eta >= 0.   
            !         * in_epsL      : Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.)                
            !
            !         NOTICE: DATA IS GIVEN IN VECTOR FORMAT ! Due to this values for observations/kits are given in these vector one after another
            !
            !         * 'x'          : the matrix X of input variables in vector form --> will be stored to in_mX (each row is one observation)
            !                               REAL, DIMENSION(nrow*ncol)  
            !         * 'y'          : the vector Y of outputs --> will be stored to in_mY 
            !                               REAL, DIMENSION(nrow) 
            !         * 'costs'      : the cost vector C for kits --> will be stored to in_mC
            !                               REAL, DIMENSION(knits)
            !         * 'kits'       : the matrix K of kit structures in vector form (values are either 0 or 1) --> will be stored to in_mK (each row is one kit)
            !                               INTEGER, DIMENSION(knits*ncol)
            !
            !
            ! OUTPUT: * 'fperk'      : The vector containing objective function values at for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(nkits) 
            !         * 'beta'       : The vector containing solution vectors beta obtained for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION((ncol+1)*nkits)    (in eah solution first 'ncol' values are beta_1,..,beta_nft and the last one is beta_0)
            !
            ! NOTICE: * The dimension of vectors 'fperk' has to be 'nkits'         ('nkits' is the maximum number of kits)
            !         * The dimension of vectors 'beta' has to be '(ncol+1)*nkits' ('ncol+1' is the dimension of the problem and 'nkits' is the number of kits)
            !
            ! OPTIONAL INPUT (CAN BE INCLUDED AS INPUT PARAMETERS IF NEEDED. AT THE MOMENT DEFAULT VALUES ARE USED FOR THEM.):
            !         * 'mit'            : The maximum number of 'main iterations', INTEGER                  
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration', INTEGER
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Escape procedure', INTEGER
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm, LOGICAL           
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm, LOGICAL
            !         * 'scale_in_use'   : If .TRUE. then the data is scaled
            !         * 'CPUtime'       : the CPU time (in seconds) REAL
            !
            ! NOTICE: * 'in_print' has to be 0, 1, 2 or 3. If it is NOT then DEFAULT value 1 is used.             
            !         * 'in_start' has to be 1, 2, 3 or 4. If it is NOT then DEFAULT value 2 is used.             
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! ! INPUTs
               ! INTEGER(KIND = c_int), INTENT(IN) :: nrow                  ! Number of rows in x (i.e. records)
               ! INTEGER(KIND = c_int), INTENT(IN) :: ncol                  ! Number of cols in x (i.e. features)
               ! INTEGER(KIND = c_int), INTENT(IN) :: nkits                 ! Number of kits for features
               
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_start              ! Starting point procedure used
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_print              ! Print used
               ! INTEGER(KIND = c_int), INTENT(IN) :: in_k_max              ! The maximum number of kits in the loop
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds              ! the maximum number of rounds in one main iteration
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mit                  ! the maximum number of main iterations
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               ! REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               ! REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               ! REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance

               ! INTEGER(KIND=c_int), INTENT(IN) :: nKitOnes                ! the number of value 1 in kit matrix
			   ! INTEGER(KIND=c_int), INTENT(IN) :: solver_id               ! the used solver (1=DBDC and 2=LMBM)
!**
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_na                  ! Size of the bundle na >= 2.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcu                 ! Upper limit for maximum number of stored corrections, mcu >= 3.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcinit              ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! ! However, the value mcinit = 7 is recommented.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf                ! Tolerance for change of function values (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf2               ! Second tolerance for change of function values.
                                                                          ! !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          ! !   criterion will be ignored. 
                                                                          ! !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg                ! Tolerance for the first termination criterion (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg2               ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               ! REAL(KIND=c_double), INTENT(IN) ::  in_eta                 ! Distance measure parameter, eta >= 0.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_epsL                ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**		

              ! INPUTs
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: nrow                  ! Number of rows in x (i.e. records)
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: ncol                  ! Number of cols in x (i.e. features)
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: nkits                 ! Number of kits for features
               
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_start              ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_print              ! Print used
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_k_max              ! The maximum number of kits in the loop
               
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mrounds              ! the maximum number of rounds in one main iteration
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mit                  ! the maximum number of main iterations
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_crit_tol             ! the stopping tolerance: criticality tolerance

               INTEGER(KIND=c_int), INTENT(IN), VALUE :: nKitOnes                ! the number of value 1 in kit matrix

			   INTEGER(KIND=c_int), INTENT(IN), VALUE :: solver_id               ! the used solver (1=DBDC and 2=LMBM)

 !**
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_na           ! Size of the bundle na >= 2.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcu          ! Upper limit for maximum number of stored corrections, mcu >= 3.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcinit       ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! However, the value mcinit = 7 is recommented.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf         ! Tolerance for change of function values (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf2        ! Second tolerance for change of function values.
                                                                          !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          !   criterion will be ignored. 
                                                                          !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg         ! Tolerance for the first termination criterion (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg2        ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_eta          ! Distance measure parameter, eta >= 0.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_epsL         ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**

            !--------------------------------------------------------------------------
            ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nrow*ncol)  :: x        ! Vector of data values
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nrow)       :: y        ! Vector of outputs
               INTEGER(KIND = c_int), INTENT(IN), DIMENSION(nkits*ncol) :: kits     ! Vector of kit indicator values (binary indicators)
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nkits)      :: costs    ! Costs associated to each kit
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------


            !--------------------------------------------------------------------------            
            ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
            !--------------------------------------------------------------------------
            !   CHARACTER(LEN=80), INTENT(IN) :: infileX         ! The name of "predictor matrix" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileY         ! The name of "observed time and label" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileC         ! The name of "kit costs" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileK         ! The name of "predictor matrix" file   
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^          
            !--------------------------------------------------------------------------
            
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION((ncol+1)*nkits)  :: beta       !Output variable for beta coefficients per k
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nkits)           :: fperk      !Output variable target function value per k     
               INTEGER(KIND = c_int), INTENT(OUT), DIMENSION(nkits*nkits)     :: betakits   !Output variable telling kits in beta coefficients per k
               
           
           !***************************** LOCAL VARIABLES ************************************      
 
               TYPE(set_info) :: info                         ! The set of information                     

               REAL(KIND=c_double) :: CPUtime                 ! the CPU time (in seconds)

               INTEGER(KIND=c_int) :: solver                  ! defines the solver used (1=DBDC, 2=LMBM)			   

               INTEGER(KIND=c_int) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int) :: nk                      ! the number of kits 
               INTEGER(KIND=c_int) :: nk_max                  ! the maximum number of kits in the loop 
   
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) :: beta_solution    ! the solution vector beta obtained for the problem

               REAL(KIND=c_double), DIMENSION(nKitOnes+1,nkits) :: points     ! the beta_solutions for problem 4 for fixed k ('nkits' different starting points) 
               REAL(KIND=c_double), DIMENSION(nkits) ::      f_points     ! the objective function values for problem 4 for fixed k ('nkits' different starting points)
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) ::       x_koe      ! The solution to  problem 4  without regularization
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller
               
               REAL(KIND=c_double), DIMENSION(nKitOnes+1,nkits) ::  x_solution      ! the beta solutions for the problem with fixed starting points and number of kits
               REAL(KIND=c_double), DIMENSION(nkits) :: f_solution              ! the objective function values at the solution 'x_solution'
               REAL(KIND=c_double) :: f_solution_DBDC                           ! the f_solution obtained from DBDC method
        
               REAL(KIND=c_double), DIMENSION(nrow,ncol) :: in_mX     ! predictor matrix (row is an observation)
               REAL(KIND=c_double), DIMENSION(nrow) :: in_mY          ! output values 
               INTEGER(KIND=c_int), DIMENSION(nkits,ncol) :: in_mK    ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nkits) :: in_mC         ! kit costs                   
              
               REAL(KIND=c_double), DIMENSION(nKitOnes,nrow) :: mXt       ! predictor matrix (column is an observation)
               REAL(KIND=c_double), DIMENSION(nrow) :: mYt                 ! Outputs
               INTEGER(KIND=c_int), DIMENSION(nkits,nKitOnes) :: mK        ! Modified kit matrix (There is only one value one in each column)  
               REAL(KIND = c_double), DIMENSION((nKitOnes+1)*nkits)  :: beta_nft   !Output variable for beta coefficients per k

               INTEGER(KIND=c_int), DIMENSION(ncol) :: KitOnes        ! The number of kits where each variables is located. If kits do not intersect then each variable is only in one kit!
               
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta     ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'

               REAL(KIND=c_double), DIMENSION(nKitOnes+1) :: x_0          ! the starting point
                             
               REAL(KIND=c_double), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=c_int) :: nstart              ! the number of starting point
               INTEGER(KIND=c_int) :: start_max           ! the number of starting point when 'start = 5'

               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                       
               INTEGER(KIND=c_int) :: mit                ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds            ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds_clarke     ! the maximum number of rounds during one 'Clarke stationary' algorithm
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

               ! Possible USER PARAMETER
               INTEGER(KIND=c_int) :: iprint   ! specifies the print
                                               !   iprint = 0 : print is suppressed
                                               !   iprint = -1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Not table form)
                                               !   iprint = 1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Table form)
                                               !   iprint = 2 : prints for each L0-norm problem the final value of Cox's proportional hazard model obtained from each starting point
                                               !   iprint = 3 : prints for each starting point and L0-norm problem all intermediate results together the final results
 
               INTEGER(KIND=c_int) :: start    !   start = 1  : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                               !                (starting point is the solution obtained for L0-norm problem with k-1 nonzero components)
                                               !
                                               !   start = -1 : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1
                                               !                (starting point is the solution obtained for L0-norm problem with k+1 nonzero components)  
                                               !
                                               !   start = 2  : for L0-norm problem with k nonzero components uses 'nkits-k+1' starting points
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
               REAL(KIND=c_double) :: small           ! The cost of solution beta
               
               REAL(KIND=c_double) :: s_time, f_time  ! The start and finish times
               REAL(KIND=c_double) :: cpu             ! The cpu time
               
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                     ! Random number
               REAL(KIND=c_double), DIMENSION(nkits) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=c_int) :: problem1              ! The DC component f1 
               INTEGER(KIND=c_int) :: problem2              ! The DC component f2 
        
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: nremoved
               INTEGER(KIND=c_int) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=c_int) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii
               INTEGER(KIND=c_int) :: l, l2
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization
               
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mPrNum   ! The number of problems for threads
               INTEGER(KIND=c_int) :: tread_num                  ! The number of thread
               INTEGER(KIND=c_int) :: jaannos                    ! The residue in division
               INTEGER(KIND=c_int) :: kokoosa                    ! The integer obtained in division
               INTEGER(KIND=c_int) :: startind, j3               ! Help variables

               REAL(KIND=c_double) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER(KIND=c_int) :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time 
               
!**			   
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.               
 !**
 
               CALL cpu_time(s_time)   ! Start CPU timing
               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time     

!**
              ! The solver is defined/selected
               IF (solver_id == 1) THEN
                   solver = 1 !DBDC
               ELSE IF (solver_id == 2) THEN
                   solver = 2 !LMBM
               ELSE 
                   solver = 1 !DBDC              
               END IF
!**
  
             ! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
               
               ! Maximum number of possible treads in parallellization
               max_threads = omp_get_max_threads()
			   !max_threads = 1
               
!**
			   IF (solver == 2) THEN  !No parallellization with LMBM can be used
                   max_threads = 1
			   END If	   
!**        
               CALL omp_set_num_threads(max_threads)
               
               
               tread_num = max_threads

               WRITE(*,*) 'The number of threads', max_threads             
               ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----   

               ! The initialization of parametrs used in DBDC methods
               CALL allocate_parameters(info, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)

!**                                          
                 ! The initialization of parameters used in LMBM method
                 CALL init_par(in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, &
                                & in_tolg2, in_eta, in_epsL)
!**
               
               ! Set the number of rows and columns inside Fortran  + kits           
               nrecord = nrow
               nft = nKitOnes
               nk = nkits
               
               ! The maximum number of kits in the loop
               nk_max = min(nk,in_k_max)               ! Cannot be greater than nk
               nk_max = max(1_c_int,nk_max)            ! Cannot be smaller than 1               
               
               start = in_start       ! Starting point generation procedure
               iprint = in_print      ! Print option
              
               ! The default print is used if user specifieed value is not acceptable
               IF (iprint < -1 .OR. iprint > 4) THEN
                  iprint = 1_c_int
               END IF   

               ! The default start is used if user specifieed value is not acceptable
               IF ((start > 4) .OR. (start <= 0)) THEN
                  start = 2_c_int
               END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5_c_int
               
               mrounds = in_mrounds            ! maximum number of rounds during one 'main iterations'
               mit = in_mit                    ! maximum number of 'main iteration'
               mrounds_clarke = in_mrounds_esc ! maximum number of rounds during one 'Clarke stationary' algorithm
          
               iprint_DBDC = 0_c_int           ! basic print of intermediate results and extended print of final results
          
               agg_used = .TRUE.         ! Aggregation is used
               stepsize_used = .FALSE.   ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.     ! The scaling of data is used
            
               ed_sol_in_pen = .FALSE.
               !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem

               !Values for parameter rho used in penalization problem 
               !mrho = (/0.1_c_double, 0.2_c_double, 0.5_c_double, 1.0_c_double &
               !      & 2.0_c_double, 5.0_c_double, 10.0_c_double, 20.0_c_double /) 
               mrho = (/0.2_c_double, 0.5_c_double, 1.0_c_double, 2.0_c_double, &
                     & 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double /) 

               ! Problem 
               problem1 = 4_c_int
               problem2 = 4_c_int                
               
               user_n = nft+1
                
               tol_zero = (10.0_c_double)**(-6)
       
               IF ( (iprint > 0) .OR. (iprint == -1) ) THEN 
                 IF (ed_sol_in_pen) THEN 
                  WRITE(*,*) 'When the penalized problems is solved we utilize&
                              & the previous solution as a starting point.'
                 ELSE
                  WRITE(*,*) 'When the penalized problems is solved we do NOT utilize&
                              & the previous solution as a starting point.'
                 END IF
                
                 IF ( start > 0 ) THEN 
                  WRITE(*,*) 'Problems are solved from smallest to largest.'
                 ELSE   
                  WRITE(*,*) 'Problems are solved from largest to smallest.'
                 END IF
               
                 WRITE(*,*) 'Starting points are generated with the procedure', start 
                 WRITE(*,*) 'The rho values in the penalized problem:', mrho
               
               END IF                
               
              !--------------------------------------------------------------------------------------
          
              ! The starting point
        
                  x_0 = 0.0_c_double
                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------

              !---------------------------------------------------------------------------             
              ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       READING DATA MATRICES
              !---------------------------------------------------------------------------
          
                ! OPEN(78,file=infileX,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) (in_mX(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileY,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) in_mY(i)
                ! END DO
                ! CLOSE(78)        

                ! OPEN(78,file=infileK,status='old',form='formatted')      
                ! DO i=1,nk
                   ! READ(78,*) (in_mK(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileC,status='old',form='formatted')
                ! DO i=1,nk
                   ! READ(78,*) in_mC(i)
                ! END DO
                ! CLOSE(78)  
                
              !--------------------------------------------------------------------------
              ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------

 
              !---------------------------------------------------------------------------            
              ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       POPULATE MATRICES
              !--------------------------------------------------------------------------- 
              
               ! Populate the input matrix 'in_mX' of dim {nft,nrecord}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nrecord
                      ind = ind + 1
                      in_mX(i,j) = x(ind)
                    END DO
                 END DO
                 
               ! Populate response matrix 'in_mY' of dim {nrecord}
                 ind = 0
                 DO i = 1, nrecord
                      ind = ind + 1
                      in_mY(i) = y(ind)
                 END DO
                 
              ! Populate kit matrix 'in_mK' of dim {nkits,nft}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nk
                      ind = ind + 1
                      in_mK(i,j) = kits(ind)
                    END DO
                 END DO
                 
              ! Populate the cost vector 'in_mC' for kits of dim {nkits}
                 ind = 0
                 DO i = 1, nk
                    ind = ind + 1
                    in_mC(i) = costs(ind)
                 END DO 
                 
              !--------------------------------------------------------------------------
              !^^^^ END: IF Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------

           
              !---------------------------------------------------------------------------
              !               ALLOCATION OF DATA MATRICES 
              !---------------------------------------------------------------------------
                   
               ! Allocation of data matrices in function.f95
               CALL allocate_data_mse(info,nft,nrecord,nk,nk)   
               
                ! The number of value ones in columns of kit matrix
                KitOnes = 0_c_int
                DO j = 1, ncol
                   DO i = 1, nk
                       KitOnes(j) = KitOnes(j) + in_mK(i,j)
                   END DO
                END DO

               
              !---------------------------------------------------------------------------
              !                     STORING DATA MATRICES 
              !---------------------------------------------------------------------------
              ! Notice: Matrices are transposed! => Each column presents either an observation or a kit!
              
              
              !Transpose of matrix mX
               DO i = 1, nrecord
                 l2 = 1
                 DO j = 1, ncol
                   DO l = 1, KitOnes(j) 
                       mXt(l2,i) = in_mX(i,j)
                       l2 = l2 + 1
                   END DO   
                 END DO
               END DO             
            
               !Transpose of matrix mY          
               DO i = 1, nrecord
                    mYt(i) = in_mY(i)
               END DO    

               ! Alters the kit matrix
               mK = 0_c_int
               l2 = 1_c_int
               DO j = 1, ncol 
                  DO i = 1, nk
                     IF (in_mK(i,j)==1_c_int) THEN  
                        mK(i,l2) = in_mK(i,j)
                        l2 = l2 + 1_c_int
                     END IF
                  END DO               
               END DO 
               
              CALL allocate_matrices_mse(info, mXt, mYt, mK, in_mC,  & 
                                      & nrecord, nft, nk)

              ! Scaling of data 
              IF (scale_in_use) THEN
                  CALL scaling_mse(info)
                  
                 ! Original data matrix mXt is changed to the scaled one 
                  CALL scaling_matrix_mse(info, mXt)
                  
              END IF  
           

              ! The initialization of beta vector
              beta_nft = 0.0_c_double
              beta = 0.0_c_double
              
              ! The best beta_solution for Cox's proportional hazard model without regularization/penalization    
              
              CALL set_k(info, nkits)           ! All kits can be used        

!**                      
			  IF (solver == 1) THEN	  
                  CALL DBDC_algorithm( f_solution_DBDC, x_koe, x_0, 0.0_c_double, 0.0_c_double, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n, &
                            & info)            
			  ELSE IF (solver == 2) THEN 
			  
                  CALL allocate_xn(user_n) 
                  CALL init_LMBMinfo(problem1, info) 
                  CALL init_x(x_0)          
                  CALL set_rho_LMBM(0.0_c_double)                 
                  CALL set_lambda_LMBM(0.0_c_double)                  
                  CALL lmbm(mc,f_solution_DBDC,iout(1),iout(2),iout(3),iout(4),LMBMstart)     
                  CALL copy_x(x_koe)
                  			  
              END IF			  
!**


                            
              ! Notice: * solution x_koe is obtained by fitting Cox's model to data without regularization
              !         * x_koe is utilized in formation of starting points   
              
               x_ed = 0.01_c_double    ! We initialize the previous solution, since we do not have a value for it 
               
              
               IF ((iprint > 0) .OR. (iprint == -1)) THEN
                 WRITE(*,*) 
                 WRITE(*,*) 'Value of mean square error model without regularization:', f_solution_DBDC
                 WRITE(*,*) 
               END IF 
               
               IF (iprint==1) THEN
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
                 WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
               END IF              
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
               
                 kit_num_ed = 0           ! The number of kits in the previous solution
                 kits_beta_ed = 0         ! The kits in the previous solution
                 small = (10_c_double)**9 ! There is no previous solution known
                 
                 IF (start == 1) THEN
                    nstart = 1
                 ELSE IF (start == 2 .OR. start == 3) THEN
                    nstart = nk
                 ELSE IF (start == 4) THEN 
                    nstart = start_max        ! The number of random starting points
                 END IF 
                 
                 !---------------------------------------------------------------------------------------------
                 ! Solves in a loop all L0-norm problems with the fixed number of k=1,..., nk_max of kits    
                 !---------------------------------------------------------------------------------------------
                 DO k = 1, nk_max              ! In this loop we calculate the solution for problem 3 wiht L0-norm such that the number of kits varies from 1 to nk_max
                  
                  IF (iprint>=2 .OR. iprint == -1) THEN
                       WRITE(*,*) '-------------------------------------' 
                       WRITE(*,*) 'PROBLEM with', k, 'kits' 
                       
                  END IF
                  
                  CALL set_k(info, k)            ! The number of nonzero kits is fixed
                                   
                  jaannos = MOD(nk,tread_num)    ! The residue in the division
                  kokoosa = nk/tread_num         ! The integer part obtained from the division
                  
                  mPrNum = kokoosa               ! The initialization of the number of problems for threads
                  DO i = 1, jaannos              ! We increase the number of problems with one for specific threads
                     mPrNum(i) = mPrNum(i)+1                  
                  END DO

                  !----------------------------------------------------------
                  ! Starting points for problem with fixed number of k kits 
                  !----------------------------------------------------------
                  
                  !$OMP PARALLEL DO PRIVATE(x_solution, f_solution, i) & 
                  !$OMP FIRSTPRIVATE(x_ed, x_koe, kits_beta_ed, kit_num_ed)  &             
                  !$OMP FIRSTPRIVATE(small, k, iprint, mrho)                 &  
                  !$OMP FIRSTPRIVATE(mit, mrounds, mrounds_clarke)           &    
                  !$OMP FIRSTPRIVATE(agg_used, stepsize_used)                &
                  !$OMP FIRSTPRIVATE(problem1, problem2)                     &  
                  !$OMP FIRSTPRIVATE(in_b1, in_b2, in_m, in_c, in_r_dec)     &  
                  !$OMP FIRSTPRIVATE(in_r_inc, in_eps1, in_b, in_m_clarke)   &  
                  !$OMP FIRSTPRIVATE(in_eps, in_crit_tol)                    &  
                  !$OMP FIRSTPRIVATE(nft, nrecord, nk, user_n)               &  
                  !$OMP FIRSTPRIVATE(mXt, mYt, mK, in_mC)                    &  
                  !$OMP SHARED(points, f_points, mPrNum)               
                  
                   DO i = 1, tread_num          ! Different starting points for threads are looked through to solve the problem 4 with fixed number of nonzero kits                   
                    
                      CALL solution_with_k_kits_mse(x_solution, f_solution, k, i, small, &
                                               & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                               & nk, start, iprint, mrho, mit, mrounds, mrounds_clarke, &
                                               & agg_used, stepsize_used, nft, problem1, problem2, &
                                               & mXt, mYt, mK, in_mC, nrecord,  & 
                                               & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                               & in_b, in_m_clarke, in_eps, in_crit_tol, mPrNum, solver)
                    
                    ! Storing of the obtained solution and the corresponding objective value
                    !$OMP CRITICAL 
                     startind = 1
                     IF (i > 1) THEN    
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                     END IF                        
                        
                     DO j3 = 1, mPrNum(i)   
                        DO j = 1, user_n
                            points(j,startind+(j3-1)) = x_solution(j,startind+(j3-1))
                        END DO     
                      f_points(startind+j3-1) = f_solution(startind+j3-1)
                     END DO  
                    !$OMP END CRITICAL                    
                        
                   END DO          
                  !$OMP END PARALLEL DO   
                            
                   IF (start == 1) THEN 
                      small = f_points(1)
                      min_ind = 1
                   END IF 
                   
                   IF (start == 2 .OR. start == 3) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nk              ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF

                   IF (start == 4) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nstart          ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF                  

                   ind = (k-1)*user_n  
                   DO i = 1, user_n
                     x_ed(i) = points(i,min_ind)            ! The beta solution yielding the smallest objective function value is set as the previous solution
                     beta_nft(ind+i) = points(i,min_ind)        ! The beta solution yielding the smallest objective function value is stored to solution vector 
                   END DO
                   
                   help_counter = 0           
                   DO j = 1, nft         ! The number of zero components in the previous beta vector is computed
                      IF ( ABS(x_ed(j)) < tol_zero ) THEN
                        help_counter = help_counter+1
                      END IF
                   END DO
                   
                   cost = 0.0_c_double         ! The initialization of the previous solution
                   kit_num_ed = 0              ! The initialization of the number of kits in the previous solution
                   kits_beta_ed = 0            ! The initialization of kits in the previous solution                    

                   !Calculation of kits in the previous solution
                   DO j1 = 1, nk    ! Each kit is looked through
                     kit_in_use = .FALSE.
                     DO j2 = 1, nft
                       IF (info%mk(j2,j1)==1) THEN 
                         IF ( ABS(x_ed(j2)) >= tol_zero ) THEN 
                           kit_in_use = .TRUE.           
                         END IF
                       END IF
                     END DO
                     IF (kit_in_use) THEN              ! Executed if kit j1 is in the previous solution
                       cost = cost + info%mC(j1)            ! The cost of kit j1 is taken into account
                       kit_num_ed = kit_num_ed + 1     ! The number of kits is updated
                       kits_beta_ed(kit_num_ed) = j1   ! The index of kit j1 is updated to table kits_beta                 
                     END IF 
                   END DO       

                   ! The storing of kits to "betakits" vector
                   ind = (k-1)*nk
                   DO i = 1, nk
                      betakits(ind+i) = kits_beta_ed(i)                
                   END DO 
 
                   IF ((iprint>=2) .OR. (iprint == -1)) THEN 
                     WRITE(*,*)
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*) 'Information about the best solution:'
                     WRITE(*,*)                  
                     WRITE(*,*)  '  f=',small 
                     WRITE(*,*)  '  number of zero elements=',help_counter
                     WRITE(*,*)  '  number of nonzero elements=',nft-help_counter 
                     WRITE(*,*)  '  cost=',cost
                     WRITE(*,*)  '  number of kits=',kit_num_ed 
                     WRITE(*,*)  '  kits=',kits_beta_ed(1:kit_num_ed)    
                     WRITE(*,*)                  
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*)
                     
                   END IF
 
                   IF (iprint==1) THEN
                     WRITE(*,*) small, help_counter, nft-help_counter, &
                     & cost, kit_num_ed, kits_beta_ed(1:kit_num_ed)                
                   END IF
                   
                 END DO

              END IF              
              
              !--------------------------------------------------------------------------
              !                 SOLVING OF L0-NORM PROBLEM ENDS 
              !---------------------------------------------------------------------------  

            
              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk
              !---------------------------------------------------------------------------                



              IF (scale_in_use) THEN     ! Rescaling
                 
                 CALL set_k(info, nk)               
                 info%user_lambda = 0.0_c_double
                 CALL rescaling_mse(info)        ! the rescaling of data 
                  
                 DO k = 1, nk                    ! Each solution is rescaled
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   CALL rescaling_beta_mse(info, beta_solution)        ! Rescaling of solution
                   f1_current = f1(info,beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem1,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current                    ! The objective function value for problem 4 with k nonzero kits 
                   WRITE(*,*) 'f',k,':', fperk(k)
                   
                   l2 = 1   
                   ind = (k-1)*(ncol+1)
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO   
                   beta(ind+ncol+1) = beta_solution(l2)                     
                   
                 END DO 
              ELSE ! No rescaling
              
                CALL set_k(info,nk)               
                 info%user_lambda = 0.0_c_double
                 
                 DO k = 1, nk-1         
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   f1_current = f1(info,beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem2,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current                    ! The objective function value for problem 4 with k nonzero kits         
                    WRITE(*,*) 'f',k,':', fperk(k)
                    
                   l2 = 1   
                   ind = (k-1)*(ncol+1)
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+1)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO   
                   beta(ind+ncol+1) = beta_solution(l2)    

                 END DO   
                 
                   ind = (nk-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = x_koe(i)
                   END DO
                   f1_current = f1(info,beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem2,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current                    ! The objective function value for problem 4 with k nonzero kits         
                    WRITE(*,*) 'f',k,':', fperk(k)

                   l2 = 1   
                   ind = (k-1)*(ncol+1)
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO   
                   beta(ind+ncol+1) = beta_solution(l2) 
                 
            
              END IF       

   
              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk COMPLETED
              !---------------------------------------------------------------------------  
               
             CALL cpu_time(f_time)   ! Finish CPU timing    
             cpu = f_time-s_time             
             
             CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
             ! Calculate the elapsed 'clock' time in seconds:
             elapsed_time=(1.0_c_double*clock_end-clock_start)/clock_rate                
             
             
             IF ((iprint >= 1) .OR. (iprint == -1)) THEN               
                 WRITE(*,*) 'Used CPU:', cpu   , 'Elapsed time:', elapsed_time             
             END IF
       
             CALL deallocate_data_mse(info) 
             CALL deallocate_x()     

         END SUBROUTINE oscar_mse


        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |            THE DOUBLE BUNDLE ALGORITHM FOR LOGISTIC REGRESSION MODEL           |  |
        !  |                                                                                |  |
        !  |                     SOLUTION FOR EVERY NUMBER OF KITS                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        ! Subroutine for oscar with the logistic regression model
        
    !--------------------------------------------------------------------------
    ! ^^^^ START: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------        
          SUBROUTINE oscar_logistic(x, y, kits, costs, nrow, ncol, nkits, beta, fperk, &
       & in_print, in_start, in_k_max, &
       & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
       & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, & 
       & nKitOnes, betakits, &
       & solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
       & in_tolg, in_tolg2, in_eta, in_epsL ) &     
        & BIND(C, name = "oscar_logistic_f_")              
    !--------------------------------------------------------------------------
    ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
    !--------------------------------------------------------------------------


          !--------------------------------------------------------------------------            
          ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------    
          ! SUBROUTINE oscar_logistic(infileX, infileY, infileK, infileC,  &
          !                     & nrow, ncol, nkits, &
          !                     & beta, fperk, in_print, in_start, in_k_max, &
          !                     & in_mrounds, in_mit, in_mrounds_esc, in_b1, in_b2, in_b, &
          !                     & in_m, in_m_clarke, in_c, in_r_dec, in_r_inc, in_eps1, in_eps, in_crit_tol, &                                
          !                     & nKitOnes, betakits, solver_id, & 
          !                     & in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, &
          !                     & in_tolg, in_tolg2, in_eta, in_epsL)    
          !--------------------------------------------------------------------------            
          ! ^^^^ END: If Fortran code is used without R-C-interface ^^^^
          !--------------------------------------------------------------------------   

            !_____________________________________________________________________________________
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            !
            !  !! NOTICE !!  The DATA is SCALED during this SUBROUTINE !
            !
            ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            !
            !
            ! * 'problem = 5'        : the used objective function is the logistic regression model with L0-norm for kits
            !                          (parametric approach since the number of nonzero elements is fixed) 
            !
            !
            !   Solves in a loop all L0-norm problems with the fixed number of k=1,...,in_k_max of kits    (if in_start > 0)
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'ncol'        : The dimension of the problem = the number of features in a predictor, INTEGER
            !         * 'nrow'        : The number of records (data points), INTEGER
            !         * 'nkits'       : the (maximum) number of kits, INTEGER
            !
            !         * 'in_print'    : specifies the print, INTEGER
            !         * 'in_start'    : specifies how starting points are selected when the L0-norm problem is solved for  
            !                           the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'in_k_max'    : specifies how many kits are used in the last problem of the loop  (NEEDS to be 1 <= 'in_k_max' <= nkits !)
            !
            !         * 'nKitOnes'    : the number of value 1 in kit matrix
            !
            !         * 'solver_id'   : the used solver (1=DBDC and 2=LMBM) INTEGER
			!
            !         PARAMETERS in DBDC method:
            !
            !         * in_mrounds     : the maximum number of rounds in one main iteratioin
            !         * in_mit         : the maximum number of main iterations
            !         * in_mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance          
            !
            !         PARAMETERS in LMBM method:
            !
            !         * in_na        : Size of the bundle na >= 2.
            !         * in_mcu       : Upper limit for maximum number of stored corrections, mcu >= 3.
            !         * in_mcinit    : Initial maximum number of stored corrections, mcu >= mcinit >= 3.
            !                          If mcinit <= 0, the default value mcinit = 3 will be used. 
            !                          However, the value mcinit = 7 is recommented.             
            !         * in_tolf      : Tolerance for change of function values (default = 1.0E-5).
            !         * in_tolf2     : Second tolerance for change of function values.
            !                          - If tolf2 < 0 the the parameter and the corresponding termination 
            !                            criterion will be ignored. 
            !                          - If tolf2 = 0 the default value 1.0E+4 will be used
            !         * in_tolg      : Tolerance for the first termination criterion (default = 1.0E-5).   
            !         * in_tolg2     : Tolerance for the second termination criterion (default = tolg). clustering code small data  
            !         * in_eta       : Distance measure parameter, eta >= 0.   
            !         * in_epsL      : Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.)    			
            !
            !         NOTICE: DATA IS GIVEN IN VECTOR FORMAT ! Due to this values for observations/kits are given in these vector one after another
            !
            !         * 'x'          : the matrix X of input variables in vector form --> will be stored to in_mX (each row is one observation)
            !                               REAL, DIMENSION(nrow*ncol)  
            !         * 'y'          : the vector Y of outputs (values are either 0 or 1, i.e. binary) --> will be stored to in_mY 
            !                               INTEGER, DIMENSION(nrow) 
            !         * 'costs'      : the cost vector C for kits --> will be stored to in_mC
            !                               REAL, DIMENSION(knits)
            !         * 'kits'       : the matrix K of kit structures in vector form (values are either 0 or 1, i.e. binary) --> will be stored to in_mK (each row is one kit)
            !                               INTEGER, DIMENSION(knits*ncol)
            !
            !
            ! OUTPUT: * 'fperk'      : The vector containing objective function values at for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION(nkits) 
            !         * 'beta'       : The vector containing solution vectors beta obtained for each L0-norm problem with the fixed number of k=1,...,nkits of kits, 
            !                              REAL, DIMENSION((ncol+1)*nkits)    (in eah solution first 'ncol' values are beta_1,..,beta_nft and the last one is beta_0)
            !
            ! NOTICE: * The dimension of vectors 'fperk' has to be 'nkits'         ('nkits' is the maximum number of kits)
            !         * The dimension of vectors 'beta' has to be '(ncol+1)*nkits' ('ncol+1' is the dimension of the problem and 'nkits' is the number of kits)
            !
            ! OPTIONAL INPUT (CAN BE INCLUDED AS INPUT PARAMETERS IF NEEDED. AT THE MOMENT DEFAULT VALUES ARE USED FOR THEM.):
            !         * 'mit'            : The maximum number of 'main iterations', INTEGER                  
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration', INTEGER
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Escape procedure', INTEGER
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm, LOGICAL           
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm, LOGICAL
            !         * 'scale_in_use'   : If .TRUE. then the data is scaled
            !         * 'CPUtime'       : the CPU time (in seconds) REAL
            !
            ! NOTICE: * 'in_print' has to be 0, 1, 2 or 3. If it is NOT then DEFAULT value 1 is used.             
            !         * 'in_start' has to be  1, 2, 3 or 4. If it is NOT then DEFAULT value 2 is used.             
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! INPUTs
               ! INTEGER(KIND = c_int), INTENT(IN)     :: nrow       ! Number of rows in x (i.e. records)
               ! INTEGER(KIND = c_int), INTENT(IN)     :: ncol       ! Number of cols in x (i.e. features)
               ! INTEGER(KIND = c_int), INTENT(IN)     :: nkits      ! Number of kits for features
               
               ! INTEGER(KIND = c_int), INTENT(IN)    :: in_start    ! Starting point procedure used
               ! INTEGER(KIND = c_int), INTENT(IN)    :: in_print    ! Print used
               ! INTEGER(KIND = c_int), INTENT(IN)    :: in_k_max    ! The maximum number of kits in the loop
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds              ! the maximum number of rounds in one main iteration
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mit                  ! the maximum number of main iterations
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               ! INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               ! REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               ! REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               ! REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               ! REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance    

               ! INTEGER(KIND=c_int), INTENT(IN) :: nKitOnes                ! the number of value 1 in kit matrix
			   ! INTEGER(KIND=c_int), INTENT(IN) :: solver_id               ! the used solver (1=DBDC and 2=LMBM)

!**
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_na                  ! Size of the bundle na >= 2.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcu                 ! Upper limit for maximum number of stored corrections, mcu >= 3.
               ! INTEGER(KIND=c_int), INTENT(IN) ::  in_mcinit              ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! ! However, the value mcinit = 7 is recommented.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf                ! Tolerance for change of function values (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolf2               ! Second tolerance for change of function values.
                                                                          ! !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          ! !   criterion will be ignored. 
                                                                          ! !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg                ! Tolerance for the first termination criterion (default = 1.0E-5).
               ! REAL(KIND=c_double), INTENT(IN) ::  in_tolg2               ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               ! REAL(KIND=c_double), INTENT(IN) ::  in_eta                 ! Distance measure parameter, eta >= 0.
               ! REAL(KIND=c_double), INTENT(IN) ::  in_epsL                ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**		

              ! INPUTs
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: nrow                  ! Number of rows in x (i.e. records)
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: ncol                  ! Number of cols in x (i.e. features)
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: nkits                 ! Number of kits for features
               
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_start              ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_print              ! Print used
               INTEGER(KIND = c_int), INTENT(IN), VALUE :: in_k_max              ! The maximum number of kits in the loop
               
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mrounds              ! the maximum number of rounds in one main iteration
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mit                  ! the maximum number of main iterations
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_mrounds_esc          ! the maximum number of rounds in one escape procedure
               
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN), VALUE :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN), VALUE :: in_crit_tol             ! the stopping tolerance: criticality tolerance

               INTEGER(KIND=c_int), INTENT(IN), VALUE :: nKitOnes                ! the number of value 1 in kit matrix

			   INTEGER(KIND=c_int), INTENT(IN), VALUE :: solver_id               ! the used solver (1=DBDC and 2=LMBM)
 
 !**
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_na           ! Size of the bundle na >= 2.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcu          ! Upper limit for maximum number of stored corrections, mcu >= 3.
               INTEGER(KIND=c_int), INTENT(IN), VALUE ::  in_mcinit       ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                                                          ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                                                          ! However, the value mcinit = 7 is recommented.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf         ! Tolerance for change of function values (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolf2        ! Second tolerance for change of function values.
                                                                          !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                                                          !   criterion will be ignored. 
                                                                          !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg         ! Tolerance for the first termination criterion (default = 1.0E-5).
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_tolg2        ! Tolerance for the second termination criterion (default = tolg). clustering code small data
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_eta          ! Distance measure parameter, eta >= 0.
               REAL(KIND=c_double), INTENT(IN), VALUE ::  in_epsL         ! Line search parameter, 0 < epsL < 0.25 (default = 1.0E-4.) 
!**
 
            !--------------------------------------------------------------------------
            ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nrow*ncol)  :: x        ! Vector of data values
               INTEGER(KIND = c_int), INTENT(IN), DIMENSION(nrow)       :: y        ! Vector of outputs (values are either 0 or 1, i.e. binary)
               INTEGER(KIND = c_int), INTENT(IN), DIMENSION(nkits*ncol) :: kits     ! Vector of kit indicator values (binary indicators)
               REAL(KIND = c_double), INTENT(IN), DIMENSION(nkits)      :: costs    ! Costs associated to each kit
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used with R-C-interface ^^^^
            !--------------------------------------------------------------------------


            !--------------------------------------------------------------------------            
            ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
            !--------------------------------------------------------------------------
            !   CHARACTER(LEN=80), INTENT(IN) :: infileX         ! The name of "predictor matrix" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileY         ! The name of "observed time and label" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileC         ! The name of "kit costs" file 
            !   CHARACTER(LEN=80), INTENT(IN) :: infileK         ! The name of "predictor matrix" file   
            !--------------------------------------------------------------------------
            ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^          
            !--------------------------------------------------------------------------
            
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION((ncol+1)*nkits)  :: beta       !Output variable for beta coefficients per k
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nkits)           :: fperk      !Output variable target function value per k     
               INTEGER(KIND = c_int), INTENT(OUT), DIMENSION(nkits*nkits)     :: betakits   !Output variable telling kits in beta coefficients per k

           
           !***************************** LOCAL VARIABLES ************************************      

               TYPE(set_info) :: info           ! The set of information                     
 
               REAL(KIND=c_double) :: CPUtime                 ! the CPU time (in seconds)

               INTEGER(KIND=c_int) :: solver                  ! defines the solver used (1=DBDC, 2=LMBM)			   

               INTEGER(KIND=c_int) :: nft                     ! the dimension of the problem = the number of features in a predictor
               INTEGER(KIND=c_int) :: nrecord                 ! the number of records (data points)
               INTEGER(KIND=c_int) :: nk                      ! the number of kits 
               INTEGER(KIND=c_int) :: nk_max                  ! the maximum number of kits in the loop
   
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) :: beta_solution    ! the solution vector beta obtained for the problem

               REAL(KIND=c_double), DIMENSION(nKitOnes+1,nkits) :: points     ! the beta_solutions for problem 3 for fixed k ('nkits' different starting points) 
               REAL(KIND=c_double), DIMENSION(nkits) ::      f_points     ! the objective function values for problem 3 for fixed k ('nkits' different starting points)
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) ::       x_koe      ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nKitOnes+1) ::       x_ed       ! the beta solution for the previous problem where the number of nonzero elements was one smaller

               REAL(KIND=c_double), DIMENSION(nKitOnes+1,nkits) ::  x_solution      ! the beta solutions for the problem with fixed starting points and number of kits (for all starting points)
               REAL(KIND=c_double), DIMENSION(nkits) :: f_solution              ! the objective function values at the solution 'x_solution' (for all starting points)
               REAL(KIND=c_double) :: f_solution_DBDC                           ! the f_solution obtained from DBDC method
               
               REAL(KIND=c_double), DIMENSION(nrow,ncol) :: in_mX     ! predictor matrix (row is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrow) :: in_mY          ! output values 
               INTEGER(KIND=c_int), DIMENSION(nkits,ncol) :: in_mK    ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nkits) :: in_mC         ! kit costs                   

               REAL(KIND=c_double), DIMENSION(nKitOnes,nrow) :: mXt        ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(nrow) :: mYt                 ! Outputs  
               INTEGER(KIND=c_int), DIMENSION(nkits,nKitOnes) :: mK        ! Modified kit matrix (There is only one value one in each column)  
               REAL(KIND = c_double), DIMENSION((nKitOnes+1)*nkits)  :: beta_nft   !Output variable for beta coefficients per k

               INTEGER(KIND=c_int), DIMENSION(ncol) :: KitOnes        ! The number of kits where each variables is located. If kits do not intersect then each variable is only in one kit!
                
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta     ! indices of kits in the solution 'beta_solution'
               INTEGER(KIND=c_int), DIMENSION(nkits) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'

               REAL(KIND=c_double), DIMENSION(nKitOnes+1) :: x_0          ! the starting point
                             
               REAL(KIND=c_double), DIMENSION(8) :: mrho        ! Vector containing the values of rho parameter used in the method 
              
               INTEGER(KIND=c_int) :: nstart              ! the number of starting point
               INTEGER(KIND=c_int) :: start_max           ! the number of starting point when 'start = 4'

               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                       
               INTEGER(KIND=c_int) :: mit                ! the maximum number of 'main iterations'
                                                         ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds            ! the maximum number of rounds during one 'main iteration'
                                                         ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int) :: mrounds_clarke     ! the maximum number of rounds during one 'Clarke stationary' algorithm
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

               ! Possible USER PARAMETER
               INTEGER(KIND=c_int) :: iprint   ! specifies the print
                                               !   iprint = 0 : print is suppressed
                                               !   iprint = -1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Not table form)
                                               !   iprint = 1 : prints for each L0-norm problem the final value of Cox's proportional hazard model (Table form)
                                               !   iprint = 2 : prints for each L0-norm problem the final value of Cox's proportional hazard model obtained from each starting point
                                               !   iprint = 3 : prints for each starting point and L0-norm problem all intermediate results together the final results
 
               INTEGER(KIND=c_int) :: start    !   start = 1  : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = 1, 2, ... , nkits
                                               !                (starting point is the solution obtained for L0-norm problem with k-1 nonzero components)
                                               !
                                               !   start = -1 : only one starting point for L0-norm problem with k nonzero components 
                                               !                L0-problems are solved in order: k = nkits, nkits-1, ... , 1
                                               !                (starting point is the solution obtained for L0-norm problem with k+1 nonzero components)  
                                               !
                                               !   start = 2  : for L0-norm problem with k nonzero components uses 'nkits-k+1' starting points
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
               REAL(KIND=c_double) :: small           ! The cost of solution beta
               
               REAL(KIND=c_double) :: s_time, f_time  ! The start and finish times
               REAL(KIND=c_double) :: cpu             ! The cpu time
               
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                     ! Random number
               REAL(KIND=c_double), DIMENSION(nkits) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nkits) :: mRandInd     ! Original indices of random numbers in matrix
               
               INTEGER(KIND=c_int) :: problem1              ! The DC component f1 
               INTEGER(KIND=c_int) :: problem2              ! The DC component f2 
        
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: nremoved
               INTEGER(KIND=c_int) :: kit_num, kit_num_ed   ! The number of kits in the current and previous solutions
               INTEGER(KIND=c_int) :: i, j, k, ind, min_ind, j1, j2, ii, i2, iii
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization

               INTEGER(KIND=c_int), DIMENSION(nkits) :: mPrNum   ! The number of problems for threads
               INTEGER(KIND=c_int) :: tread_num                  ! The number of thread
               INTEGER(KIND=c_int) :: jaannos                    ! The residue in division
               INTEGER(KIND=c_int) :: kokoosa                    ! The integer obtained in division
               INTEGER(KIND=c_int) :: startind, j3               ! Help variables
               INTEGER(KIND=c_int) :: l, l2                      ! Help variables

               REAL(KIND=c_double) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER(KIND=c_int) :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time 
                
!**			   
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.               
 !**
 
               CALL cpu_time(s_time)   ! Start CPU timing
               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time     

!**
              ! The solver is defined/selected
               IF (solver_id == 1) THEN
                   solver = 1 !DBDC
               ELSE IF (solver_id == 2) THEN
                   solver = 2 !LMBM
               ELSE 
                   solver = 1 !DBDC              
               END IF
!**  
               
             ! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
               
             ! Maximum number of possible treads in parallellization
               max_threads = omp_get_max_threads()                
              !max_threads = 1
			   
!**
			   IF (solver == 2) THEN  !No parallellization with LMBM can be used
                   max_threads = 1
			   END If	   
!**			   
               
               CALL omp_set_num_threads(max_threads)
               tread_num = max_threads

               WRITE(*,*) 'The number of threads', max_threads             
               ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----   

               ! The initialization of parametrs used in DBDC methods
               CALL allocate_parameters(info, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)
               
!**                                          
                 ! The initialization of parameters used in LMBM method
                 CALL init_par(in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, &
                                & in_tolg2, in_eta, in_epsL)
!**
			   
               ! Set the number of rows and columns inside Fortran  + kits           
               nrecord = nrow
               nft = nKitOnes
               nk = nkits
               
               ! The maximum number of kits in the loop
               nk_max = min(nk,in_k_max)         ! Cannot be greater than nk
               nk_max = max(1,nk_max)            ! Cannot be smaller than 1
              
               start = in_start       ! Starting point generation procedure
               iprint = in_print      ! Print option
              
               ! The default print is used if user specifieed value is not acceptable
               IF (iprint < -1 .OR. iprint > 4) THEN
                  iprint = 1_c_int
               END IF   

               ! The default start is used if user specifieed value is not acceptable
               IF ((start> 4) .OR. (start <= 0)) THEN
                  start = 2_c_int
               END IF              
               
               ! If start = 5 then start_max is the number of starting points in each L0-norm problem
               start_max = 5_c_int
               
               mrounds = in_mrounds            ! maximum number of rounds during one 'main iterations'
               mit = in_mit                    ! maximum number of 'main iteration'
               mrounds_clarke = in_mrounds_esc ! maximum number of rounds during one 'Clarke stationary' algorithm
                    
               agg_used = .TRUE.          ! Aggregation is used
               stepsize_used = .FALSE.    ! Simple stepsize determination is not used
          
               scale_in_use = .TRUE.      ! The scaling of data is used
            
               ed_sol_in_pen = .FALSE.
               !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem

               !Values for parameter rho used in penalization problem 
               !mrho = (/0.1_c_double, 0.2_c_double, 0.5_c_double, 1.0_c_double &
               !      & 2.0_c_double, 5.0_c_double, 10.0_c_double, 20.0_c_double /) 
               mrho = (/0.2_c_double, 0.5_c_double, 1.0_c_double, 2.0_c_double, &
                     & 5.0_c_double, 10.0_c_double, 20.0_c_double, 50.0_c_double /) 

                ! Problem 
                problem1 = 5_c_int
                problem2 = 5_c_int                
               
                user_n = nft+1
                
                tol_zero = (10.0_c_double)**(-6)
       
               IF ( (iprint > 0) .OR. (iprint == -1) ) THEN 
                 IF (ed_sol_in_pen) THEN 
                  WRITE(*,*) 'When the penalized problems is solved we utilize&
                              & the previous solution as a starting point.'
                 ELSE
                  WRITE(*,*) 'When the penalized problems is solved we do NOT utilize&
                              & the previous solution as a starting point.'
                 END IF
                
                 IF ( start > 0 ) THEN 
                  WRITE(*,*) 'Problems are solved from smallest to largest.'
                 ELSE   
                  WRITE(*,*) 'Problems are solved from largest to smallest.'
                 END IF
               
                 WRITE(*,*) 'Starting points are generated with the procedure', start 
                 WRITE(*,*) 'The rho values in the penalized problem:', mrho
               
               END IF                  
               
              !--------------------------------------------------------------------------------------
          
              ! The starting point
        
                  x_0 = 0.0_c_double
                 
              !---------------------------------------------------------------------------
              !                       POPULATING DATA MATRICES
              !---------------------------------------------------------------------------

              !---------------------------------------------------------------------------             
              ! ^^^^ START: If Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       READING DATA MATRICES
              !---------------------------------------------------------------------------
          
                ! OPEN(78,file=infileX,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) (in_mX(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileY,status='old',form='formatted')
                ! DO i=1,nrecord
                   ! READ(78,*) in_mY(i)
                ! END DO
                ! CLOSE(78)        

                ! OPEN(78,file=infileK,status='old',form='formatted')      
                ! DO i=1,nk
                   ! READ(78,*) (in_mK(i,j),j=1,ncol)
                ! END DO
                ! CLOSE(78)
            
                ! OPEN(78,file=infileC,status='old',form='formatted')
                ! DO i=1,nk
                   ! READ(78,*) in_mC(i)
                ! END DO
                ! CLOSE(78)  
                
              !--------------------------------------------------------------------------
              ! ^^^^ END: IF Fortran code is used without R-C-interface ^^^^
              !---------------------------------------------------------------------------

 
              !---------------------------------------------------------------------------            
              ! ^^^^ START: If Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------
              !                       POPULATE MATRICES
              !--------------------------------------------------------------------------- 
              
               ! Populate the input matrix 'in_mX' of dim {nrecord,nft}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nrecord
                      ind = ind +1
                      in_mX(i,j) = x(ind)
                    END DO
                 END DO
                 
               ! Populate the output vector 'in_mY' of dim {nrecord}
                 ind = 0
                 DO i = 1, nrecord
                    ind = ind + 1   
                    in_mY(i) = y(ind)
                 END DO
                 
              ! Populate the kit matrix 'in_mK' of dim {nkits,nft}
                 ind = 0
                 DO j = 1, ncol
                    DO i = 1, nk
                      ind = ind + 1
                      in_mK(i,j) = kits(ind)
                    END DO
                 END DO

              ! Populate the cost vector 'in_mC' for kits of dim {nkits}
                 ind = 0
                 DO i = 1, nk
                    ind = ind + 1
                    in_mC(i) = costs(ind)
                 END DO 
                 
              !--------------------------------------------------------------------------
              !^^^^ END: IF Fortran code is used with R-C-interface ^^^^
              !---------------------------------------------------------------------------
               
           
              !---------------------------------------------------------------------------
              !               ALLOCATION OF DATA MATRICES 
              !---------------------------------------------------------------------------
                   
               ! Allocation of data matrices in function.f95
               CALL allocate_data_log(info,nft,nrecord,nk,nk)   

                ! The number of value ones in columns of kit matrix
                KitOnes = 0_c_int
                DO j = 1, ncol
                   DO i = 1, nk
                       KitOnes(j) = KitOnes(j) + in_mK(i,j)
                   END DO
                END DO
               
              !---------------------------------------------------------------------------
              !                     STORING DATA MATRICES 
              !---------------------------------------------------------------------------
              ! Notice: Matrices are transposed! => Each column presents either an observation or a kit!
             
              !Transpose of matrix mX
               DO i = 1, nrecord
                 l2 = 1
                 DO j = 1, ncol
                   DO l = 1, KitOnes(j) 
                       mXt(l2,i) = in_mX(i,j)
                       l2 = l2 + 1
                   END DO   
                 END DO
               END DO   
            
               !Transpose of matrix mY          
               DO i = 1, nrecord
                    mYt(i) = in_mY(i)
               END DO

               ! Alters the kit matrix
               mK = 0_c_int
               l2 = 1_c_int
               DO j = 1, ncol 
                  DO i = 1, nk
                     IF (in_mK(i,j)==1_c_int) THEN  
                        mK(i,l2) = in_mK(i,j)
                        l2 = l2 + 1_c_int
                     END IF
                  END DO               
               END DO              

               CALL allocate_matrices_log(info, mXt, mYt, mK, in_mC,  & 
                                      & nrecord, nft, nk)

              ! Scaling of data 
              IF (scale_in_use) THEN
                  CALL scaling_log(info)
                  
                 ! Original data matrix mXt is changed to the scaled one 
                  CALL scaling_matrix_log(info, mXt)
                  
              END IF  


              ! The initialization of beta vector
              beta = 0.0_c_double
              beta_nft = 0.0_c_double
              
              
              ! The best beta_solution for Cox's proportional hazard model without regularization/penalization    
              
              CALL set_k(info, nkits)           ! All kits can be used

               iprint_DBDC = 0_c_int           ! basic print of intermediate results and extended print of final results
                            
!**                      
			  IF (solver == 1) THEN	  
                  CALL DBDC_algorithm( f_solution_DBDC, x_koe, x_0, 0.0_c_double, 0.0_c_double, &
                            & mit, mrounds, mrounds_clarke, termination, counter, CPUtime,  &
                            & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n, &
                            & info)            
			  ELSE IF (solver == 2) THEN 
			  
                  CALL allocate_xn(user_n) 
                  CALL init_LMBMinfo(problem1, info) 
                  CALL init_x(x_0)          
                  CALL set_rho_LMBM(0.0_c_double)                 
                  CALL set_lambda_LMBM(0.0_c_double)                  
                  CALL lmbm(mc,f_solution_DBDC,iout(1),iout(2),iout(3),iout(4),LMBMstart)     
                  CALL copy_x(x_koe)
                  			  
              END IF			  
!**           
 
              ! Notice: * solution x_koe is obtained by fitting Cox's model to data without regularization
              !         * x_koe is utilized in formation of starting points   
              
               x_ed = 0.01_c_double    ! We initialize the previous solution, since we do not have a value for it 
               
               
               IF ((iprint > 0) .OR. (iprint == -1)) THEN
                 WRITE(*,*) 
                 WRITE(*,*) 'Value of logistic model without regularization:', f_solution_DBDC
                 WRITE(*,*) 
               END IF 
               
               IF (iprint==1) THEN
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
                 WRITE(*,*)  '  f  ', '  zero_elements  ', '  nonzero_elements  ', '  cost  ', '  num_kits  ', 'kits'   ! Prints the solution to the file 'ratkaisu.txt'
                 WRITE(*,*) '------------------------------------------------------------------------------------------'
               END IF              
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
               
                 kit_num_ed = 0           ! The number of kits in the previous solution
                 kits_beta_ed = 0         ! The kits in the previous solution
                 small = (10_c_double)**9 ! There is no previous solution known
                 
                 IF (start == 1) THEN
                    nstart = 1
                 ELSE IF (start == 2 .OR. start == 3) THEN
                    nstart = nk
                 ELSE IF (start == 4) THEN 
                    nstart = start_max        ! The number of random starting points
                 END IF 
                 
                 !---------------------------------------------------------------------------------------------
                 ! Solves in a loop all L0-norm problems with the fixed number of k=1,..., nk_max of kits    
                 !---------------------------------------------------------------------------------------------
                 DO k = 1, nk_max              ! In this loop we calculate the solution for problem 3 wiht L0-norm such that the number of kits varies from 1 to nk_max
                  
                  IF (iprint>=2 .OR. iprint == -1) THEN
                       WRITE(*,*) '-------------------------------------' 
                       WRITE(*,*) 'PROBLEM with', k, 'kits' 
                       
                  END IF
                  
                  CALL set_k(info, k)            ! The number of nonzero kits is fixed
                  
                  jaannos = MOD(nk,tread_num)    ! The residue in the division
                  kokoosa = nk/tread_num         ! The integer part obtained from the division
                  
                  mPrNum = kokoosa               ! The initialization of the number of problems for threads
                  DO i = 1, jaannos              ! We increase the number of problems with one for specific threads
                     mPrNum(i) = mPrNum(i)+1                  
                  END DO

                  !----------------------------------------------------------
                  ! Starting points for problem with fixed number of k kits 
                  !----------------------------------------------------------         
 
                  !$OMP PARALLEL DO PRIVATE(x_solution, f_solution, i) & 
                  !$OMP FIRSTPRIVATE(x_ed, x_koe, kits_beta_ed, kit_num_ed)  &             
                  !$OMP FIRSTPRIVATE(small, k, iprint, mrho)                 &  
                  !$OMP FIRSTPRIVATE(mit, mrounds, mrounds_clarke)           &    
                  !$OMP FIRSTPRIVATE(agg_used, stepsize_used)                &
                  !$OMP FIRSTPRIVATE(problem1, problem2)                     &  
                  !$OMP FIRSTPRIVATE(in_b1, in_b2, in_m, in_c, in_r_dec)     &  
                  !$OMP FIRSTPRIVATE(in_r_inc, in_eps1, in_b, in_m_clarke)   &  
                  !$OMP FIRSTPRIVATE(in_eps, in_crit_tol)                    &  
                  !$OMP FIRSTPRIVATE(nft, nrecord, nk, user_n)               &  
                  !$OMP FIRSTPRIVATE(mXt, mYt, mK, in_mC)                    &  
                  !$OMP SHARED(points, f_points, mPrNum)               
                  
                   DO i = 1, tread_num          ! Different starting points for threads are looked through to solve the problem 3 with fixed number of nonzero kits                   
                    
                      CALL solution_with_k_kits_log(x_solution, f_solution, k, i, small, &
                                               & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                               & nk, start, iprint, mrho, mit, mrounds, mrounds_clarke, &
                                               & agg_used, stepsize_used, nft, problem1, problem2, &
                                               & mXt, mYt, mK, in_mC, nrecord,  & 
                                               & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                               & in_b, in_m_clarke, in_eps, in_crit_tol, mPrNum, solver)
                    
                    ! Storing of the obtained solution and the corresponding objective value
                    !$OMP CRITICAL 
                     startind = 1
                     IF (i > 1) THEN    
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                     END IF                        
                        
                     DO j3 = 1, mPrNum(i)   
                        DO j = 1, user_n
                            points(j,startind+(j3-1)) = x_solution(j,startind+(j3-1))
                        END DO     
                      f_points(startind+j3-1) = f_solution(startind+j3-1)
                     END DO  
                    !$OMP END CRITICAL                    
                        
                   END DO          
                  !$OMP END PARALLEL DO   
                   
                   IF (start == 1) THEN 
                      small = f_points(1)
                      min_ind = 1
                   END IF 
                   
                   IF (start == 2 .OR. start == 3) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nk           ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF

                   IF (start == 4) THEN 
                     ! The selection of the best solution for problem 3 with k kits
                     small = (10.0_c_double)**9      ! The initialization of the smallest objective function values in table 'f_points'
                     min_ind = 0               ! The index of 'f_points' yielding the smallest value 
                     DO i = 1, nstart          ! The determination of the smallest objective function value
                       IF (small >= f_points(i)) THEN
                         small = f_points(i)
                         min_ind = i
                       END IF
                     END DO
                   END IF                  
                  
                   ind = (k-1)*user_n  
                   DO i = 1, user_n
                     x_ed(i) = points(i,min_ind)            ! The beta solution yielding the smallest objective function value is set as the previous solution
                     beta_nft(ind+i) = points(i,min_ind)        ! The beta solution yielding the smallest objective function value is stored to solution vector 
                   END DO
              
                   IF (iprint >= 1) THEN
                     !WRITE(*,*) 'Objective function value for problem 3 with', k, 'nonzero kits:' , small
                   END IF
              
                   help_counter = 0           
                   DO j = 1, nft         ! The number of zero components in the previous beta vector is computed
                      IF ( ABS(x_ed(j)) < tol_zero ) THEN
                        help_counter = help_counter+1
                      END IF
                   END DO

                   cost = 0.0_c_double         ! The initialization of the previous solution
                   kit_num_ed = 0              ! The initialization of the number of kits in the previous solution
                   kits_beta_ed = 0            ! The initialization of kits in the previous solution                    

                   !Calculation of kits in the previous solution
                   DO j1 = 1, nk    ! Each kit is looked through
                     kit_in_use = .FALSE.
                     DO j2 = 1, nft
                       IF (info%mk(j2,j1)==1) THEN 
                         IF ( ABS(x_ed(j2)) >= tol_zero ) THEN 
                           kit_in_use = .TRUE.           
                         END IF
                       END IF
                     END DO
                     IF (kit_in_use) THEN              ! Executed if kit j1 is in the previous solution
                       cost = cost + info%mC(j1)       ! The cost of kit j1 is taken into account
                       kit_num_ed = kit_num_ed + 1     ! The number of kits is updated
                       kits_beta_ed(kit_num_ed) = j1   ! The index of kit j1 is updated to table kits_beta                 
                     END IF 
                   END DO       

                   ! The storing of kits to "betakits" vector
                   ind = (k-1)*nk
                   DO i = 1, nk
                      betakits(ind+i) = kits_beta_ed(i)                
                   END DO
 
                   IF ((iprint>=2) .OR. (iprint == -1)) THEN 
                     WRITE(*,*)
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*) 'Information about the best solution:'
                     WRITE(*,*)                  
                     WRITE(*,*)  '  f=',small 
                     WRITE(*,*)  '  number of zero elements=',help_counter
                     WRITE(*,*)  '  number of nonzero elements=',nft-help_counter 
                     WRITE(*,*)  '  cost=',cost
                     WRITE(*,*)  '  number of kits=',kit_num_ed 
                     WRITE(*,*)  '  kits=',kits_beta_ed(1:kit_num_ed)    
                     WRITE(*,*)                  
                     WRITE(*,*) '-**--**--**--**--**--**--**--**--**--**--**--**-'
                     WRITE(*,*)
                     
                   END IF
 
                   IF (iprint==1) THEN
                     WRITE(*,*) small, help_counter, nft-help_counter, &
                     & cost, kit_num_ed, kits_beta_ed(1:kit_num_ed)                
                   END IF
                   
                 END DO


              END IF              
              
              !--------------------------------------------------------------------------
              !                 SOLVING OF L0-NORM PROBLEM ENDS 
              !---------------------------------------------------------------------------  

            
              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk
              !---------------------------------------------------------------------------                

              IF (scale_in_use) THEN  ! Resclaing  
                 
                 CALL set_k(info,nk)               
                 info%user_lambda = 0.0_c_double
                 CALL rescaling_mse(info)        ! the rescaling of data 
                 
                 DO k = 1, nk         ! Each solution is rescaled
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   CALL rescaling_beta_mse(info, beta_solution)         ! Rescaling of solution
                   f1_current = f1(info,beta_solution,problem1,user_n)  ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem1,user_n)  ! The f_2 value
                   fperk(k) = f1_current-f2_current                     ! The objective function value for problem 5 with k nonzero kits 
                   WRITE(*,*) 'f',k,':', fperk(k)
                  
                   l2 = 1   
                   ind = (k-1)*(ncol+1)
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO   
                   beta(ind+ncol+1) = beta_solution(l2)      
                   
                 END DO 
              ELSE   ! No rescaling
              
                  CALL set_k(info, nk)               
                 
                  DO k = 1, nk         
                   ind = (k-1)*user_n
                   ! The solution under consideration is stored to 'beta_solution' vector
                   DO i = 1, user_n
                        beta_solution(i) = beta_nft(ind+i)
                   END DO
                   f1_current = f1(info,beta_solution,problem1,user_n) ! The f_1 value 
                   f2_current = f2(info,beta_solution,problem2,user_n) ! The f_2 value
                   fperk(k) = f1_current-f2_current                    ! The objective function value for problem 5 with k nonzero kits         
                    WRITE(*,*) 'f',k,':', fperk(k)

                   l2 = 1   
                   ind = (k-1)*(ncol+1)
                   DO i = 1, ncol 
                      DO j = 1, KitOnes(i)
                        beta(ind+i) = beta(ind+i)+beta_solution(l2)        ! the beta vector for problem 3 with k nonzero kits
                        l2 = l2 + 1 
                      END DO    
                   END DO   
                   beta(ind+ncol+1) = beta_solution(l2)    
                   
                 END DO               
              
              END IF  

              !--------------------------------------------------------------------------
              !                 SOLUTION VECTORS beta and OBJECTIVE VALUES fperk COMPLETED
              !---------------------------------------------------------------------------  
               
             CALL cpu_time(f_time)   ! Finish CPU timing    
             cpu = f_time-s_time             
             
             CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
             ! Calculate the elapsed 'clock' time in seconds:
             elapsed_time=(1.0_c_double*clock_end-clock_start)/clock_rate                
             
             
             IF ((iprint >= 1) .OR. (iprint == -1)) THEN               
                 WRITE(*,*) 'Used CPU:', cpu   , 'Elapsed time:', elapsed_time             
             END IF
       
             CALL deallocate_data_log(info)
	     CALL deallocate_x() 
             
         END SUBROUTINE oscar_logistic           
        !---------------------------------------------------------------------------------------
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |             SOLVING THE PROBLEM WITH k KITS AND STARTING POINT i               |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
         SUBROUTINE solution_with_k_kits_cox( x_solution, f_solution, k, i, f_prev, & 
                                    & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                    & nk, start, iprint, mrho, mit, mrounds, mrounds_esc, &
                                    & agg_used, stepsize_used, nft, problem1, problem2,  &      
                                    & in_mX, in_mY, in_mK, in_mC, nrecord, & 
                                    & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                    & in_b, in_m_clarke, in_eps, in_crit_tol, mPrNum, solver) 
            !_____________________________________________________________________________________
            ! 
            !           
            ! * 'problem1  and problem2   : defines the used objective function is Cox's proportional hazards model with L0-norm for kits
            !                              (parametric approach since the number of nonzero elements is fixed) 
            !
            !   Solves the L0-norm problem from a fixed starting points for thread "i" with "k" kits
            !                                                 
            !
            ! INPUT:  * 'nft'           : The dimension of the problem = the number of features in a predictor, INTEGER
            !         * 'nk'            : the (maximum) number of kits, INTEGER
            !         * 'nrecord'       : the number of observations, INTEGER
            !
            !         * 'k'             : the number of kits in the problem, INTEGER
            !         * 'i'             : the number of thread, INTEGER
            !
            !         * 'x_ed'          : the beta solution for the previous problem where the number of nonzero kits was one smaller (INTEGER)
            !         * 'x_koe'         : The solution to Cox's proportional hazard model without regularization
            !         * 'kit_num_ed'    : The number of kits in the previous solution
            !         * 'kit_beta_ed'   : indices of kits in the previous solution 'x_ed'
            ! 
            !         * 'print'         : specifies the print, INTEGER
            !         * 'start'         : specifies how starting points are selected when the L0-norm problem is solved for  
            !                             the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'mrho'          : List of used rho values
            !
            !         * 'problem1'      : Defines the DC component f1 of the objective
            !         * 'problem2'      : Defines the DC component f2 of the objective
            !
            !         * 'in_mX', 'in_mY', 'in_mK', 'in_mC'    : Are the data matrices 
            !         * 'mPrNum'        : the information about starting points for the thread
            !         
            !         * 'solver'   : defines the solver used (1=DBDC, 2=LMBM)
            !         
            !         PARAMETERS in DBDC method:
            !
            !         * mrounds     : the maximum number of rounds in one main iteratioin
            !         * mit         : the maximum number of main iterations
            !         * mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * agg_used       : .TRUE. if aggregation is used 
            !         * stepsize_used  : .TRUE. if Simple stepsize determination is not used
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance               !
            !
            ! OUTPUT: * 'f_solution'      : The objective function values for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nk) 
            !         * 'x_solution'      : The vector containing solution vector betas for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nft,nk) 
            !
            ! INOUT:  * 'info'            : The set of information
            !
            !
            ! NOTICE: * The dimension of vectors 'x_solution', 'x_ed' and 'x_koe' has to be 'nft' ('nft' is the dimension of the problem)
            !
            !
            ! NOTICE: * 'print' has to be 0, 1, 2 or 3.           
            !         * 'start' has to be 1, 2, 3 or 4.              
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! INPUTs
               REAL(KIND=c_double), DIMENSION(nft), INTENT(IN) :: x_koe    ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nft), INTENT(IN) :: x_ed     ! the beta solution for the previous problem where the number of nonzero kits was one smaller   
               
               REAL(KIND=c_double), INTENT(IN) :: f_prev                   ! The objective value for the previous problem where the number of nonzero kits was one smaller       
               
               INTEGER(KIND=c_int), INTENT(IN) :: kit_num_ed                   ! The number of kits in the previous solution
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: mPrNum        ! the starting points for the thread

               REAL(KIND=c_double), DIMENSION(8), INTENT(IN) :: mrho      ! Vector containing the values of rho parameter used in the method 
            
               INTEGER(KIND = c_int), INTENT(IN) :: nft                   ! Number of features in x (also the dimension of the problem)
               INTEGER(KIND = c_int), INTENT(IN) :: nk                    ! Number of kits for features
               INTEGER(KIND = c_int), INTENT(IN) :: nrecord               ! Number of observations in data 
               
               INTEGER(KIND = c_int), INTENT(IN) :: start                 ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN) :: iprint                ! Print used

               INTEGER(KIND = c_int), INTENT(IN) :: solver                ! the solver used
               
               INTEGER(KIND = c_int), INTENT(IN) :: k                     ! The number of kits in the problem 
               INTEGER(KIND = c_int), INTENT(IN) :: i                     ! The number of the thread
               
               INTEGER(KIND = c_int), INTENT(IN) :: problem1              ! The DC component f1 of the objective               
               INTEGER(KIND = c_int), INTENT(IN) :: problem2              ! The DC component f2 of the objective               
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds                 ! the maximum number of rounds in one main iteration
                                                                          ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mit                     ! the maximum number of main iterations
                                                                          ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
                                                                          
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds_esc             ! the maximum number of rounds in one escape procedure
                                                                          ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
                
               LOGICAL, INTENT(IN) :: agg_used                            ! .TRUE. if aggregation is used in DBDC
               LOGICAL, INTENT(IN) :: stepsize_used                       ! .TRUE. if Simple stepsize determination is not used
               
              ! Parameters in DBDC method
               INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance

              ! The data matrices
               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(2,nrecord) :: in_mY         ! observed times and labels matrix (column is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          
                     
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nft,nk)  :: x_solution   !Output variable for beta coefficients with k kits and starting points for the thread i
               REAL(KIND = c_double), DIMENSION(nk), INTENT(OUT)      :: f_solution   !Output variable target function value with k kits and starting points for the thread i
               
           
           !***************************** LOCAL VARIABLES ************************************      
 
               TYPE(set_info) :: set                                    ! The set of information                     

               REAL(KIND=c_double) :: CPUtime                           ! the CPU time (in seconds)
               REAL(KIND=c_double) :: f_solution_DBDC                   ! the f_solution obtained from DBDC method
  
               REAL(KIND=c_double), DIMENSION(nft) :: beta_solution     ! the solution vector beta obtained during the execution      
          
               INTEGER(KIND=c_int), DIMENSION(nk) :: kits_beta          ! indices of kits in the solution 'beta_solution'

               REAL(KIND=c_double), DIMENSION(nft) :: x_0               ! the starting point
               
               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                                      
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
               
               LOGICAL :: kit_in_use            ! If .TRUE. kit is used in the solution
               LOGICAL :: run_stop              ! If .TRUE. run of this subroutine is stopped for selected k
               LOGICAL :: mukana                ! If .TRUE. specific kit is in the solution
               
               LOGICAL :: ed_sol_in_pen         ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
               LOGICAL :: new_start             ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
       
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
                            
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                  ! Random number
               REAL(KIND=c_double), DIMENSION(nk) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nk) :: mRandInd     ! Original indices of random numbers in matrix
        
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: kit_num               ! The number of kits in the current solution
               INTEGER(KIND=c_int) :: j, j1, j2, ii, i2
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization         

               INTEGER(KIND=c_int) :: startind, j3, kk      ! the start index for starting point    
 
!**
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.
 !** 
                    !--------------------------------------------        
                    !                INITIALIZATION
                    !--------------------------------------------
                    
                    user_n = nft              ! The dimension of the problem
                    iprint_DBDC = 0_c_int     ! The print is supressed 
                  
                    ed_sol_in_pen = .FALSE.
                   !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem  

                    tol_zero = (10.0_c_double)**(-6)
 

                   ! The initialization of parametrs used in DBDC methods
                   CALL allocate_parameters(set, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)
  
                   ! Allocation of sizes od data matrices in function.f95
                    CALL allocate_data_cox(set, nft, nrecord, nk, k)   
                    
                   ! Allocation of data matrices    
                    CALL allocate_matrices_cox(set, in_mX, in_mY, in_mK, in_mC,  & 
                                      & nrecord, nft, nk)
                    
                    CALL  set_k(set, k)                        ! The number of nonzero kits is fixed
                    
                     ! Determination of the first starting point for the thread i                    
                    startind = 1
                    IF (i > 1) THEN 
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                    END IF

                    ! The starting points for the thread are looked through
                    DO kk = startind, startind + mPrNum(i) - 1 
                        !--------------------------------------------
                        !    THE GENERATION OF THE STARTING POINT
                        !--------------------------------------------
                         x_0 = x_ed                  ! The base of the starting point is the previous solution (k-1 nonzero kits)
                         mukana = .FALSE.           
                         IF (start == 2 .OR. start == 3) THEN     ! Only nkits-k+1 starting points used             
                           DO ii = 1, kit_num_ed                  ! Then we check if the kit kk belongs to the previous solution
                             IF (kk == kits_beta_ed(ii)) THEN   
                               mukana = .TRUE.                    ! 'mukana' tells if the kit kk is in solution     
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
                        
                         
                       IF (new_start) THEN       ! .TRUE. if a new starting point is generated
                        
                        IF ((start == 2) .OR. (start == 3)) THEN 
                          DO ii = 1, nft               ! the features are looke
                            IF (set%mK(ii,kk)==1) THEN
                               x_0(ii) = x_koe(ii)    ! In starting poitnt x_0, we initialize features in kit kk with values from solution x_koe
                            END IF   
                          END DO
                        END IF  
                        
                        IF (start == 4) THEN ! A random starting point is generated
                           x_0 = 0.01_c_double
                           DO ii = 1, nk
                             CALL RANDOM_NUMBER(random_num)
                             mRand(ii) = random_num
                             mRandInd(ii) = ii
                           END DO 
                           CALL heapsort_ind(mRand,mRandInd)                    
                           DO ii = 1, k
                              DO i2 = 1, nft 
                                IF (set%mK(i2,mRandInd(ii))==1) THEN
                                  x_0(i2) = x_koe(i2)    ! In starting point x_0, we initialize features in kit i with values from solution x_koe
                                END IF   
                              END DO     
                           END DO
                        END IF                  
                          
                        run_stop = .FALSE.          ! Initialization of 'run_stop' to .FALSE. since we cannot stop
                        num_rho = 0                 ! The selection of the first value of penalization parameter rho
                      
                        IF (iprint >= 2) THEN
                           WRITE(*,*) '-------------------------------------' 
                           WRITE(*,*) 'New start point used.' 
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
                          rho = 10.0_c_double * rho
                          
                          IF (num_rho>=24) THEN 
                          run_stop = .TRUE.         ! Forced stop in optimization
                          END IF 
                        END IF
                    
                        ! The optimization problem is solved fot the current rho and x_0
!**						
						IF (solver==1) THEN 
         
                            CALL DBDC_algorithm( f_solution_DBDC, beta_solution, x_0, rho, 0.0_c_double, &
                                & mit, mrounds, mrounds_esc, termination, counter, CPUtime,  &
                                & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, nft, &
                                & set) 
								
						ELSE IF (solver == 2) THEN
						
                           CALL init_LMBMinfo(problem1, set)   
                           CALL init_x(x_0)
                           CALL set_rho_LMBM(rho)                 
                           CALL set_lambda_LMBM(0.0_c_double)              
                           CALL lmbm(mc, f_solution_DBDC, iout(1),iout(2),iout(3),iout(4),LMBMstart)      
                           CALL copy_x(beta_solution)   
						   
                        END IF
!**						
                
                        IF (ed_sol_in_pen) THEN 
                          x_0 = beta_solution   ! Starting point for the next round
                        END IF    
                    
                        help_counter = 0           
                        DO j = 1, nft         ! The number of zero components in beta vector is computed
                          IF ( ABS(beta_solution(j)) < tol_zero) THEN
                            help_counter = help_counter +1
                            beta_solution(j) = 0.0_c_double
                          END IF
                        END DO
                  
                        cost = 0.0_c_double      ! The initialization of the cost of 'beta_solution'
                        kit_num = 0              ! The initialization of the number of kits in 'beta_solution'
                        kits_beta = 0            ! The initialization of kits in 'beta_solution'
                      
                        !Calculation of kits in solution
                        DO j1 = 1, nk        ! Each kit is looked through
                          kit_in_use = .FALSE.  
                          DO j2 = 1, nft
                            IF (set%mK(j2,j1)==1) THEN 
                              IF ( ABS(beta_solution(j2)) >= tol_zero) THEN    ! If the condition is satisfied the kit j1 is in solution 'beta_solution'
                                 kit_in_use = .TRUE.                                  
                              END IF
                            END IF
                          END DO
                          IF (kit_in_use) THEN       ! Executed if kit j1 is in solution 'beta_solution'
                            cost = cost + set%mC(j1)     ! The cost of kit j1 is taken into account
                            kit_num = kit_num + 1    ! The number of kits is updated
                            kits_beta(kit_num) = j1  ! The index of kit j1 is updated to table kits_beta                      
                          END IF 
                        END DO
        
                        ! We check if the optimization of problem 3 with k kits can be stopped            
                        IF (kit_num <= k .OR. run_stop) THEN 
                           IF (run_stop) THEN 
                             f_solution(kk) = (10.0_c_double)**10
                           END IF
                           run_stop = .TRUE.     ! The optimization can be stopped
                           DO j = 1, nft
                             x_solution(j,kk) = beta_solution(j)        ! We store the solution to x_solution
                           END DO       
                           
                           f1_current = f1(set,beta_solution,problem1,user_n)
                           f2_current = f2(set,beta_solution,problem2,user_n)
                           f_solution(kk) = f1_current-f2_current       ! We store the objective funtion value without penalization to f_solution
                        END IF               
                       
                        f1_current = f1(set,beta_solution,problem1,user_n) 
                        f2_current = f2(set,beta_solution,problem2,user_n)
        
                        IF (iprint > 2) THEN
                           WRITE(*,*) 'rho', rho, 'f',f1_current-f2_current, 'kits', kit_num   
                        END IF
                        
                        IF (run_stop) THEN
                          IF ((iprint >= 2) .AND. (kit_num <= k)) THEN   
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num
                          END IF  
                          IF ((iprint >=2) .AND. (kit_num > k)) THEN                                      
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num, 'should be equal to', k 
                          END IF
                        END IF
                        
                        END DO
                        
                       ELSE
                       
                         DO j = 1, nft
                            x_solution(j,kk) = x_ed(j)
                         END DO     
                         f_solution(kk) = f_prev
                       END IF           
                   
                   END DO
                   
                   CALL deallocate_data_cox(set)             
             
         
         END SUBROUTINE  solution_with_k_kits_cox
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|



       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |             SOLVING THE PROBLEM WITH k KITS AND STARTING POINT i               |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
         SUBROUTINE solution_with_k_kits_mse( x_solution, f_solution, k, i, f_prev, & 
                                    & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                    & nk, start, iprint, mrho, mit, mrounds, mrounds_esc, &
                                    & agg_used, stepsize_used, nft, problem1, problem2,  &      
                                    & in_mX, in_mY, in_mK, in_mC, nrecord, & 
                                    & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                    & in_b, in_m_clarke, in_eps, in_crit_tol ,mPrNum, solver) 
            !_____________________________________________________________________________________
            ! 
            !           
            ! * 'problem1  and problem2   : defines the used objective function is the mean square error model with L0-norm for kits
            !                              (parametric approach since the number of nonzero elements is fixed) 
            !
            !   Solves the L0-norm problem from a fixed starting points for thread "i" with "k" kits
            !                                                 
            !
            ! INPUT:  * 'nft'           : Number of rows. The dimension of the problem = nft + 1, INTEGER
            !         * 'nk'            : the (maximum) number of kits, INTEGER
            !         * 'nrecord'       : the number of observations, INTEGER
            !
            !         * 'k'             : the number of kits in the problem, INTEGER
            !         * 'i'             : the number of thread, INTEGER
            !
            !         * 'x_ed'          : the beta solution for the previous problem where the number of nonzero kits was one smaller (INTEGER)
            !         * 'x_koe'         : The solution to Cox's proportional hazard model without regularization
            !         * 'kit_num_ed'    : The number of kits in the previous solution
            !         * 'kit_beta_ed'   : indices of kits in the previous solution 'x_ed'
            ! 
            !         * 'print'         : specifies the print, INTEGER
            !         * 'start'         : specifies how starting points are selected when the L0-norm problem is solved for  
            !                             the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'mrho'          : List of used rho values
            !
            !         * 'problem1'      : Defines the DC component f1 of the objective
            !         * 'problem2'      : Defines the DC component f2 of the objective
            !
            !         * 'in_mX', 'in_mY', 'in_mK', 'in_mC'    : Are the data matrices 
            !         * 'mPrNum'        : the information about starting points for the thread
            !         
            !         * 'solver'   : defines the solver used (1=DBDC, 2=LMBM)
            !         
            !         PARAMETERS in DBDC method:
            !
            !         * mrounds     : the maximum number of rounds in one main iteratioin
            !         * mit         : the maximum number of main iterations
            !         * mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * agg_used       : .TRUE. if aggregation is used 
            !         * stepsize_used  : .TRUE. if Simple stepsize determination is not used
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance               !
            !
            ! OUTPUT: * 'f_solution'      : The objective function values for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nk) 
            !         * 'x_solution'      : The vector containing solution vector betas for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nft,nk) 
            !
            ! INOUT:  * 'info'            : The set of information
            !
            !
            ! NOTICE: * The dimension of vectors 'x_solution', 'x_ed' and 'x_koe' has to be 'nft' ('nft' is the dimension of the problem)
            !
            !
            ! NOTICE: * 'print' has to be 0, 1, 2 or 3.           
            !         * 'start' has to be 1, 2, 3 or 4.              
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! INPUTs
               REAL(KIND=c_double), DIMENSION(nft+1), INTENT(IN) :: x_koe    ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nft+1), INTENT(IN) :: x_ed     ! the beta solution for the previous problem where the number of nonzero kits was one smaller   
               
               REAL(KIND=c_double), INTENT(IN) :: f_prev                     ! The objective value for the previous problem where the number of nonzero kits was one smaller       
               
               INTEGER(KIND=c_int), INTENT(IN) :: kit_num_ed                   ! The number of kits in the previous solution
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: mPrNum        ! the starting points for thread

               REAL(KIND=c_double), DIMENSION(8), INTENT(IN) :: mrho      ! Vector containing the values of rho parameter used in the method 
            
               INTEGER(KIND = c_int), INTENT(IN) :: nft                   ! Number of features in x (also the dimension of the problem)
               INTEGER(KIND = c_int), INTENT(IN) :: nk                    ! Number of kits for features
               INTEGER(KIND = c_int), INTENT(IN) :: nrecord               ! Number of observations in data 
               
               INTEGER(KIND = c_int), INTENT(IN) :: start                 ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN) :: iprint                ! Print used

               INTEGER(KIND = c_int), INTENT(IN) :: solver                ! the solver used
               
               INTEGER(KIND = c_int), INTENT(IN) :: k                     ! The number of kits in the problem 
               INTEGER(KIND = c_int), INTENT(IN) :: i                     ! The number of the thread
               
               INTEGER(KIND = c_int), INTENT(IN) :: problem1              ! The DC component f1 of the objective               
               INTEGER(KIND = c_int), INTENT(IN) :: problem2              ! The DC component f2 of the objective               
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds                 ! the maximum number of rounds in one main iteration
                                                                          ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mit                     ! the maximum number of main iterations
                                                                          ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
                                                                          
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds_esc             ! the maximum number of rounds in one escape procedure
                                                                          ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
                
               LOGICAL, INTENT(IN) :: agg_used                            ! .TRUE. if aggregation is used in DBDC
               LOGICAL, INTENT(IN) :: stepsize_used                       ! .TRUE. if Simple stepsize determination is not used
               
              ! Parameters in DBDC method
               INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance

              ! The data matrices
               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               REAL(KIND=c_double), DIMENSION(2,nrecord) :: in_mY         ! observed times and labels matrix (column is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          
                     
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nft+1,nk)  :: x_solution   !Output variable for beta coefficients with k kits and starting points for the thread i
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nk)        :: f_solution   !Output variable target function value with k kits and starting point for the thread i             
           
           !***************************** LOCAL VARIABLES ************************************      
 
               TYPE(set_info) :: set                                    ! The set of information                     

               REAL(KIND=c_double) :: CPUtime                           ! the CPU time (in seconds)
               REAL(KIND=c_double) :: f_solution_DBDC                   ! the solution of DBDC method
  
               REAL(KIND=c_double), DIMENSION(nft+1) :: beta_solution   ! the solution vector beta obtained during the execution      
          
               INTEGER(KIND=c_int), DIMENSION(nk) :: kits_beta          ! indices of kits in the solution 'beta_solution'

               REAL(KIND=c_double), DIMENSION(nft+1) :: x_0             ! the starting point
               
               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                                      
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
               
               LOGICAL :: kit_in_use            ! If .TRUE. kit is used in the solution
               LOGICAL :: run_stop              ! If .TRUE. run of this subroutine is stopped for selected k
               LOGICAL :: mukana                ! If .TRUE. specific kit is in the solution
               
               LOGICAL :: ed_sol_in_pen         ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
               LOGICAL :: new_start             ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
       
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
                            
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                  ! Random number
               REAL(KIND=c_double), DIMENSION(nk) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nk) :: mRandInd     ! Original indices of random numbers in matrix
        
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: kit_num               ! The number of kits in the current solution
               INTEGER(KIND=c_int) :: j, j1, j2, ii, i2
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization  

               
               INTEGER(KIND=c_int) :: startind, j3, kk      ! the start index for starting point    
			   
 !**
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.
 !** 
 
                    !--------------------------------------------        
                    !                INITIALIZATION
                    !--------------------------------------------
                    
                    user_n = nft+1            ! The dimension of the problem
                    iprint_DBDC = 0_c_int     ! The print is supressed 
                  
                    ed_sol_in_pen = .FALSE.
                   !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem  

                   tol_zero = (10.0_c_double)**(-6)
 

                   ! The initialization of parametrs used in DBDC methods
                   CALL allocate_parameters(set, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)
  
                   ! Allocation of sizes od data matrices in function.f95
                    CALL allocate_data_mse(set, nft, nrecord, nk, k)   
                    
                   ! Allocation of data matricse    
                    CALL allocate_matrices_mse(set, in_mX, in_mY, in_mK, in_mC,  & 
                                      & nrecord, nft, nk)
                    
                    CALL  set_k(set, k)                        ! The number of nonzero kits is fixe                 
                    
                    ! Determination of the first starting point for the thread i                   
                    startind = 1
                    IF (i > 1) THEN 
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                    END IF
                    
                    ! Starting points for the thread i are looked through
                    DO kk = startind, startind + mPrNum(i) - 1 
                        !--------------------------------------------
                        !    THE GENERATION OF THE STARTING POINT
                        !--------------------------------------------
                         x_0 = x_ed                  ! The base of the starting point is the previous solution (k-1 nonzero kits)
                         mukana = .FALSE.           
                         IF (start == 2 .OR. start == 3) THEN          ! Only nkits-k+1 starting points used             
                           DO ii = 1, kit_num_ed       ! Then we check if the kit 'kk' belongs to the previous solution
                             IF (kk == kits_beta_ed(ii)) THEN   
                               mukana = .TRUE.         ! 'mukana' tells if the kit 'kk' is in solution     
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
                                             
                         
                       IF (new_start) THEN       ! .TRUE. if a new starting point is generated
                        
                        IF ((start == 2) .OR. (start == 3)) THEN 
                          DO ii = 1, nft               ! the features are looked through
                            IF (set%mK(ii,kk)==1) THEN
                               x_0(ii) = x_koe(ii)     ! In starting poitnt x_0, we initialize features in kit 'kk' with values from solution x_koe
                            END IF   
                          END DO
                        END IF  
                        
                        IF (start == 4) THEN ! A random starting point is generated
                           x_0 = 0.01_c_double
                           DO ii = 1, nk
                             CALL RANDOM_NUMBER(random_num)
                             mRand(ii) = random_num
                             mRandInd(ii) = ii
                           END DO 
                           CALL heapsort_ind(mRand,mRandInd)               
                           DO ii = 1, k
                              DO i2 = 1, nft 
                                IF (set%mK(i2,mRandInd(ii))==1) THEN
                                  x_0(i2) = x_koe(i2)    ! In starting point x_0, we initialize features in kit i with values from solution x_koe
                                END IF   
                              END DO     
                           END DO
                        END IF                  
                          
                        run_stop = .FALSE.          ! Initialization of 'run_stop' to .FALSE. since we cannot stop
                        num_rho = 0                 ! The selection of the first value of penalization parameter rho
                      
                        IF (iprint >= 2) THEN
                           WRITE(*,*) '-------------------------------------' 
                           WRITE(*,*) 'New start point used.' 
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
                          rho = 10.0_c_double * rho
                          
                          IF (num_rho>=24) THEN 
                          run_stop = .TRUE.         ! Forced stop in optimization
                          END IF 
                        END IF
                     
            
                        ! The optimization problem is solved fot the current rho and x_0
!**						
						IF (solver==1) THEN 
         
                            CALL DBDC_algorithm( f_solution_DBDC, beta_solution, x_0, rho, 0.0_c_double, &
                                & mit, mrounds, mrounds_esc, termination, counter, CPUtime,  &
                                & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n, &
                                & set) 
								
						ELSE IF (solver == 2) THEN
						
                           CALL init_LMBMinfo(problem1, set)   
                           CALL init_x(x_0)
                           CALL set_rho_LMBM(rho)                 
                           CALL set_lambda_LMBM(0.0_c_double)              
                           CALL lmbm(mc, f_solution_DBDC, iout(1),iout(2),iout(3),iout(4),LMBMstart)      
                           CALL copy_x(beta_solution)   
						   
                        END IF
!**										
								
                
                        IF (ed_sol_in_pen) THEN 
                          x_0 = beta_solution   ! Starting point for the next round
                        END IF    
                    
                        help_counter = 0           
                        DO j = 1, nft         ! The number of zero components in beta vector is computed
                          IF ( ABS(beta_solution(j)) < tol_zero) THEN
                            help_counter = help_counter +1  
                            beta_solution(j) = 0.0_c_double
                          END IF
                        END DO
                  
                        cost = 0.0_c_double      ! The initialization of the cost of 'beta_solution'
                        kit_num = 0              ! The initialization of the number of kits in 'beta_solution'
                        kits_beta = 0            ! The initialization of kits in 'beta_solution'
                      
                        !Calculation of kits in solution
                        DO j1 = 1, nk        ! Each kit is looked through
                          kit_in_use = .FALSE.  
                          DO j2 = 1, nft
                            IF (set%mK(j2,j1)==1) THEN 
                              IF ( ABS(beta_solution(j2)) >= tol_zero) THEN    ! If the condition is satisfied the kit j1 is in solution 'beta_solution'
                                 kit_in_use = .TRUE.                                  
                              END IF
                            END IF
                          END DO
                          IF (kit_in_use) THEN       ! Executed if kit j1 is in solution 'beta_solution'
                            cost = cost + set%mC(j1)     ! The cost of kit j1 is taken into account
                            kit_num = kit_num + 1    ! The number of kits is updated
                            kits_beta(kit_num) = j1  ! The index of kit j1 is updated to table kits_beta                      
                          END IF 
                        END DO
        
                        ! We check if the optimization of problem 3 with k kits can be stopped            
                        IF (kit_num <= k .OR. run_stop) THEN 
                           IF (run_stop) THEN 
                             f_solution(kk) = (10.0_c_double)**10
                           END IF
                           run_stop = .TRUE.     ! The optimization can be stopped
                           DO j = 1, user_n
                             x_solution(j,kk) = beta_solution(j)        ! We store the solution to x_solution
                           END DO       
                           f1_current = f1(set,beta_solution,problem1,user_n)
                           f2_current = f2(set,beta_solution,problem2,user_n)
                           f_solution(kk) = f1_current-f2_current       ! We store the objective funtion value without penalization to f_solution
                        END IF               
                       
                        f1_current = f1(set,beta_solution,problem1,user_n) 
                        f2_current = f2(set,beta_solution,problem2,user_n)
        
                        IF (iprint > 2) THEN
                           WRITE(*,*) 'rho', rho, 'f',f1_current-f2_current, 'kits', kit_num   
                        END IF
                        
                        IF (run_stop) THEN
                          IF ((iprint >= 2) .AND. (kit_num <= k)) THEN   
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num
                          END IF  
                          IF ((iprint >=2) .AND. (kit_num > k)) THEN                                      
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num, 'should be equal to', k 
                          END IF
                        END IF
                        
                        END DO
                        
                       ELSE
                       
                         DO j = 1, user_n
                            x_solution(j,kk) = x_ed(j)
                         END DO     
                         f_solution(kk) = f_prev
                       END IF           
                       
                   END DO           
                       
             
                   CALL deallocate_data_mse(set)             
             
         
         END SUBROUTINE  solution_with_k_kits_mse
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|




       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |             SOLVING THE PROBLEM WITH k KITS AND STARTING POINT i               |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
         SUBROUTINE solution_with_k_kits_log( x_solution, f_solution, k, i, f_prev, & 
                                    & x_ed, x_koe, kit_num_ed, kits_beta_ed, &
                                    & nk, start, iprint, mrho, mit, mrounds, mrounds_esc, &
                                    & agg_used, stepsize_used, nft, problem1, problem2,  &      
                                    & in_mX, in_mY, in_mK, in_mC, nrecord, & 
                                    & in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                    & in_b, in_m_clarke, in_eps, in_crit_tol, mPrNum, solver) 
            !_____________________________________________________________________________________
            ! 
            !           
            ! * 'problem1  and problem2   : defines the used objective function is the logistic regression model with L0-norm for kits
            !                              (parametric approach since the number of nonzero elements is fixed) 
            !
            !   Solves the L0-norm problem from a fixed starting points for thread "i" with "k" kits
            !                                                 
            !
            ! INPUT:  * 'nft'           : Number of rows. The dimension of the problem = nft + 1, INTEGER
            !         * 'nk'            : the (maximum) number of kits, INTEGER
            !         * 'nrecord'       : the number of observations, INTEGER
            !
            !         * 'k'             : the number of kits in the problem, INTEGER
            !         * 'i'             : the number of thread, INTEGER
            !
            !         * 'x_ed'          : the beta solution for the previous problem where the number of nonzero kits was one smaller (INTEGER)
            !         * 'x_koe'         : The solution to Cox's proportional hazard model without regularization
            !         * 'kit_num_ed'    : The number of kits in the previous solution
            !         * 'kit_beta_ed'   : indices of kits in the previous solution 'x_ed'
            ! 
            !         * 'print'         : specifies the print, INTEGER
            !         * 'start'         : specifies how starting points are selected when the L0-norm problem is solved for  
            !                             the fixed number of kits and in which order the L0-norm problems are solved, INTEGER
            !         * 'mrho'          : List of used rho values
            !
            !         * 'problem1'      : Defines the DC component f1 of the objective
            !         * 'problem2'      : Defines the DC component f2 of the objective
            !
            !         * 'in_mX', 'in_mY', 'in_mK', 'in_mC'    : Are the data matrices 
            !         * 'mPrNum'        : the information about starting points for the thread i
            !         
            !         * 'solver'   : defines the solver used (1=DBDC, 2=LMBM)
            !         
            !         PARAMETERS in DBDC method:
            !
            !         * mrounds     : the maximum number of rounds in one main iteratioin
            !         * mit         : the maximum number of main iterations
            !         * mrounds_esc : the maximum number of rounds in escape procedure
            !
            !         * agg_used       : .TRUE. if aggregation is used 
            !         * stepsize_used  : .TRUE. if Simple stepsize determination is not used
            !
            !         * in_b1          : the size of bundle B1
            !         * in_b2          : the size of bundle B2
            !         * in_b           : the size of bundle in escape procedure
            !
            !         * in_m           : the descent parameter in main iteration
            !         * in_m_clarke    : the descent parameter in escape procedure
            !         * in_c           : the extra decrease parameter in main iteration
            !         * in_r_dec       : the decrease parameter in main iteration
            !         * in_r_inc       : the increase parameter in main iteration
            !         * in_eps1        : the enlargement parameter
            !         * in_eps         : the stopping tolerance: proximity measure  
            !         * in_crit_tol    : the stopping tolerance: criticality tolerance               !
            !
            ! OUTPUT: * 'f_solution'      : The objective function values for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nk) 
            !         * 'x_solution'      : The vector containing solution vector betas for the L0-norm problem with "k" kits from fixed strating points in thread "i", 
            !                              REAL, DIMENSION(nft,nk) 
            !
            ! INOUT:  * 'info'            : The set of information
            !
            !
            ! NOTICE: * The dimension of vectors 'x_solution', 'x_ed' and 'x_koe' has to be 'nft' ('nft' is the dimension of the problem)
            !
            !
            ! NOTICE: * 'print' has to be 0, 1, 2 or 3.           
            !         * 'start' has to be 1, 2, 3 or 4.              
            
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER (INPUT/OUTPUT) *************************************   
            ! INPUTs
               REAL(KIND=c_double), DIMENSION(nft+1), INTENT(IN) :: x_koe    ! The solution to Cox's proportional hazard model without regularization
               REAL(KIND=c_double), DIMENSION(nft+1), INTENT(IN) :: x_ed     ! the beta solution for the previous problem where the number of nonzero kits was one smaller   
               
               REAL(KIND=c_double), INTENT(IN) :: f_prev                   ! The objective value for the previous problem where the number of nonzero kits was one smaller       
               
               INTEGER(KIND=c_int), INTENT(IN) :: kit_num_ed                   ! The number of kits in the previous solution
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: kits_beta_ed  ! indices of kits in the previous solution 'x_ed'
               INTEGER(KIND=c_int), DIMENSION(nk), INTENT(IN) :: mPrNum        ! the starting points for thread

               REAL(KIND=c_double), DIMENSION(8), INTENT(IN) :: mrho      ! Vector containing the values of rho parameter used in the method 
            
               INTEGER(KIND = c_int), INTENT(IN) :: nft                   ! Number of features in x (also the dimension of the problem)
               INTEGER(KIND = c_int), INTENT(IN) :: nk                    ! Number of kits for features
               INTEGER(KIND = c_int), INTENT(IN) :: nrecord               ! Number of observations in data 
               
               INTEGER(KIND = c_int), INTENT(IN) :: start                 ! Starting point procedure used
               INTEGER(KIND = c_int), INTENT(IN) :: iprint                ! Print used

               INTEGER(KIND = c_int), INTENT(IN) :: solver                ! the solver used
               
               INTEGER(KIND = c_int), INTENT(IN) :: k                     ! The number of kits in the problem 
               INTEGER(KIND = c_int), INTENT(IN) :: i                     ! The number of the thread
               
               INTEGER(KIND = c_int), INTENT(IN) :: problem1              ! The DC component f1 of the objective               
               INTEGER(KIND = c_int), INTENT(IN) :: problem2              ! The DC component f2 of the objective               
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds                 ! the maximum number of rounds in one main iteration
                                                                          ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
               
               INTEGER(KIND=c_int), INTENT(INOUT) :: mit                     ! the maximum number of main iterations
                                                                          ! If 'mit' <=0 then DEFAULT value 'mit'=5000 is used
                                                                          
               INTEGER(KIND=c_int), INTENT(INOUT) :: mrounds_esc             ! the maximum number of rounds in one escape procedure
                                                                          ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
                
               LOGICAL, INTENT(IN) :: agg_used                            ! .TRUE. if aggregation is used in DBDC
               LOGICAL, INTENT(IN) :: stepsize_used                       ! .TRUE. if Simple stepsize determination is not used
               
              ! Parameters in DBDC method
               INTEGER(KIND=c_int), INTENT(IN) :: in_b1                   ! the size of bundle B1
               INTEGER(KIND=c_int), INTENT(IN) :: in_b2                   ! the size of bundle B2
               INTEGER(KIND=c_int), INTENT(IN) :: in_b                    ! the size of bundle in escape procedure
               
               REAL(KIND=c_double), INTENT(IN) :: in_m                    ! the descent parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_m_clarke             ! the descent parameter in escape procedure
               REAL(KIND=c_double), INTENT(IN) :: in_c                    ! the extra decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_dec                ! the decrease parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_r_inc                ! the increase parameter in main iteration
               REAL(KIND=c_double), INTENT(IN) :: in_eps1                 ! the enlargement parameter
               REAL(KIND=c_double), INTENT(IN) :: in_eps                  ! the stopping tolerance: proximity measure  
               REAL(KIND=c_double), INTENT(IN) :: in_crit_tol             ! the stopping tolerance: criticality tolerance

              ! The data matrices
               REAL(KIND=c_double), DIMENSION(nft,nrecord) :: in_mX       ! predictor matrix (column is an observation)
               INTEGER(KIND=c_int), DIMENSION(2,nrecord) :: in_mY         ! observed times and labels matrix (column is an observation)  
               INTEGER(KIND=c_int), DIMENSION(nk,nft) :: in_mK            ! kit matrix (row is a kit)
               REAL(KIND=c_double), DIMENSION(nk) :: in_mC                ! kit costs          
                     
            ! OUTPUTs
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nft+1,nk)  :: x_solution   !Output variable for beta coefficients with k kits and starting points for the thread i
               REAL(KIND = c_double), INTENT(OUT), DIMENSION(nk)        :: f_solution   !Output variable target function value with k kits and starting point for the thread i
               
           
           !***************************** LOCAL VARIABLES ************************************      
 
               TYPE(set_info) :: set                                    ! The set of information                     

               REAL(KIND=c_double) :: CPUtime                           ! the CPU time (in seconds)
               REAL(KIND=c_double) :: f_solution_DBDC                   ! f value at the solution of DBDC method
  
               REAL(KIND=c_double), DIMENSION(nft+1) :: beta_solution     ! the solution vector beta obtained during the execution      
          
               INTEGER(KIND=c_int), DIMENSION(nk) :: kits_beta          ! indices of kits in the solution 'beta_solution'

               REAL(KIND=c_double), DIMENSION(nft+1) :: x_0             ! the starting point
               
               INTEGER(KIND=c_int) :: termination         ! The reason for termination in DBDC method
                                                          ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER(KIND=c_int), DIMENSION(8) :: counter   ! contains the values of different counteres for DBDC method: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'             
        
               INTEGER(KIND=c_int) :: user_n             ! the dimension of the problem
                                      
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
               
               LOGICAL :: kit_in_use            ! If .TRUE. kit is used in the solution
               LOGICAL :: run_stop              ! If .TRUE. run of this subroutine is stopped for selected k
               LOGICAL :: mukana                ! If .TRUE. specific kit is in the solution
               
               LOGICAL :: ed_sol_in_pen         ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
               LOGICAL :: new_start             ! If .TRUE. previous solution is utilized during the solution of the penalized problem 
       
               REAL(KIND=c_double) :: rho             ! The parameter rho used in L0-norm
               REAL(KIND=c_double) :: cost            ! The cost of solution beta
                            
               REAL(KIND=c_double) :: f1_current      ! The value of f1
               REAL(KIND=c_double) :: f2_current      ! The value of f2
               
               REAL(KIND=c_double) :: tol_zero        ! The tolerance for value zero (i.e. if value is smaller than 'tol_zero' -> it is set to be zero
               
               REAL(KIND=c_double) :: random_num                  ! Random number
               REAL(KIND=c_double), DIMENSION(nk) :: mRand        ! Random number matrix
               INTEGER(KIND=c_int), DIMENSION(nk) :: mRandInd     ! Original indices of random numbers in matrix
        
               INTEGER(KIND=c_int) :: help_counter
               INTEGER(KIND=c_int) :: num_rho
               INTEGER(KIND=c_int) :: kit_num               ! The number of kits in the current solution
               INTEGER(KIND=c_int) :: j, j1, j2, ii, i2
               INTEGER(KIND=c_int) :: max_threads           ! the maximum number of threads that can be used in parallellization         
               INTEGER(KIND=c_int) :: startind, j3, kk      ! the start index for starting point    

 !**
               ! Scalars in LMBM
               INTEGER(KIND=c_int) ::  mc        ! Number of corrections.
               REAL(KIND=c_double) :: LMBMstart  ! The starting time
               INTEGER, DIMENSION(4) :: iout     ! Output integer parameters.
                                                 !   iout(1)   Number of used iterations.
                                                 !   iout(2)   Number of used function evaluations.
                                                 !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                 !               or number of outer iterations (LDGBM).
                                                 !   iout(4)   Cause of termination:
                                                 !               1  - The problem has been solved
                                                 !                    with desired accuracy.
                                                 !               2  - Changes in function values < tolf in mtesf
                                                 !                    subsequent iterations.
                                                 !               3  - Changes in function value < 
                                                 !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                 !                    where small is the smallest positive number
                                                 !                    such that 1.0 + small > 1.0.
                                                 !               4  - Number of function calls > mfe.
                                                 !               5  - Number of iterations > mit.
                                                 !               6  - Time limit exceeded. 
                                                 !               7  - f < tolb.
                                                 !              -1  - Two consecutive restarts.
                                                 !              -2  - Number of restarts > maximum number
                                                 !                    of restarts.
                                                 !              -3  - Failure in function or subgradient
                                                 !                    calculations (assigned by the user).
                                                 !              -4  - Failure in attaining the demanded
                                                 !                    accuracy.
                                                 !              -5  - Invalid input parameters.
                                                 !              -6  - Unspecified error.
 !** 
 
                    !--------------------------------------------        
                    !                INITIALIZATION
                    !--------------------------------------------
                    
                    user_n = nft+1            ! The dimension of the problem
                    iprint_DBDC = 0_c_int     ! The print is supressed 
                  
                    ed_sol_in_pen = .FALSE.
                   !ed_sol_in_pen = .TRUE.    ! the previous solution is utilized during the solution of the penalized problem  

                   tol_zero = (10.0_c_double)**(-6)
 

                   ! The initialization of parametrs used in DBDC methods
                   CALL allocate_parameters(set, in_b1, in_b2, in_m, in_c, in_r_dec, in_r_inc, in_eps1, &
                                           & in_b, in_m_clarke, in_eps, in_crit_tol)
  
                   ! Allocation of sizes od data matrices in function.f95
                    CALL allocate_data_log(set, nft, nrecord, nk, k)   
                    
                   ! Allocation of data matricse    
                    CALL allocate_matrices_log(set, in_mX, in_mY, in_mK, in_mC,  & 
                                      & nrecord, nft, nk)
                    
                    CALL  set_k(set, k)                        ! The number of nonzero kits is fixed

                    ! Determination of the first starting point for the thread i                   
                    startind = 1
                    IF (i > 1) THEN 
                       DO j3 = 1, i-1
                          startind = startind + mPrNum(j3)
                       END DO   
                    END IF
                    
                    ! Starting points for the thread i are looked through
                    DO kk = startind, startind + mPrNum(i) - 1                     
                        !--------------------------------------------
                        !    THE GENERATION OF THE STARTING POINT
                        !--------------------------------------------
                         x_0 = x_ed                  ! The base of the starting point is the previous solution (k-1 nonzero kits)
                         mukana = .FALSE.           
                         IF (start == 2 .OR. start == 3) THEN          ! Only nkits-k+1 starting points used             
                           DO ii = 1, kit_num_ed       ! Then we check if the kit 'kk' belongs to the previous solution
                             IF (kk == kits_beta_ed(ii)) THEN   
                               mukana = .TRUE.         ! 'mukana' tells if the kit 'kk' is in solution     
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
                                             
                         
                       IF (new_start) THEN       ! .TRUE. if a new starting point is generated
                        
                        IF ((start == 2) .OR. (start == 3)) THEN 
                          DO ii = 1, nft               ! the features are looked through
                            IF (set%mK(ii,kk)==1) THEN
                               x_0(ii) = x_koe(ii)     ! In starting poitnt x_0, we initialize features in kit 'kk' with values from solution x_koe
                            END IF   
                          END DO
                        END IF  
                        
                        IF (start == 4) THEN ! A random starting point is generated
                           x_0 = 0.01_c_double
                           DO ii = 1, nk
                             CALL RANDOM_NUMBER(random_num)
                             mRand(ii) = random_num
                             mRandInd(ii) = ii
                           END DO 
                           CALL heapsort_ind(mRand,mRandInd)                   
                           DO ii = 1, k
                              DO i2 = 1, nft 
                                IF (set%mK(i2,mRandInd(ii))==1) THEN
                                  x_0(i2) = x_koe(i2)    ! In starting point x_0, we initialize features in kit i with values from solution x_koe
                                END IF   
                              END DO     
                           END DO
                        END IF                  
                          
                        run_stop = .FALSE.          ! Initialization of 'run_stop' to .FALSE. since we cannot stop
                        num_rho = 0                 ! The selection of the first value of penalization parameter rho
                      
                        IF (iprint >= 2) THEN
                           WRITE(*,*) '-------------------------------------' 
                           WRITE(*,*) 'New start point used.' 
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
                          rho = 10.0_c_double * rho
                          
                          IF (num_rho>=24) THEN 
                          run_stop = .TRUE.         ! Forced stop in optimization
                          END IF 
                        END IF
                    
                        ! The optimization problem is solved fot the current rho and x_0
         
!**						
						IF (solver==1) THEN 
         
                            CALL DBDC_algorithm( f_solution_DBDC, beta_solution, x_0, rho, 0.0_c_double, &
                                & mit, mrounds, mrounds_esc, termination, counter, CPUtime,  &
                                & agg_used, stepsize_used, iprint_DBDC, problem1, problem2, user_n, &
                                & set) 
								
						ELSE IF (solver == 2) THEN
						
                           CALL init_LMBMinfo(problem1, set)   
                           CALL init_x(x_0)
                           CALL set_rho_LMBM(rho)                 
                           CALL set_lambda_LMBM(0.0_c_double)              
                           CALL lmbm(mc, f_solution_DBDC, iout(1),iout(2),iout(3),iout(4),LMBMstart)      
                           CALL copy_x(beta_solution)   
						   
                        END IF
!**		 
 
                        IF (ed_sol_in_pen) THEN 
                          x_0 = beta_solution   ! Starting point for the next round
                        END IF    
                    
                        help_counter = 0           
                        DO j = 1, nft         ! The number of zero components in beta vector is computed
                          IF ( ABS(beta_solution(j)) < tol_zero) THEN
                            help_counter = help_counter +1
                            beta_solution(j) = 0.0_c_double
                          END IF
                        END DO
                  
                        cost = 0.0_c_double      ! The initialization of the cost of 'beta_solution'
                        kit_num = 0              ! The initialization of the number of kits in 'beta_solution'
                        kits_beta = 0            ! The initialization of kits in 'beta_solution'
                      
                        !Calculation of kits in solution
                        DO j1 = 1, nk        ! Each kit is looked through
                          kit_in_use = .FALSE.  
                          DO j2 = 1, nft
                            IF (set%mK(j2,j1)==1) THEN 
                              IF ( ABS(beta_solution(j2)) >= tol_zero) THEN    ! If the condition is satisfied the kit j1 is in solution 'beta_solution'
                                 kit_in_use = .TRUE.                                  
                              END IF
                            END IF
                          END DO
                          IF (kit_in_use) THEN       ! Executed if kit j1 is in solution 'beta_solution'
                            cost = cost + set%mC(j1)     ! The cost of kit j1 is taken into account
                            kit_num = kit_num + 1    ! The number of kits is updated
                            kits_beta(kit_num) = j1  ! The index of kit j1 is updated to table kits_beta                      
                          END IF 
                        END DO
        
                        ! We check if the optimization of problem 3 with k kits can be stopped            
                        IF (kit_num <= k .OR. run_stop) THEN 
                           IF (run_stop) THEN 
                             f_solution(kk) = (10.0_c_double)**10
                           END IF
                           run_stop = .TRUE.     ! The optimization can be stopped
                           DO j = 1, user_n
                             x_solution(j,kk) = beta_solution(j)        ! We store the solution to x_solution
                           END DO       
                           f1_current = f1(set,beta_solution,problem1,user_n)
                           f2_current = f2(set,beta_solution,problem2,user_n)
                           f_solution(kk) = f1_current-f2_current       ! We store the objective funtion value without penalization to f_solution
                        END IF               
                       
                        f1_current = f1(set,beta_solution,problem1,user_n) 
                        f2_current = f2(set,beta_solution,problem2,user_n)
        
                        IF (iprint > 2) THEN
                           WRITE(*,*) 'rho', rho, 'f',f1_current-f2_current, 'kits', kit_num   
                        END IF
                        
                        IF (run_stop) THEN
                          IF ((iprint >= 2) .AND. (kit_num <= k)) THEN   
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num
                          END IF  
                          IF ((iprint >=2) .AND. (kit_num > k)) THEN                                      
                            WRITE(*,*)
                            WRITE(*,*) 'f=',f1_current-f2_current, 'and kits', kit_num, 'should be equal to', k 
                          END IF
                        END IF
                        
                        END DO
                        
                       ELSE
                         DO j = 1, user_n
                            x_solution(j,kk) = x_ed(j)
                         END DO     
                         f_solution(kk) = f_prev
                       END IF        

                   END DO                      
             
                   CALL deallocate_data_log(set)             
             
         
         END SUBROUTINE  solution_with_k_kits_log
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|


        
  
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
  

  END MODULE oscar





