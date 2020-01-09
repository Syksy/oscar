module blasso
    use, intrinsic :: iso_c_binding
 
    implicit none
    private
    public :: testblasso
    public :: blassomat
    public :: blassocoxfit
 
contains
 
    subroutine testblasso(x, n, res) bind(C, name = "testblasso_f_")
 
        integer(kind = c_int), intent(in), value        :: n    !Length of x
        real(kind = c_double), intent(in), dimension(n) :: x    !Vector of values
        real(kind = c_double), intent(out)              :: res  !Output variable
        integer                                         :: i    !Internal count
 
        res = 0.0_c_double
        do i = 1, n
            res = res + x(i)
        end do
 
    end subroutine testblasso
 
    ! Subroutine testing for computing column averages for an input matrix from R
    subroutine blassomat(x, nrow, ncol, res) bind(C, name = "blassomat_f_")
 
        integer(kind = c_int), intent(in), value		:: nrow 	!Number of rows in x
        integer(kind = c_int), intent(in), value		:: ncol 	!Number of cols in x
        real(kind = c_double), intent(in), dimension(nrow*ncol) :: x    	!Vector of values
        real(kind = c_double), intent(out), dimension(ncol)     :: res  	!Output variable
        integer                                         	:: i,j,ind	!Internal count
	real(kind = c_double), dimension(nrow,ncol)		:: xmat 	!Matrix inside Fortran 
 
 	write(*,*) 'Function starts'
 	write(*,*) nrow, ncol
 	write(*,*) x
 
 	! Populate the matrix xmat
        ind = 0
	DO j = 1, ncol
          DO i = 1, nrow
            ind = ind +1
            xmat(i,j) = x(ind)
          END DO
        END DO

 	write(*,*) 'Matrix populated'

	! Compute column averages as an example of a matrix operation
        res = 0.0_c_double
	DO j=1, ncol
	  DO i=1, nrow
	    res(j) = res(j) + xmat(i,j)
	  END DO
	  res(j) = res(j) / nrow
	END DO
 
  	write(*,*) 'End of Fortran subroutine'

    end subroutine blassomat

 
    ! Subroutine placeholder for actual budget-aware LASSO
    subroutine blassocoxfit(yevent, ytime, x, xnrow, xncol, k, knrow, kncol, costvec, beta) bind(C, name = "blassocoxfit_f_")

	! Data input (response and data)
	integer(kind = c_int), intent(in), dimension(xnrow)		:: yevent	!Integer of death event (1) or censoring (0)
	integer(kind = c_int), intent(in), dimension(xnrow)		:: ytime	!Integer of follow-up times (presumably in full days/months/etc)
        real(kind = c_double), intent(in), dimension(xnrow*xncol)	:: x    	!Vector of data values
        integer(kind = c_int), intent(in), dimension(knrow*kncol) 	:: k    	!Vector of indicator matrix values
        real(kind = c_double), intent(in), dimension(xncol)		:: costvec    	!Vector of data values
	! Helper variables, i.e. populating x and k matrices
        integer(kind = c_int), intent(in), value			:: xnrow 	!Number of rows in data matrix X
        integer(kind = c_int), intent(in), value			:: xncol 	!Number of cols in data matrix X
        integer(kind = c_int), intent(in), value			:: knrow 	!Number of rows in kit indicator matrix K
        integer(kind = c_int), intent(in), value			:: kncol 	!Number of cols in kit indicator matrix K
	! Variables for internal use in Fortran
	real(kind = c_double), dimension(xnrow,xncol)			:: xmat 	!Data matrix inside Fortran 
	integer(kind = c_int), dimension(knrow,kncol)			:: kmat		!Indicator matrix inside Fortran 
        integer                                         		:: i,j,ind	!Internal variables
	! Output back to C and R
        real(kind = c_double), intent(out), dimension(ncol)     	:: beta  	!Output variable
 
 	write(*,*) 'Function starts'
 	!write(*,*) nrow, ncol
 	!write(*,*) x
 
 	! Populate the matrix xmat
        ind = 0
	DO j = 1, xncol
          DO i = 1, xnrow
            ind = ind +1
            xmat(i,j) = x(ind)
          END DO
        END DO
	! Populate the indicator matrix kmat
	ind = 0
	DO j = 1, kncol
          DO i = 1, knrow
            ind = ind +1
            kmat(i,j) = k(ind)
          END DO
        END DO
	

 	write(*,*) 'Matrix populated'

	!!!!
	!
	! Kaisa's algorithm procedures here!
	! Ought to produce suitable beta fit based on input variables
	!
	!!!!

	! Compute column averages as an example of a matrix operation (here beta would be the fitted model coefficients)
        beta = 0.0_c_double
	DO j=1, xncol
	  DO i=1, xnrow
	    beta(j) = beta(j) + xmat(i,j)
	  END DO
	  beta(j) = beta(j) / xnrow
	END DO
 
  	write(*,*) 'End of Fortran subroutine'

    end subroutine blassocoxfit
  
end module blasso
