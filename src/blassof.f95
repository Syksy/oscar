module blasso
    use, intrinsic :: iso_c_binding
 
    implicit none
    private
    public :: testblasso
    public :: blassomat
 
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
  
end module blasso
