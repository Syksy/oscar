module blasso
    use, intrinsic :: iso_c_binding
 
    implicit none
    private
    public :: blassocox
 
contains
 
    ! Subroutine for blasso
    subroutine blassocox(x, y, kits, costs, nrow, ncol, nkits, beta, fperk) &
    	& bind(C, name = "blassocox_f_")

	!!!
	!
 	! In/Out from R and C
 	!
 	! INPUT
 	!
        integer(kind = c_int), intent(in), value		:: nrow 	!Number of rows in x (i.e. records)
        integer(kind = c_int), intent(in), value		:: ncol 	!Number of cols in x (i.e. features)
        integer(kind = c_int), intent(in), value		:: nkits 	!Number of kits for features
        real(kind = c_double), intent(in), dimension(nrow*ncol) :: x    	!Vector of data values
        real(kind = c_double), intent(in), dimension(nrow*2) 	:: y    	!Vector of response values (2-column survival, time + event)
        integer(kind = c_int), intent(in), dimension(nkits*ncol):: kits		!Vector of kit indicator values (binary indicators)
        real(kind = c_double), intent(in), dimension(nkits)	:: costs	!Costs associated to each kit
        ! OUTPUT
        real(kind = c_double), intent(out), dimension(ncol)     :: beta  	!Output variable for beta coefficients
        real(kind = c_double), intent(out), dimension(ncol*nkits):: fperk  	!Output variable target function value per k
	!
	!!!

        !!!
        !
        ! Local variables
        !
        integer                                         	:: i,j,ind	! Internal count
	real(kind = c_double), dimension(nrow,ncol)		:: xmat 	! Data matrix
	integer(kind = c_int), dimension(nrow,2)		:: ymat		! Survival response matrix/vector
 	integer(kind = c_int), dimension(nkits,ncol)		:: kmat		! Kit indicator matrix
 	real(kind = c_double), dimension(nkits)			:: cvec		! Cost vector for kits
	!
	!!!
 
 	write(*,*) 'Function starts'
 	write(*,*) nrow, ncol
 	write(*,*) x
  	write(*,*) y
  	write(*,*) kits
  	
 	! Populate the matrix xmat of dim {nrow,ncol}
        ind = 0
	DO j = 1, ncol
          DO i = 1, nrow
            ind = ind +1
            xmat(i,j) = x(ind)
          END DO
        END DO
        ! Populate response matrix of dim {nrow,2}
	ind = 0
	DO j = 1, 2
	  DO i = 1, nrow
	    ind = ind + 1
	    ymat(i,j) = y(ind)
	  END DO
	END DO
        ! Populate response matrix of dim {nkits,ncol}
	ind = 0
	DO j = 1, ncol
	  DO i = 1, nkits
	    ind = ind + 1
	    kmat(i,j) = kits(ind)
	  END DO
	END DO
	! Populate the cost vector for kits
	ind = 0
	DO i = 1, nkits
	  ind = ind + 1
	  cvec(i) = costs(ind)
	END DO

 	write(*,*) 'Populated matrix values:'
 	write(*,*) 'xmat:'
 	write(*,*) xmat
 	write(*,*) 'ymat:'
 	write(*,*) ymat
 	write(*,*) 'kmat:'
 	write(*,*) kmat
 	write(*,*) 'cvec:'
 	write(*,*) cvec

	! Format variables for coefficients, function value per k, etc
        beta = 1.0_c_double ! Test non-zero double formatting
        fperk = 0.0_c_double

	! 
	! ... ACTUAL ALGORITHM ...
	!
 
  	write(*,*) 'End of Fortran subroutine'

    end subroutine blassocox

end module blasso
