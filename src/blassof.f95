module blasso
    use, intrinsic :: iso_c_binding
 
    implicit none
    private
    public :: testblasso
 
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
 
end module blasso
