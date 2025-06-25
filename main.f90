module my_functions
implicit none
integer, parameter :: dp = selected_real_kind(16, 307) !double precision

abstract interface
    function jacobian_function(x) result(y)
        import dp
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(:,:), allocatable :: y 
    end function jacobian_function
end interface

contains
function jacobian(f, x, h_offset) result(J)
    procedure(jacobian_function) :: f
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), allocatable :: J, eps
    real(dp), intent(in), optional :: h_offset
    real(dp), dimension(:,:), allocatable :: fp, fm
    integer :: k, i, dim_x, dim_fx
    real(dp) :: h

    !set default value for h
    if (present(h_offset)) then
        h = h_offset 
    else
        h = 1e-3
    end if

    !get the dimension of input and output vectors, allocate jacobian size
    dim_x = size(x)
    dim_fx = size(f(x))
    allocate(J(dim_fx, dim_x))
    allocate(eps(dim_x, 1))
    J = 0
    eps = 0

    !numerically solve for the jacobian
    fx: do k=1, dim_fx
        in: do i=1, dim_x
            eps = 0
            eps(i,1) = h
            fp = f(x+eps)
            fm = f(x-eps)
            J(k,i) = (fp(k,1) - fm(k,1))/(2*h)
        end do in
    end do fx
    
end function jacobian

end module my_functions



program Main
use my_functions




end program Main
