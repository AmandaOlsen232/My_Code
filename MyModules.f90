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

function test(x) result(y)
    real(dp), dimension(:,:), allocatable :: y
    real(dp), dimension(:,:), intent(in) :: x
    allocate(y(3,1))
    y(1,1) = x(1,1) + x(2,1) + x(3,1)
    y(2,1) = x(1,1)*2 + x(2,1)*2 + x(3,1)*2
    y(3,1) = x(1,1)**2 + x(2,1)**2 + x(3,1)**2
end function test

function jacobian(f, x, h_offset) result(J)
    procedure(jacobian_function) :: f !function names
    real(dp), dimension(:,:), intent(in) :: x !nx1 array describing current state
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
implicit none

!Jacobian Function example using the test function and 3x1 array
!Prints out the Jacobian to the console
real(dp), dimension(:,:), allocatable :: J, x
integer :: i

allocate(x(3,1))

x = reshape([1.0, 2.0, 3.0], shape(x))
J = jacobian(test, x)

do i=1, 3
    write(*,*) J(i,:)
end do



end program Main
