module my_functions
    use math_m
implicit none


abstract interface
    function jacobian_function(x) result(y)
        real, dimension(:), intent(in) :: x
        real, dimension(:), allocatable :: y 
    end function jacobian_function
end interface

contains

function test(x) result(y)
    real, dimension(:), allocatable :: y
    real, dimension(:), intent(in) :: x
    allocate(y(3))
    y(1) = x(1) + x(2) + x(3)**3
    y(2) = x(2)**2 + x(3)
    y(3) = x(3)**3 + x(2)**2 + x(1)**3
end function test



function jacobian(f, x, h_offset) result(J)
    procedure(jacobian_function) :: f !function names
    real, dimension(:), intent(in) :: x !nx1 array describing current state
    real, dimension(:), allocatable :: eps
    real, dimension(:,:), allocatable :: J
    real, intent(in), optional :: h_offset
    real, dimension(:), allocatable :: fp, fm
    integer :: k, i, dim_x, dim_fx
    real :: h

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
    allocate(eps(dim_x))
    J = 0
    eps = 0

    !numerically solve for the jacobian
    fx: do k=1, dim_fx
        in: do i=1, dim_x
            eps = 0
            eps(i) = h
            fp = f(x+eps)
            fm = f(x-eps)
            J(k,i) = (fp(k) - fm(k))/(2*h)
        end do in
    end do fx
    
end function jacobian

function multivariable_newtons_method(f, x) result(G)
    implicit none 
    procedure(jacobian_function) :: f
    real, dimension(:,:), allocatable :: J
    real, dimension(:), intent(in) :: x 
    real, dimension(:), allocatable :: G
    ! real, dimension(:,:), intent(in) :: x
    ! real, dimension(:,:), allocatable :: J
    integer, dimension(:), allocatable :: o
    integer :: n, er, k
    real, dimension(:), allocatable :: delta_G, R
    
    n = size(x)
    
    allocate(o(n))
    allocate(delta_G(n))

    G = x
    
    do k=1, 10
        ! R = f(x)
        J = jacobian(f, x)
        R = -1*f(x)
        
        call ludecomp(J, n, 1e-5, o, er)
        
        if (er == 0) then
            call substitute(J, o, n, R, delta_G)
        else 
            write(*,*) "Unable to compute. Try lowering the tolerance"
        end if 

        G = G + delta_G
    end do
end function multivariable_newtons_method

end module my_functions



program Main
use my_functions
implicit none

!Jacobian Function example using the test function and 3x1 array
!Prints out the Jacobian to the console
real, dimension(:,:), allocatable :: J
real, dimension(:), allocatable :: hello, x
integer :: i, n

allocate(x(3))
x = [1, 2, 3]
! J = jacobian(test, x)
hello = multivariable_newtons_method(test, x)


do i=1, 3
    write(*,*) hello(i)
end do



end program Main
