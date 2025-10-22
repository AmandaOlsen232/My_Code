module my_functions
    ! use math_m
    use linalg_mod 

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
    y(1) = x(1) + x(2) + x(3)**3 -3
    y(2) = x(2)**2 + x(3) -5 
    y(3) = x(3)**3 + x(2)**2 + x(1)**3 +2
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
        h = 1.0e-3
    end if

    !get the dimension of input and output vectors, allocate jacobian size
    dim_x = size(x)
    dim_fx = size(f(x))
    
    allocate(J(dim_fx, dim_x))
    allocate(eps(dim_x))
    J = 0.0
    eps = 0.0

    !numerically solve for the jacobian
    in: do i=1, dim_x
        eps = 0.0
        eps(i) = h
        fp = f(x+eps)
        fm = f(x-eps)
        J(:,i) = (fp - fm)/(2.*h)
    end do in

    
end function jacobian

function multivariable_newtons_method(f, x, tol, max_it, relax_factor, verbose) result(G)
    !for a given function and initial input, solve for the x required to make the function equal zero
    implicit none 
    procedure(jacobian_function) :: f
    real, dimension(:), intent(in) :: x 
    real, intent(in), optional :: tol 
    real :: tolerance 
    integer, intent(in), optional :: max_it 
    integer :: m_it
    real, intent(in), optional :: relax_factor
    real :: relax 
    logical, intent(in), optional :: verbose
    logical :: verb
    real, dimension(:), allocatable :: G

    real, dimension(:,:), allocatable :: J, J_inv
    integer :: dim_x, dim_fx
    integer, dimension(:), allocatable :: o
    integer :: n, er
    real, dimension(:), allocatable :: delta_G, R, test_dG
    integer :: iterations 
    real :: max_r 

    !set defaults for the optional values
    if (present(tol)) then
        tolerance = tol 
    else
        tolerance = 1.0e-9
    end if

    if (present(max_it)) then
        m_it = max_it 
    else
        m_it = 50
    end if

    if (present(relax_factor)) then 
        relax = relax_factor
    else 
        relax = 1.0 
    end if 

    if (present(verbose)) then 
        verb = verbose
    else 
        verb = .false. 
    end if 

    !allocate memory for o and delta_G
    n = size(x)
    dim_x = size(x)
    dim_fx = size(f(x))

    allocate(o(n))
    allocate(delta_G(n))
    allocate(G(n))
    allocate(J(dim_fx, dim_x))
    allocate(J_inv(dim_fx, dim_x))
    allocate(R(dim_x))
    allocate(test_dG(dim_x))

    if (verb .eqv. .true.) then
        write(*,*) "Newton Solver: iteration   max_R"
        write(*,*) "___________________________________"
    end if
    
    !run the while loop
    G = x
    max_r = 1000.
    iterations = 1
    do while ((max_r > tolerance) .and. (iterations < m_it))
        ! write(*,*) "iteration: ", iterations
        J = jacobian(f, G)
        R = -1.*f(G)

        !! These lines use Zach's math_m module
        ! call ludecomp(J, n, 1.0e-5, o, er) !could raise up tolerance on this
        ! if (er == 0) then
        !     call substitute(J, o, n, R, delta_G)
        ! else 
        !     write(*,*) "Unable to compute LU Decomposition"
        ! end if 

        !this 2 lines use linalg.f90
        call matinv(dim_fx, J, J_inv)
        delta_G = matmul(J_inv, R)
 
        G = G + relax*delta_G
        
        R = abs(R)
        max_r = maxval(R)
        iterations = iterations + 1

        if (verb .eqv. .true.) then
            write(*,*) iterations, max_r
        end if
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
x = [1., 4., 3.]
! J = jacobian(test, x)
hello = multivariable_newtons_method(test, x, verbose=.true.)
x = test(hello)
write(*,*) hello
write(*,*) x
! do i=1, 3
!     write(*,*) J(i,:)
! end do



end program Main
