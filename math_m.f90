module math_m
    use point_m
    use vector_m
    implicit none
    
    real, parameter :: pi = 3.1415926535897932384626, &
                           d2r = pi / 180.0, &
                           r2d = 180.0 / pi
    
    
    real, parameter :: Id_matrix(3,3) = reshape([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], shape=[3,3])
    
    type rotation_matrix
        
        real, dimension(3,3) :: mat
        
    contains
        
        procedure :: init => rotation_matrix_init
        procedure :: T => rotation_matrix_transpose
        
    end type rotation_matrix
    
    interface operator (*)
        procedure matrix_times_matrix
        procedure matrix_times_vector
    end interface operator (*)
    
    interface rotate
        procedure :: rotate_point
        procedure :: rotate_vector
    end interface rotate
    
contains
    
    subroutine rotation_matrix_init(t, x, y, z, order)
        class(rotation_matrix), intent(inout) :: t
        real, intent(in) :: x, y, z
        character(len=*), intent(in) :: order
        integer :: i
        
        t%mat = Id_matrix
        
        do i=1, len(order)
            if (order(i:i) == 'x') then
                t%mat = matmul(rot_x(x), t%mat)
            else if (order(i:i) == 'y') then
                t%mat = matmul(rot_y(y), t%mat)
            else if (order(i:i) == 'z') then
                t%mat = matmul(rot_z(z), t%mat)
            else
                write(*,*) 'Error! Unknown rotation direction: ', order(i:i)
                stop
            end if
        end do
        
    end subroutine rotation_matrix_init
    
    function   rotation_matrix_transpose(t) result(mt)
        class(rotation_matrix), intent(in) :: t
        type(rotation_matrix) :: mt
        integer :: row, col
        do row=1, 3
            do col=1, 3
                mt%mat(col, row) = t%mat(row, col)
            end do
        end do
    end function rotation_matrix_transpose
    
    function   matrix_times_vector(m, v)
        type(rotation_matrix), intent(in) :: m
        type(vector), intent(in) :: v
        real, dimension(3) :: b
        type(vector) :: matrix_times_vector
        b = [v%x, v%y, v%z]
        b = matmul(m%mat, b)
        call matrix_times_vector%init(b(1), b(2), b(3))
    end function matrix_times_vector
    
    function   matrix_times_matrix(m, n)
        type(rotation_matrix), intent(in) :: m, n
        type(rotation_matrix) :: matrix_times_matrix
        matrix_times_matrix%mat = matmul(m%mat, n%mat)
    end function matrix_times_matrix
    
    !===================================================================
    
    function   rot_x(a) result(m)
        real, intent(in) :: a
        real, dimension(3,3) :: m
        real :: s, c
        s = sin(a)
        c = cos(a)
        m = reshape([1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c], shape=[3,3])
    end function rot_x
    
    function   rot_y(a) result(m)
        real, intent(in) :: a
        real, dimension(3,3) :: m
        real :: s, c
        s = sin(a)
        c = cos(a)
        m = reshape([c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c], shape=[3,3])
    end function rot_y
    
    function   rot_z(a) result(m)
        real, intent(in) :: a
        real, dimension(3,3) :: m
        real :: s, c
        s = sin(a)
        c = cos(a)
        m = reshape([c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0], shape=[3,3])
    end function rot_z
    
    !C
    !C--------------------------------------------------------------------C
    !C The following subroutine computes the LU decomposition for a       C
    !C diagonally dominant matrix (no pivoting is done).                  C
    !C   Inputs:  n = number of equations/unknowns                        C
    !C            a = nxn coefficient matrix                              C
    !C   Outputs: a = nxn matrix containing the LU matrices               C
    !C                                                                    C
    !C Deryl Snyder, 10-16-98                                             C
    !C modified: Zach Montgomery 03-03-2025
    !C--------------------------------------------------------------------C
    subroutine math_snyder_ludcmp(a)
        real, intent(inout), dimension(:,:) :: a
        integer :: n, i, j, k
        real :: z
        n = size(a, 1)
        do k=1, n-1
            do i=k+1, n
                z = a(i,k) / a(k,k)                     !compute gauss factor 
                a(i,k) = z                            !store gauss factor in matrix
                do j=k+1, n
                    a(i,j) = a(i,j) - z*a(k,j)            !apply row operation
                end do
            end do
        end do
        return
    end subroutine math_snyder_ludcmp
    !C--------------------------------------------------------------------C
    !C The following subroutine solves for the unknowns (x) given the LU  C
    !C matrix and the right hand side.                                    C
    !C   Inputs:  n = number of equations/unknowns                        C
    !C            a = nxn matrix containing the L and U values            C
    !C            b = n vector containing right hand side values          C
    !C   Outputs: x = n vector containing solution                        C
    !C                                                                    C
    !C Deryl Snyder, 10-16-98                                             C
    !C modified: Zach Montgomery, 03-03-2025
    !C--------------------------------------------------------------------C
    subroutine math_snyder_lusolv(a,b,x)
        real, intent(in), dimension(:,:) :: a
        real, intent(in), dimension(:) :: b
        real, intent(out), dimension(:) :: x
        integer :: n, i, j, k
        n = size(b)
        !~ real :: a(n,n),b(n),x(n)
        do i=1, n
            x(i) = b(i)
        end do
        do k=1, n-1                                 !do forward substitution
            do i=k+1, n
                x(i) = x(i) - a(i,k)*x(k)
            end do
        end do
        do i=n, 1, -1                                !do back subsitution
            do j=i+1, n
                x(i) = x(i) - a(i,j)*x(j)
            end do
            x(i) = x(i) / a(i,i)
        end do
        !~ return
    end subroutine math_snyder_lusolv
    
    
!***************************************************************************
subroutine ludecomp(a, n, tol, o, er)
integer, intent(in) :: n
integer, intent(inout), dimension(n) :: o
real, intent(in) :: tol
real, intent(inout), dimension(n,n) :: a
real, dimension(n) :: s
integer, intent(inout) :: er
er = 0
call decompose(a, n, tol, o, s, er)
end subroutine ludecomp
!***************************************************************************
subroutine decompose(a, n, tol, o, s, er)
integer, intent(in) :: n
integer, intent(inout) :: er
integer, intent(inout), dimension(n) :: o
integer :: k, i, j
real, intent(in) :: tol
real, intent(inout), dimension(n) :: s
real, intent(inout), dimension(n,n) :: a
real :: factor
do i = 1, n
    o(i) = i
    s(i) = abs(a(i,1))
    do j = 2, n
        if (abs(a(i,j)) > s(i)) s(i) = abs(a(i,j))
    end do
end do
do k = 1, n-1
    call pivot(a, o, s, n, k)
    if (abs(a(o(k),k)/s(o(k))) < tol) then
        er = -1
        write(*,*) a(o(k),k)/s(o(k))
        exit
    end if
    do i = k+1, n
        factor = a(o(i),k)/a(o(k),k)
        a(o(i),k) = factor
        do j = k+1,n
            a(o(i),j) = a(o(i),j) - factor*a(o(k),j)
        end do
    end do
end do
if (abs(a(o(k),k)/s(o(k)))<tol) then
    er = -1
    write(*,*) a(o(k),k)/s(o(k))
end if
end subroutine decompose
!***************************************************************************
subroutine pivot(a, o, s, n, k)
integer, intent(in) :: n, k
integer, intent(inout), dimension(n) :: o
real, intent(inout), dimension(n) :: s
real, intent(inout), dimension(n,n) :: a
real :: big, dummy
integer :: p, dumb, i
p = k
big = abs(a(o(k),k)/s(o(k)))
do i = k+1,n
    dummy = abs(a(o(i),k)/s(o(i)))
    if (dummy > big) then
        big = dummy
        p = i
    end if
end do
dumb = o(p)
o(p) = o(k)
o(k) = dumb
end subroutine pivot
!***************************************************************************
subroutine substitute(a, o, n, b, gamma)
integer, intent(in) :: n
integer, intent(inout), dimension(n) :: o
real, intent(inout), dimension(n) :: b, gamma
real, intent(inout), dimension(n,n) :: a
real :: total
integer :: i, j
do i = 2, n
    total = b(o(i))
    do j = 1, i-1
        total = total - a(o(i),j)*b(o(j))
    end do
    b(o(i)) = total
end do
gamma(n) = b(o(n))/a(o(n),n)
do i = n-1,1,-1
    total = 0.0
    do j = i+1,n
        total = total + a(o(i),j)*gamma(j)
    end do
    gamma(i) = (b(o(i)) - total)/a(o(i),i)
end do
end subroutine substitute
!***************************************************************************
    
    function m_norm(v) result(x)
        real, dimension(:) :: v
        real :: x
        integer :: i
        x = 0.0
        do i=1, size(v)
            x = x + v(i) ** 2.0
        end do
        x = x ** 0.50
    end function m_norm
    
    subroutine rotate_point(rot, p, z)
        type(rotation_matrix), intent(in) :: rot
        type(point), intent(inout) :: p
        type(point), intent(in), optional :: z
        type(point) :: zero
        type(vector) :: r
        call zero%init(0.0, 0.0, 0.0)
        r = zero .vec. p
        r = rot * r
        p = r%convert2point()
        if (present(z)) p = p + z
    end subroutine rotate_point
    
    subroutine rotate_vector(rot, v)
        type(rotation_matrix), intent(in) :: rot
        type(vector), intent(inout) :: v
        v = rot * v
    end subroutine rotate_vector
    
    
    
end module math_m
