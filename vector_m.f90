module vector_m
    use point_m
    implicit none
    
    type vector
        
        real x, y, z
        
    contains
        
        procedure :: mag => vector_mag
        procedure :: init => vector_init
        procedure :: mywrite => vector_write
        procedure :: convert2point => vector_convert2point
        
    end type vector
    
    interface vector
        procedure vector_constructor
    end interface vector
    
    interface operator (+)
        procedure vector_add_vector
        procedure point_add_vector
        procedure vector_add_point
    end interface operator (+)
    
    interface operator (-)
        procedure vector_sub_vector
        procedure point_sub_vector
        procedure vector_sub_point
        procedure vector_negate
    end interface operator (-)
    
    interface operator (*)
        procedure scalar_times_vector
        procedure vector_times_scalar
    end interface operator (*)
    
    interface operator (/)
        procedure vector_over_scalar
    end interface operator (/)
    
    interface operator (.vec.)
        procedure create_vector
    end interface operator (.vec.)
    
    interface operator (.dot.)
        procedure dot_prod
    end interface operator (.dot.)
    
    interface operator (.cross.)
        procedure cross_product
    end interface operator (.cross.)
    
    type vectorptr
        type(vector), pointer :: v
    end type vectorptr
    
contains
    
    function vector_constructor(x, y, z) result(v)
        real, intent(in) :: x, y, z
        type(vector) :: v
        v%x = x
        v%y = y
        v%z = z
    end function vector_constructor
    
    subroutine vector_write(this)
        class(vector), intent(in) :: this
        write(*,*) this%x, this%y, this%z
    end subroutine vector_write
    
    subroutine vector_init(this, x, y, z)
        class(vector), intent(inout) :: this
        real :: x, y, z
        this%x = x
        this%y = y
        this%z = z
    end subroutine vector_init
    
    function vector_mag(this)
        class(vector), intent(in) :: this
        real :: vector_mag
        vector_mag = sqrt(this%x ** 2.0 + this%y ** 2.0 + this%z ** 2.0)
    end function vector_mag
    
    function vector_add_vector(v1,v2)
        type(vector), intent(in) :: v1, v2
        type(vector) vector_add_vector
        vector_add_vector%x = v1%x + v2%x
        vector_add_vector%y = v1%y + v2%y
        vector_add_vector%z = v1%z + v2%z
    end function vector_add_vector
    
    function vector_add_point(v1,p2)
        type(vector), intent(in) :: v1
        type(point), intent(in) :: p2
        type(point) vector_add_point
        vector_add_point%x = v1%x + p2%x
        vector_add_point%y = v1%y + p2%y
        vector_add_point%z = v1%z + p2%z
    end function vector_add_point
    
    function point_add_vector(p1,v2)
        type(vector), intent(in) :: v2
        type(point), intent(in) :: p1
        type(point) point_add_vector
        point_add_vector = v2 + p1
    end function point_add_vector
    
    function vector_sub_vector(v1,v2)
        type(vector), intent(in) :: v1, v2
        type(vector) vector_sub_vector
        vector_sub_vector%x = v1%x - v2%x
        vector_sub_vector%y = v1%y - v2%y
        vector_sub_vector%z = v1%z - v2%z
    end function vector_sub_vector
    
    function vector_sub_point(v1,p2)
        type(vector), intent(in) :: v1
        type(point), intent(in) :: p2
        type(vector) vector_sub_point
        vector_sub_point%x = v1%x - p2%x
        vector_sub_point%y = v1%y - p2%y
        vector_sub_point%z = v1%z - p2%z
    end function vector_sub_point
    
    function point_sub_vector(p1,v2)
        type(vector), intent(in) :: v2
        type(point), intent(in) :: p1
        type(vector) point_sub_vector
        point_sub_vector%x = p1%x - v2%x
        point_sub_vector%y = p1%y - v2%y
        point_sub_vector%z = p1%z - v2%z
    end function point_sub_vector
    
    function vector_negate(v)
        type(vector), intent(in) :: v
        type(vector) :: vector_negate
        vector_negate%x = -v%x
        vector_negate%y = -v%y
        vector_negate%z = -v%z
    end function vector_negate
    
    function scalar_times_vector(s,v)
        real, intent(in) :: s
        type(vector), intent(in) :: v
        type(vector) :: scalar_times_vector
        scalar_times_vector%x = s * v%x
        scalar_times_vector%y = s * v%y
        scalar_times_vector%z = s * v%z
    end function scalar_times_vector
    
    function vector_times_scalar(v,s)
        real, intent(in) :: s
        type(vector), intent(in) :: v
        type(vector) :: vector_times_scalar
        vector_times_scalar = s * v
    end function vector_times_scalar
    
    function vector_over_scalar(v,s)
        real, intent(in) :: s
        type(vector), intent(in) :: v
        type(vector) :: vector_over_scalar
        real :: sinv
        sinv = 1.0 / s
        vector_over_scalar = sinv * v
    end function vector_over_scalar
    
    function create_vector(p1,p2)
        type(point), intent(in) :: p1, p2
        type(vector) :: create_vector
        create_vector%x = p2%x - p1%x
        create_vector%y = p2%y - p1%y
        create_vector%z = p2%z - p1%z
    end function create_vector
    
    function dot_prod(v1,v2)
        type(vector), intent(in) :: v1, v2
        real :: dot_prod
        dot_prod = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
    end function dot_prod
    
    function cross_product(v1,v2)
        type(vector), intent(in) :: v1, v2
        type(vector) :: cross_product
        cross_product%x = v1%y * v2%z - v1%z * v2%y
        cross_product%y = v1%z * v2%x - v1%x * v2%z
        cross_product%z = v1%x * v2%y - v1%y * v2%x
    end function cross_product
    
    function vector_convert2point(t)
        class(vector), intent(in) :: t
        type(point) :: vector_convert2point
        call vector_convert2point%init(t%x, t%y, t%z)
    end function vector_convert2point
    
end module vector_m
