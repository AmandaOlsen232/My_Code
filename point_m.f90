module point_m
    implicit none
    
    type point
        
        real x, y, z
        !~ real y
        !~ real z
        
    contains
        
        procedure :: init => point_init
        procedure :: mywrite => point_write
        
    end type point
    
    interface point
        procedure point_constructor
    end interface point
    
    interface operator (+)
        procedure point_add
        !~ procedure point_add_vector
    end interface operator (+)
    
    interface operator (-)
        procedure point_sub
        procedure point_negate
    end interface operator (-)
    
    type pointptr
        type(point), pointer :: p
    end type pointptr
    
contains
    
    function point_constructor(x, y, z) result(p)
        real, intent(in) :: x, y, z
        type(point) :: p
        p%x = x
        p%y = y
        p%z = z
    end function point_constructor
    
    subroutine point_write(this)
        class(point), intent(in) :: this
        write(*,*) this%x, this%y, this%z
    end subroutine point_write
    
    subroutine point_init(this, x, y, z)
        class(point), intent(inout) :: this
        real, intent(in) :: x, y, z
        this%x = x
        this%y = y
        this%z = z
    end subroutine point_init
    
    function point_add(p1,p2)
        type(point), intent(in) :: p1, p2
        type(point) :: point_add
        point_add%x = p1%x + p2%x
        point_add%y = p1%y + p2%y
        point_add%z = p1%z + p2%z
    end function point_add
    
    !~ function point_add_vector(p1,v2)
        !~ type(point), intent(in) :: p1
        !~ type(vector), intent(in) :: v2
        !~ type(point) :: point_add_vector
        !~ point_add_vector%x = p1%x + v2%x
        !~ point_add_vector%y = p1%y + v2%y
        !~ point_add_vector%z = p1%z + v2%z
    !~ end function point_add_vector
    
    function point_sub(p1,p2)
        type(point), intent(in) :: p1, p2
        type(point) :: point_sub
        point_sub%x = p1%x - p2%x
        point_sub%y = p1%y - p2%y
        point_sub%z = p1%z - p2%z
    end function point_sub
    
    function point_negate(p)
        type(point), intent(in) :: p
        type(point) :: point_negate
        point_negate%x = -p%x
        point_negate%y = -p%y
        point_negate%z = -p%z
    end function point_negate

end module point_m
