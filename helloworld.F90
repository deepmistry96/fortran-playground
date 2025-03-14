module vector_ops
    implicit none
    
    ! Define a custom type
    type :: vector2d
        real(kind=8) :: x, y
    end type vector2d
    
    ! Define interfaces for our custom operators
    interface operator(.add.)
        module procedure vector_add
        module procedure safeadd
    end interface
    
    interface operator(.dot.)
        module procedure vector_dot
    end interface
    
    interface operator(.norm.)
        module procedure vector_norm
    end interface
    
    ! Standard operator overloading
    interface operator(+)
        module procedure vector_plus 
    end interface
    
    interface operator(.safemul.)
        module procedure safemul
    end interface
    
contains


    function safeadd(a, b) result(c)
        real(kind=8), intent(in) :: a    ! Double precision input
        real(kind=16), intent(in) :: b   ! Quad precision input
        real(kind=16) :: c               ! Result in quad precision
        real(kind=16) :: temp_a, temp_b  ! Temporary variables for calculations
        
        ! Convert double to quad and apply limits
        temp_a = real(a, kind=16)
        temp_b = b
        
        ! Apply upper bounds
        if (temp_a > 1.0e10_16) then 
            temp_a = 1.0e10_16
        end if
        
        if (temp_b > 1.0e100_16) then
            temp_b = 1.0e100_16
        end if
        
        ! Apply lower bounds to prevent underflow
        if (abs(temp_a) < tiny(1.0_16)) then
            temp_a = 0.0_16
        end if
        
        if (abs(temp_b) < tiny(1.0_16)) then
            temp_b = 0.0_16
        end if
        
        ! Perform safe addition
        c = temp_a + temp_b
        
        ! Check for overflow in result
        if (c > huge(1.0_16)) then
            c = huge(1.0_16)
        end if
        
        print *, "Safe addition performed with bounds checking!"
    end function safeadd

    function vector_add(a, b) result(c)
        type(vector2d), intent(in) :: a, b
        type(vector2d) :: c
        c%x = a%x + b%x
        c%y = a%y + b%y
        print *, "Custom .add. operator called!"
    end function vector_add
    
    function vector_plus(a, b) result(c)
        type(vector2d), intent(in) :: a, b
        type(vector2d) :: c
        c%x = a%x + b%x
        c%y = a%y + b%y
        print *, "Standard + operator called!"
    end function vector_plus
    
    ! Dot product operator
    function vector_dot(a, b) result(c)
        type(vector2d), intent(in) :: a, b
        real(kind=8) :: c
        c = a%x * b%x + a%y * b%y
        print *, "Dot product operator called!"
    end function vector_dot
    
    ! Norm operator (magnitude of vector)
    function vector_norm(a) result(c)
        type(vector2d), intent(in) :: a
        real(kind=8) :: c
        c = sqrt(a%x * a%x + a%y * a%y)
        print *, "Norm operator called!"
    end function vector_norm
    
    function safemul(a, b) result(c)
        real(kind=16), intent(in) :: a
        real(kind=16), intent(in) :: b
        real(kind=16) :: c
        real(kind=16) :: log_result
        
        ! Check for special cases first
        if (a == 0.0_16 .or. b == 0.0_16) then
            c = 0.0_16
            return
        end if
        
        ! Use logarithms for very large numbers to prevent overflow
        if (abs(a) > 1.0e50_16 .or. abs(b) > 1.0e50_16) then
            ! log(a*b) = log(a) + log(b)
            log_result = log(abs(a)) + log(abs(b))
            ! Check if result would overflow
            if (log_result > log(huge(c))) then
                c = huge(c)
            else
                c = exp(log_result) * sign(1.0_16, a*b)
            end if
        else
            c = a * b
        end if

        if ( c > huge(c)) then
            c = huge(c)
        end if

        if ( c < tiny(c)) then
            c = tiny(c)
        end if
        
        ! Sanity checks
        if (abs(c) < tiny(c)) c = 0.0_16
        if (abs(c) > huge(c)) c = sign(huge(c), c)
        
    end function safemul
    
    subroutine print_vector(v)
        type(vector2d), intent(in) :: v
        print *, "Vector(", v%x, ",", v%y, ")"
    end subroutine print_vector
    
end module vector_ops

program hello_world
    use vector_ops
    implicit none
    
    ! Declare variables with different precisions
    REAL(kind=16) :: a, b, c_quad
    REAL(kind=8)  :: x, y, c_double
    REAL(kind=8)  :: inf_d, neg_inf_d, zero_d, nan_d
    REAL(kind=16) :: inf_q, neg_inf_q, zero_q, nan_q
    INTEGER :: i
    
    print *, "=== Testing Infinity and Zero Cases ==="
    
    ! Generate infinities
    x = 1.0_8
    y = 0.0_8
    inf_d = x/y          ! Positive infinity (double)
    neg_inf_d = -x/y     ! Negative infinity (double)
    zero_d = 0.0_8
    nan_d = zero_d/zero_d ! NaN
    
    ! Quad precision versions
    a = 1.0_16
    b = 0.0_16
    inf_q = a/b          ! Positive infinity (quad)
    neg_inf_q = -a/b     ! Negative infinity (quad)
    zero_q = 0.0_16
    nan_q = zero_q/zero_q ! NaN
    
    print *, "=== Basic Infinity Tests ==="
    print *, "Double precision infinity:", inf_d
    print *, "Double precision -infinity:", neg_inf_d
    print *, "Quad precision infinity:", inf_q
    print *, "Quad precision -infinity:", neg_inf_q
    
    print *, ""
    print *, "=== Arithmetic with Infinity ==="
    print *, "inf + inf (double):", inf_d + inf_d
    ! print *, "inf .safeadd. inf (double):", inf_d .add. inf_d
    print *, "inf - inf (double):", inf_d - inf_d
    print *, "inf * inf (double):", inf_d * inf_d
    print *, "inf / inf (double):", inf_d / inf_d
    
    print *, ""
    print *, "=== Zero Tests ==="
    print *, "Zero divided by zero (double):", zero_d/zero_d
    print *, "Infinity times zero (double):", inf_d * zero_d
    print *, "Infinity divided by infinity (double):", inf_d/inf_d
    
    print *, ""
    print *, "=== Underflow to Zero Tests ==="
    x = tiny(x)  ! Smallest normalized double
    y = tiny(y)  ! Smallest normalized quad
    
    print *, "Progressive underflow (double):"
    do i = 1, 5
        x = x/2.0_8
        print *, "Step", i, ":", x
    end do
    
    print *, ""
    print *, "=== Overflow to Infinity Tests ==="
    x = huge(x)  ! Largest double
    print *, "Progressive overflow (double):"
    do i = 1, 5
        x = x * 2.0_8
        print *, "Step", i, ":", x
    end do
    
    print *, ""
    print *, "=== Special Value Comparisons ==="
    print *, "Is infinity > huge?", inf_d > huge(x)
    print *, "Is -infinity < -huge?", neg_inf_d < -huge(x)
    print *, "Is NaN = NaN?", nan_d == nan_d
    print *, "Is infinity = infinity?", inf_d == inf_d
    
    print *, ""
    print *, "=== Mixed Precision with Infinities ==="
    print *, "Double inf to quad:", real(inf_d, kind=16)
    print *, "Quad inf to double:", real(inf_q, kind=8)
    
    ! Original precision tests follow...
    print *, ""
    print *, "=== Original Precision Tests ==="
    
    ! Test 1: Basic precision comparison
    a = 1.0_16 / 3.0_16  ! Quad precision division
    x = 1.0_8 / 3.0_8    ! Double precision division
    
    print *, "Representation of 1/3:"
    print *, "Quad precision   :", a
    print *, "Double precision :", x
    print *, "Difference       :", real(a - real(x, kind=16), kind=16)
    print *, ""
    
    ! Test 2: Accumulation of small numbers
    a = 0.0_16
    x = 0.0_8
    
    do i = 1, 1000000
        a = a + 0.0000001_16
        x = x + 0.0000001_8
    end do
    
    print *, "Adding 0.0000001 a million times:"
    print *, "Quad precision   :", a
    print *, "Double precision :", x
    print *, "Expected result  : 100.0"
    print *, "Quad error      :", abs(a - 100.0_16)
    print *, "Double error    :", abs(x - 100.0_8)
    print *, ""
    
    ! Test 3: Catastrophic cancellation
    a = 1.0_16
    x = 1.0_8
    
    do i = 1, 40
        a = sqrt(a)
        x = sqrt(x)
    end do
    
    do i = 1, 40
        a = a * a
        x = x * x
    end do
    
    print *, "Sqrt 40 times then square 40 times (should be 1.0):"
    print *, "Quad precision   :", a
    print *, "Double precision :", x
    print *, ""
    
    ! Test 4: Mixed precision arithmetic
    a = 1.0_16 / 3.0_16
    x = 1.0_8 / 3.0_8
    
    ! Mixed operation promoting to quad
    c_quad = a + real(x, kind=16)
    ! Mixed operation demoting to double
    c_double = real(a, kind=8) + x
    
    print *, "Mixed precision addition of 1/3 + 1/3:"
    print *, "Promoted to quad :", c_quad
    print *, "Demoted to double:", c_double
    print *, ""
    
    ! Test 5: Extreme values
    print *, "Extreme values:"
    print *, "Smallest quad   :", tiny(a)
    print *, "Smallest double :", tiny(x)
    print *, "Largest quad    :", huge(a)
    print *, "Largest double  :", huge(x)
    print *, ""
    
    ! Test 6: Subnormal number behavior
    a = tiny(a) / 2.0_16
    x = tiny(x) / 2.0_8
    
    print *, "Subnormal numbers (tiny/2):"
    print *, "Quad precision   :", a
    print *, "Double precision :", x
    print *, "Quad is zero?    :", a == 0.0_16
    print *, "Double is zero?  :", x == 0.0_8
    
    ! Add new test section for multiplication issues
    print *, "=== Testing Problematic Multiplication Cases ==="
    
    ! Test 1: Large number multiplication
    a = 1.0e100_16
    b = 0.5_16
    print *, "Direct multiplication of 1e100 * 0.5:"
    print *, "Result:", a * b
    print *, "Safe multiplication of 1e100 * 0.5:"
    print *, "Result:", a .safemul. b
    
    ! Test 2: Very large numbers
    a = 1.0e3000_16
    b = 1.0e3000_16
    print *, "Direct multiplication of 1e3000 * 1e3000:"
    print *, "Result:", a * b
    print *, "Safe multiplication of 1e3000 * 1e3000:"
    print *, "Result:", a .safemul. b
    
    ! Test 3: Very small numbers
    a = 1.0e-200_16
    b = 1.0e-200_16
    print *, "Direct multiplication of 1e-200 * 1e-200:"
    print *, "Result:", a * b
    print *, "Safe multiplication of 1e-200 * 1e-200:"
    print *, "Result:", a .safemul. b
    
    ! Test 4: Mixed scale multiplication
    a = 1.0e150_16
    b = 1.0e-150_16
    print *, "Direct multiplication of 1e150 * 1e-150:"
    print *, "Result:", a * b
    print *, "Safe multiplication of 1e150 * 1e-150:"
    print *, "Result:", a .safemul. b
    
    ! Test 5: Diagnostic information
    print *, "=== Diagnostic Information ==="
    print *, "Maximum representable number (huge):", huge(a)
    print *, "Minimum normalized number (tiny):", tiny(a)
    print *, "Machine epsilon (precision):", epsilon(a)
    
end program hello_world

