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

    interface operator(.safeexp.)
        module procedure safeexp
    end interface

    interface operator(.safediv.)
        module procedure safediv
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

    function safeexp(base, exponent) result(c)
        real(kind=16), intent(in) :: base
        real(kind=16), intent(in) :: exponent
        real(kind=16) :: c
        real(kind=16) :: log_base, total_log
        
        ! Handle special cases first
        if (base == 0.0_16) then
            if (exponent > 0.0_16) then
                c = 0.0_16
            else if (exponent == 0.0_16) then
                c = 1.0_16
            else
                c = huge(c)  ! Infinity for negative exponent of zero
            end if
            return
        end if
        
        if (exponent == 0.0_16) then
            c = 1.0_16
            return
        end if
        
        if (base == 1.0_16) then
            c = 1.0_16
            return
        end if
        
        ! For negative bases, only allow integer exponents
        if (base < 0.0_16) then
            if (abs(exponent - nint(exponent)) > tiny(1.0_16)) then
                c = 0.0_16  ! Return 0 for non-integer exponents with negative base
                return
            end if
        end if
        
        ! Use logarithms for very large or very small numbers
        if (abs(base) > 1.0_16) then
            ! For positive bases, use log method
            log_base = log(abs(base))
            total_log = log_base * exponent
            
            ! Check if result would overflow
            if (total_log > log(huge(c))) then
                c = huge(c)
            else if (total_log < log(tiny(c))) then
                c = tiny(c)
            else
                c = exp(total_log)
                ! Handle negative bases with odd exponents
                if (base < 0.0_16 .and. mod(nint(exponent), 2) == 1) then
                    c = -c
                end if
            end if
        else if (abs(base) < 1.0_16 .and. base /= 0.0_16) then
            ! Similar approach for small bases
            log_base = log(abs(base))
            total_log = log_base * exponent
            
            if (total_log < log(tiny(c))) then
                c = tiny(c)
            else if (total_log > log(huge(c))) then
                c = huge(c)
            else
                c = exp(total_log)
                ! Handle negative bases with odd exponents
                if (base < 0.0_16 .and. mod(nint(exponent), 2) == 1) then
                    c = -c
                end if
            end if
        else
            ! For bases close to 1, use direct computation
            c = base ** exponent
        end if
        
        ! Final sanity checks
        if (abs(c) < tiny(c)) c = 0.0_16
        if (abs(c) > huge(c)) c = sign(huge(c), c)
        
        print *, "Safe exponentiation performed with bounds checking!"
    end function safeexp

    function safediv(numerator, denominator) result(c)
        real(kind=16), intent(in) :: numerator
        real(kind=16), intent(in) :: denominator
        real(kind=16) :: c
        real(kind=16) :: log_num, log_den, log_result
        real(kind=16), parameter :: MIN_DENOMINATOR = tiny(1.0_16) * 1.0e10_16
        
        ! Handle division by zero or very small numbers
        if (abs(denominator) < MIN_DENOMINATOR) then
            if (denominator == 0.0_16) then
                if (numerator > 0.0_16) then
                    c = huge(c)
                else if (numerator < 0.0_16) then
                    c = -huge(c)
                else
                    c = 0.0_16  ! 0/0 = 0 by convention in our safe math
                end if
            else
                ! Very small denominator - clamp to avoid underflow
                c = sign(huge(c), numerator * denominator)
            end if
            return
        end if
        
        ! Handle zero numerator
        if (numerator == 0.0_16) then
            c = 0.0_16
            return
        end if
        
        ! For very large or very small numbers, use logarithmic division
        if (abs(numerator) > 1.0e50_16 .or. abs(denominator) > 1.0e50_16 .or. &
            abs(numerator) < 1.0e-50_16 .or. abs(denominator) < 1.0e-50_16) then
            
            log_num = log(abs(numerator))
            log_den = log(abs(denominator))
            log_result = log_num - log_den
            
            ! Check for overflow/underflow in result
            if (log_result > log(huge(c))) then
                c = sign(huge(c), numerator * denominator)
            else if (log_result < log(tiny(c))) then
                c = sign(tiny(c), numerator * denominator)
            else
                c = exp(log_result) * sign(1.0_16, numerator * denominator)
            end if
        else
            ! Normal range - direct division
            c = numerator / denominator
        end if
        
        ! Final sanity checks
        if (abs(c) < tiny(c)) c = 0.0_16
        if (abs(c) > huge(c)) c = sign(huge(c), c)
        
        print *, "Safe division performed with bounds checking!"
    end function safediv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    



    subroutine Basic_precision_comparison
        REAL(kind=16) :: a = 1.0_16 / 3.0_16  ! Quad precision division
        REAL(kind=8) :: x = 1.0_8 / 3.0_8    ! Double precision division
        
        print *, "Representation of 1/3:"
        print *, "Quad precision   :", a
        print *, "Double precision :", x
        print *, "Difference       :", real(a - real(x, kind=16), kind=16)
        print *, ""
    end subroutine Basic_precision_comparison
        
    subroutine Accumulation_of_small_numbers
        REAL(kind=16) :: a = 0.0_16
        REAL(kind=8) :: x = 0.0_8
        integer i 
        
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
    end subroutine Accumulation_of_small_numbers
        
    subroutine Catastrophic_cancellation
        REAL(kind=16) :: a = 1.0_16
        REAL(kind=8) :: x = 1.0_8
        integer i
        
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
    end subroutine Catastrophic_cancellation
    
    subroutine Mixed_precision_arithmetic
        REAL(kind=16) :: a = 1.0_16 / 3.0_16
        REAL(kind=8) :: x = 1.0_8 / 3.0_8
        REAL(kind=16) :: c_quad
        REAL(kind=8) :: c_double
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
    end subroutine Mixed_precision_arithmetic     
        
    subroutine Subnormal_number_behavior
        REAL(kind=16) :: a = tiny(a) / 2.0_16
        REAL(kind=8) :: x = tiny(x) / 2.0_8
        
        print *, "Subnormal numbers (tiny/2):"
        print *, "Quad precision   :", a
        print *, "Double precision :", x
        print *, "Quad is zero?    :", a == 0.0_16
        print *, "Double is zero?  :", x == 0.0_8
        print *, ""
    end subroutine Subnormal_number_behavior
        
    subroutine Problematic_multiplication_cases
        ! Add new test section for multiplication issues
        print *, "=== Testing Problematic Multiplication Cases ==="
    end subroutine Problematic_multiplication_cases
            
    subroutine Large_number_multiplication
        REAL(kind=16) :: a = 1.0e100_16
        real(kind=16) :: b = 0.5_16
        print *, "Direct multiplication of 1e100 * 0.5:"
        print *, "Result:", a * b
        print *, "Safe multiplication of 1e100 * 0.5:"
        print *, "Result:", a .safemul. b
    end subroutine Large_number_multiplication

    subroutine Very_large_numbers
        REAL(kind=16) :: a = 1.0e3000_16
        real(kind=16) :: b = 1.0e3000_16
        print *, "Direct multiplication of 1e3000 * 1e3000:"
        print *, "Result:", a * b
        print *, "Safe multiplication of 1e3000 * 1e3000:"
        print *, "Result:", a .safemul. b
    end subroutine Very_large_numbers
        
    subroutine Very_small_numbers
        REAL(kind=16) :: a = 1.0e-200_16
        real(kind=16) :: b = 1.0e-200_16
        print *, "Direct multiplication of 1e-200 * 1e-200:"
        print *, "Result:", a * b
        print *, "Safe multiplication of 1e-200 * 1e-200:"
        print *, "Result:", a .safemul. b
    end subroutine Very_small_numbers

    subroutine Safe_exponentiation_tests
        REAL(kind=16) :: base, exponent, result
        
        print *, "=== Testing Safe Exponentiation Cases ==="
        print *, ""
        
        ! Test 1: Basic positive cases
        base = 2.0_16
        exponent = 3.0_16
        result = base .safeexp. exponent
        print *, "Basic case (2^3):"
        print *, "Expected: 8.0"
        print *, "Result  :", result
        print *, ""
        
        ! Test 2: Zero base cases
        base = 0.0_16
        print *, "Zero base cases:"
        print *, "0^3  :", base .safeexp. 3.0_16
        print *, "0^0  :", base .safeexp. 0.0_16
        print *, "0^-2 :", base .safeexp. -2.0_16
        print *, ""
        
        ! Test 3: Very large exponents
        base = 1.1_16
        exponent = 1000.0_16
        print *, "Large exponent (1.1^1000):"
        print *, "Direct:", base ** exponent
        print *, "Safe  :", base .safeexp. exponent
        print *, ""
        
        ! Test 4: Negative base with integer and non-integer exponents
        base = -2.0_16
        print *, "Negative base cases:"
        print *, "(-2)^3    :", base .safeexp. 3.0_16
        print *, "(-2)^2    :", base .safeexp. 2.0_16
        print *, "(-2)^2.5  :", base .safeexp. 2.5_16  ! Should return 0 for non-integer exp
        print *, ""
        
        ! Test 5: Very small numbers
        base = 0.1_16
        exponent = 1000.0_16
        print *, "Very small result (0.1^1000):"
        print *, "Direct:", base ** exponent
        print *, "Safe  :", base .safeexp. exponent
        print *, ""
        
        ! Test 6: Overflow cases
        base = 2.0_16
        exponent = 1000.0_16
        print *, "Overflow case (2^1000):"
        print *, "Safe  :", base .safeexp. exponent
        print *, "Should be capped at huge()"
        print *, ""
        
        ! Test 7: Special values
        print *, "Special values:"
        print *, "1^anything:", 1.0_16 .safeexp. 1000.0_16
        print *, "anything^0:", 12345.0_16 .safeexp. 0.0_16
        print *, ""
        
    end subroutine Safe_exponentiation_tests
        
    subroutine Mixed_scale_multiplication
        REAL(kind=16) :: a = 1.0e150_16
        REAL(kind=16) :: b = 1.0e-150_16
        print *, "Direct multiplication of 1e150 * 1e-150:"
        print *, "Result:", a * b
        print *, "Safe multiplication of 1e150 * 1e-150:"
        print *, "Result:", a .safemul. b
    end subroutine Mixed_scale_multiplication
        
    subroutine Diagnostic_information
    REAL(kind=16) :: a = 0.0
        print *, "=== Diagnostic Information ==="
        print *, "Maximum representable number (huge):", huge(a)
        print *, "Minimum normalized number (tiny):", tiny(a)
        print *, "Machine epsilon (precision):", epsilon(a)
    end subroutine Diagnostic_information

    subroutine Safe_division_tests
        REAL(kind=16) :: num, den, result
        
        print *, "=== Testing Safe Division Cases ==="
        print *, ""
        
        ! Test 1: Basic division
        num = 10.0_16
        den = 2.0_16
        print *, "Basic division (10/2):"
        print *, "Direct:", num / den
        print *, "Safe  :", num .safediv. den
        print *, ""
        
        ! Test 2: Division by zero
        num = 1.0_16
        den = 0.0_16
        print *, "Division by zero cases:"
        print *, "1/0    :", num .safediv. den
        print *, "-1/0   :", (-num) .safediv. den
        print *, "0/0    :", 0.0_16 .safediv. den
        print *, ""
        
        ! Test 3: Very small denominators
        den = tiny(1.0_16) * 0.1_16
        print *, "Very small denominator:"
        print *, "1/tiny:", 1.0_16 .safediv. den
        print *, ""
        
        ! Test 4: Very large numbers
        num = 1.0e200_16
        den = 1.0e-200_16
        print *, "Large number division:"
        print *, "Direct:", num / den
        print *, "Safe  :", num .safediv. den
        print *, ""
        
        ! Test 5: Very small numbers
        num = 1.0e-200_16
        den = 1.0e200_16
        print *, "Small number division:"
        print *, "Direct:", num / den
        print *, "Safe  :", num .safediv. den
        print *, ""
        
        ! Test 6: Division resulting in 1
        num = 1.0e100_16
        den = 1.0e100_16
        print *, "Large/Large = 1 case:"
        print *, "Direct:", num / den
        print *, "Safe  :", num .safediv. den
        print *, ""
        
        ! Test 7: Negative number cases
        print *, "Negative number cases:"
        print *, "-10/2  :", (-10.0_16) .safediv. 2.0_16
        print *, "10/-2  :", 10.0_16 .safediv. (-2.0_16)
        print *, "-10/-2 :", (-10.0_16) .safediv. (-2.0_16)
        print *, ""
        
    end subroutine Safe_division_tests

    subroutine print_vector(v)
        type(vector2d), intent(in) :: v
        print *, "Vector(", v%x, ",", v%y, ")"
    end subroutine print_vector

    subroutine check_precision
        INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33)
        REAL(KIND=qp) :: x
        PRINT *, "Quad Precision KIND: ", qp
    end subroutine check_precision

    subroutine misc_tests

        ! Declare variables with different precisions
        REAL(kind=16) :: a, b, c_quad
        REAL(kind=8)  :: x, y, c_double
        REAL(kind=8)  :: inf_d, neg_inf_d, zero_d, nan_d
        REAL(kind=16) :: inf_q, neg_inf_q, zero_q, nan_q
        INTEGER :: i



        ! Generate infinities
        x = 1.0_8
        y = 0.0_8
        inf_d = x/y          ! Positive infinity (double)
        neg_inf_d = -x/y     ! Negative infinity (double)
        zero_d = 0.0_8
        nan_d = zero_d/zero_d ! NaN

        print *, "=== Testing Infinity and Zero Cases ==="

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
    
    end subroutine misc_tests



    
end module vector_ops

program hello_world
    use vector_ops
    implicit none
    
    call misc_tests

    call  Accumulation_of_small_numbers
    call  Basic_precision_comparison
    call  Catastrophic_cancellation
    call  Diagnostic_information
    call  Large_number_multiplication
    call  Mixed_precision_arithmetic
    call  Mixed_scale_multiplication
    call  Problematic_multiplication_cases
    call  Subnormal_number_behavior
    call  Very_large_numbers
    call  Very_small_numbers
    call  Safe_exponentiation_tests
    call  Safe_division_tests
    call  check_precision
    call  misc_tests
    ! call  print_vector(v)


    

    
        
end program hello_world

