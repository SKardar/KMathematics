import string
from fractions import Fraction
from math import ceil, atan, sin, cos, pi, e

def trunc(number : float, digits : int = 0) -> float | str:
    """truncate will cut a float number at a given decimal without rounding. (round toward zero)

    Args:
        number (float): the number to be cut
        digits (int): number of decimals, Defaults to 0 decimals

    Returns:
        float: truncated version of number with n decimals
    """
    try:
        digits = int(float(digits))
        pow10 = 10 ** digits
        return number * pow10 // 1 / pow10
    except ValueError:
        return f"Couldn't convert the second argument '{digits}' to an Integer number"


def sieve(n : int | float, p : int = 0) -> list[int] | str:
    """Sieve of Eratosthenes:
        It is an ancient algorithm for finding all prime numbers up to any given limit. 
        It does so by iteratively marking as composite (i.e., not prime) the multiples of each prime,
        starting with the first prime number, 2.

    Args:
        n (int | float): Any number greater than or equal to 2.
        Note: input only 1 argument (n) to output the list of primes up to 'n'.
        
        p (Optional) : should be a prime number.
        Note: if you input a prime number for 'p', it will return the list of composite numbers eliminated in the stage of removing the multiples of 'p'.

    Returns:
        list:
            for 1 argument: sieve(n) -> A list of all prime numbers less than or equal to 'n'.
            for 2 arguments: sieve(n,p) -> A list of all composite numbers eliminated in removing the multiples of 'p'.
    """
    try:
        assert n >= 2
        num_list = list(range(2,int(n+1)))
        composites_p = []
        for i in num_list:
            active_prime = i
            if i**2 in num_list:
                for k in list(filter(lambda m: m>=i**2, num_list)):
                    if k % active_prime == 0:
                        if p == active_prime:
                            composites_p.append(k)
                        num_list.remove(k)
        if p:
            return composites_p
        else:
            return num_list
    except AssertionError:
        return "AssertionError: sieve (Sieve of Eratosthenes) function needs a positive number greater than or equal to 2 as input"


def prime_factorise(n: int) -> list[int] | str:
    """Prime Factorise (Trial Division): Gives a list of all Prime Factors of 'n'.

    Args:
        n (int): Number to Prime Factorise.

    Returns:
        list[int]: List of all Prime Factors of the given number.
    """
    try:
        assert isinstance(n,int)
        # Handling zero
        if n == 0:
            return []
        
        prime_factors = []
        # Handling negative integers
        if n < 0:
            prime_factors.append(-1)
            n = abs(n)
            
        # Factor out all occurenses of 2
        while n % 2 == 0:
            prime_factors.append(2)
            n //= 2
            
        # Factoring out other odd primes
        f = 3
        while f * f <= n:
            if n % f == 0:
                prime_factors.append(f)
                n //= f
            else:
                f += 2
        
        # Handling last (prime)^2 number:
        if n != 1: prime_factors.append(n)
        
        return prime_factors
    
    except AssertionError:
        return "AssertionError: Please input an integer number to prime_factorise() function"
        

def factors(n: int) -> list[int] | str:
    """Gives a list of all Factors of 'n'.

    Args:
        n (int): Number to Factorise.

    Returns:
        list[int]: List of all Factors of the given number.
    """
    try:
        assert isinstance(n,int)
        neg = 0
        if n == 0:
            return [0] 
        else:
            if n < 0:
                neg = 1
                n = abs(n)
            factors = [1,]
            for i in range(2,n+1):
                if n % i == 0:
                    factors.append(i)
            if neg == 1:
                for l in factors.copy():
                    factors.append(-l)
            return factors
    except AssertionError:
        return "AssertionError: Please input an integer number to factors() function"

perfects = [6, 28, 496, 8128, 33550336, 8589869056, 137438691328,
            2305843008139952128, 2658455991569831744654692615953842176,
            191561942608236107294793378084303638130997321548169216]

def perfect(n: int) -> list:
    """Returns the list of perfect numbers less than given number 'n'.
        In number theory, a perfect number is a positive integer that is equal to
        the sum of its positive divisors, excluding the number itself.
        For instance, 6 has divisors 1, 2 and 3 (excluding itself), and 1 + 2 + 3 = 6,
        so 6 is a perfect number.

    Args:
        n (int): A number to find perfect numbers less than it.

    Returns:
        list: The list of perfect numbers less than given number 'n'.
    """
    pfs = perfects
    if n > pfs[-1]:
        for i in range(pfs[-1] , int(float(n)+1)):
            if sum(factors(i)[:-1]) == i: # type:ignore
                pfs.append(i)
        return pfs
    else:
        pfs = list(filter(lambda x: x <= n, perfects))
        return pfs
            
def isperfect(n: int) -> bool:
    """Checks whether a number is Perfect or not and returns a boolean.

    Args:
        n (int): Number to be checked.

    Returns:
        bool: True if Perfect, False if not.
    """
    if n > perfects[-1]:
        if sum(factors(n)[:-1]) == n: # type:ignore
            return True
    return True if n in perfects else False

def isprime(n : int) -> bool:
    """Checks whether a number is Prime or Not. In case of Prime number returns 'True', otherwise returns 'False'.

    Args:
        n (int): The number to be checked for being prime or not.

    Returns:
        bool: True for Prime numbers, False for non-Prime numbers.
    """
    try:
        assert isinstance(n,int) and n > 1
        for i in range(2, int((n**(0.5))+1)):
            if n % i == 0:
                return False
        return True
    except AssertionError:
        return False
    
def isemirp(n : int) -> bool:
    """Checks whether a number and its reversed form is Prime or Not.
    In case both of them are Prime numbers returns 'True', otherwise returns 'False'.

    Args:
        n (int): The number to be checked itself and its reversed form for being prime or not.

    Returns:
        bool: True if both of them are Prime numbers, False otherwise.
    """
    try:
        assert isinstance(n,int) and n > 1
        return isprime(int(str(n)[::-1])) and isprime(n)
    except AssertionError:
        return False

def get_primes(l : int | float, u : int | float = 0) -> list[int]:
    """Returns the list of Prime numbers in interval [l,u].
        If you enter only one argument (l) then the interval will be: [0,l].

    Args:
        l (int | float): Lower bound of interval. (or if one argument passed, it will be the Upper bound)
        u (int | float, optional): Upper bound of interval. (or if one argument passed, it will be set Default to 0 and becomes the Lower bound)

    Returns:
        list[int]: List of Prime numbers in the given interval.
    """
    if u == 0:
        l,u = 0,l
    if u < l:
        u,l = l,u
    primes = []
    for num in range(ceil(float(l)), int(float(u))+1):
        if isprime(num):
            primes.append(num)
    return primes

def gcd(*integers : int) -> int | str:
    """Greatest Common Divisor (Factor)

    Returns:
        int: Greatest Common Divisor of integers
    """
    try:
        for value in integers:
            assert isinstance(value,int)
        GCDset = set(factors(integers[0]))
        for integer in integers[1:]:
            GCDset = GCDset & set(factors(integer))
        return max(GCDset)
    except AssertionError:
        return "AssertionError: Given numbers must be Integers"
hcf = gcf = gcd
""" Greatest (Highest) Common Factor (Divisor)
    """
    # Euclid's algorithm:
    # a, b = max(a,b), min(a,b)
    # while a != b:
    #     return gcd(a-b, b)
    # return a
    
    # Euclidean algorithm:
    # while b != 0:
    #     a, b = b, a%b
    # return a
    
    # Euclidean algorithm using recurion:
    # if b == 0:
    #     return a
    # else:
    #     r = a % b
    # return gcd(b,r)
    

def iscoprime(*integers : int) -> bool | str:
    """Two numbers are called relatively prime, or coprime, if their greatest common divisor equals 1.

    Returns:
        bool: True if the numbers are Co-Prime, False if they're not.
    """
    try:
        for value in integers:
            assert isinstance(value,int)
        return True if gcd(*integers) == 1 else False
    except AssertionError:
        return "AssertionError: Given numbers must be Integers"
        

def lcm(*integers : int | float) -> int | float:
    """Least Common Multiple

    Returns:
        int: Least Common Multiple of given numbers
    """
    greater = max(integers)
    i = 0
    while True:
        i += 1
        if all(greater*i % integer == 0 for integer in integers):
            lcm = greater*i
            break
    return lcm

    # greater = max(a,b)
    # i = 0
    # while True:
    #     i += 1
    #     if ((greater*i % a == 0) and (greater*i % b == 0)):
    #         lcm = greater*i
    #         break
    # return lcm
    

def bezout(a : int, b : int):
    """Bézout's identity — Let a and b be integers with greatest common divisor d.
    Then there exist integers x and y such that ax + by = d.
    Moreover, the integers of the form az + bt are exactly the multiples of d.

    Args:
        a (int): First Integer Number
        b (int): Second Integer Number

    Returns:
        dict: A dictionary containing the: ("GCF" , "Bezout_identity", a , b) as 'KEYS' and their corresponding values.
        The values for 'a' and 'b' are their corresponding coefficients of Bézout's identity.
    """
    try:
        assert (isinstance(a,int) and isinstance(b,int))
        a , b = max(a,b) , min(a,b)
        first , second = a , b
        s = [1,0]
        t = [0,1]
        while b != 0:
            q = a // b
            s.append(s[-2] - q * s[-1])
            t.append(t[-2] - q * t[-1])
            a , b = b , a%b
        # Testing the Bézout's identity (Not necessary because the algorithm works fault-free and never Fails)
        if s[-2] * first + t[-2] * second != a: 
            return "Failed"
        else:
            return {"GCF": a, 
                    "Bezout_identity": f"({s[-2]})x({first}) + ({t[-2]})x({second}) = {a}", 
                    second: t[-2],
                    first: s[-2]}
    except AssertionError:
        return "AssertionError: Given numbers must be Integers"
    
def goldbach(n : int) -> list | str:
    """Every even integer greater than 2 can be written as the sum of two prime numbers.

    Args:
        n (int): An even integer greater than 2.

    Returns:
        list | str: A list of tuples containing prime numbers whose sum is equal to 'n'
    """
    try:
        assert (isinstance(n,int) and not n % 2 and n != 2)
        gbach = []
        for i in range(2, n//2 + 1):
            if isprime(i) and isprime(n-i):
                gbach.append((i,n-i))
        if not gbach:
            return f"Goldbach's Conjecture has been rejected by number {n}"
        return gbach
    except AssertionError:
        return "Please enter an even integer greater than 2 to check for Goldbach's Conjecture."
    
def legendre(n : int):
    try:
        assert (isinstance(n, int) and n > 0)
        plist = []
        for i in range(n**2, (n+1)**2):
            if isprime(i):
                plist.append(i)
        return plist if plist else f"The Legendre's Conjecture has been disproved by number {n}"
    except AssertionError:
        return "Please input a positive integer to check for Legendre's Conjecture"
    
def firstprimes(n : int, start : int | float = 2):
    try:
        assert (isinstance(n, int) and n > 0)
        if start <= 2:
            fplist = [2]
            start = 3
            n -= 1
        else:
            if not isinstance(start, int):
                start = ceil(start)
            if not start % 2:
                start += 1
            fplist = []
        while n:
            if isprime(start):
                fplist.append(start)
                n -= 1
            start += 2
        return fplist
    except AssertionError:
        return "Please input a positive integer for 'n'"

def distance(x1 : int | float, y1 : int | float , x2 : int | float = 0, y2: int | float = 0) -> float:
    """Returns the Pythagorean distance between two points (x1,y1) and (x2,y2)

    Args:
        x1 (int | float): x-coordinate of first point
        y1 (int | float): y-coordinate of first point
        x2 (int | float, optional): x-coordinate of second point. Defaults to 0.
        y2 (int | float, optional): y-coordinate of second point. Defaults to 0.

    Returns:
        float: The Pythagorean distance
    """
    return ((x2-x1)**2 + (y2-y1)**2)**(1/2)

def ispythagorean(s1 : int | float, s2 : int | float, s3 : int | float) -> bool:
    """Checks whether the inputted triplet forms a right triangle or not.

    Args:
        s1, s2, s3: 3 sides of a triangle

    Returns:
        True if yes, False otherwise.
    """
    return True if ((s1**2 + s2**2 == s3**2) or (s1**2 + s3**2 == s2**2) or (s2**2 + s3**2 == s1**2)) else False

def digitsum(n : int | float | str) -> int:
    """Adds the digits of the given number

    Args:
        n (int | float | str): A number either in int, float or inside a string.

    Returns:
        int: Sum of the digits.
    """
    if isinstance(n,float):
        n = "".join(str(n).split("."))
    elif isinstance(n,str):
        n = float(n)
        n = "".join(str(n).split("."))
    sum = 0
    n = abs(int(n))
    while n > 0:
        sum += n % 10
        n = n//10
    return sum

def yielddigitsum(n : int | float | str):
    """Continuously adds the digits of the given number and repeats the process till one digit number

    Args:
        n (int | float | str): A number either in int, float or inside a string.

    Returns:
        int: Sum of the sums of digits to one digit number.
    """
    while len(str(n)) != 1:
        yield n
        n = digitsum(n)
    yield n

def isfactorialable(n : int) -> int | bool:
    """Checks whether the input number is the factorial of any Natural number or not. if it is, will return the number, else returns False.

    Args:
        n (int): The number to be checked if it is the result of factorial of any natural number.

    Returns:
        int : the number that it's factorial is equal to 'n'.
        bool: False if 'n' isn't a result of factorial of any natural number.
    """
    i = 2
    while n % i == 0:
        n = n // i
        i += 1
        if n == i:
            return i
    return False

digs = string.digits + string.ascii_letters

def int2base(x : int, base : int = 2) -> str:
    """Converts the base 10 integer to desired base.
    desired base default is 2.
    (inverse of builtin "int" function) : int(km.int2base(x, base), base) = x
    
    Args:
        x (int): Base 10 integer input.
        base (int): Desired base (Default is 2)

    Returns:
        str: Equivalent number in desired base
    """
    if x < 0:
        sign = -1
    elif x == 0:
        return digs[0]
    else:
        sign = 1

    x *= sign
    digits = []

    while x:
        digits.append(digs[x % base])
        x //= base

    if sign < 0:
        digits.append('-')

    digits.reverse()

    return ''.join(digits)

def base2base(x : int | str , x_base : int , to_base : int) -> str:
    """Will be removed in future updates. Replaced with new function convert() which handles the floats as well.
        Converts the input number "x" from the given base "x_base" to desired base "to_base".

    Args:
        x (int | str): The input number in any base
        x_base (int): The base of x
        to_base (int): Desired base

    Returns:
        str: Equivalent number in desired base
    """
    if isinstance(x,int):
        x = str(x)
    return int2base(int(x, x_base),to_base)

# "base2base" only converts the whole numbers, while "convert" is the general form which converts any floating point number in any base to
# desired base.

def convert(x : int | float | str, x_base : int, to_base : int, precision : int = 0) -> str:
    """Converts the input number "x" from the given base "x_base" to desired base "to_base".

    Args:
        x (int | float | str): The input number in any base
        x_base (int): The base of x
        to_base (int): Desired base
        precision (int, optional): Number of digits after floating point. Defaults to number of the input number "x" digits after floating point.

    Returns:
        str: Equivalent number string in desired base
    """
    # Converting the number into String:
    if not isinstance(x, str):
        x = str(x)
        
    # Joining the Integer part and Fractional part together and then Converting From x_base to base 10:
    integral, point, fractional = x.strip().partition('.')
    num = int(integral + fractional, x_base) * x_base ** -len(fractional)

    # Converting From base 10 to new base:
    # # Number of digits after floating point:
    precision = len(fractional) if not precision else precision
    # # Using int2base() Function to convert from base 10 to new base:
    s = int2base(int(round(num / to_base ** -precision)), to_base)
    if precision:
        return s[:-precision] + '.' + s[-precision:]
    else:
        return s

def deconvert(x : int | float | str, y : int | float | str) -> list:
    """Brute Forces to find any occurences of number x in any base from 2 to 16 being equal to number y in any base from 2 to 16.

    Args:
        x (int | float | str): First input number in any base from 2 to 16
        y (int | float | str): Second input number in any base from 2 to 16

    Returns:
        list: A list of literals showing 'x_(from_base) = y_(to_base)'
    """
    if not isinstance(x, str):
        x = str(x)
    if not isinstance(y, str):
        y = str(y)
    dc = []
    for i in range(2,17):
        for j in range(2, 17):
            try:
                if convert(x,i,j) == y:
                    dc.append(f"({x})_{i} = ({y})_{j}")
            except ValueError:
                pass
    return dc
    
def collatz(num : int | float) -> list:
    """Collatz Conjecture Extended to all Real Numbers.
    
    The Collatz Conjecture: (Algorithm)  
        ● if the number is even, divide it by two.  
        ● if the number is odd, triple it and add one.  
        
    The Collatz Conjecture: (Description)  
        ● The Conjecture is that if the input number is a POSITIVE INTEGER: these sequences always reach 1,
          no matter which positive integer is chosen to start the sequence.
          
    The Collatz Conjecture Extention to larger domain (All REAL NUMBERS):  
        ● There will be five different cycles:  
            (1) 0, 0, ...: Zero Cycle  
            (2) 1, 4, 2, 1, ...: Positive Real Numbers Cycle  
            (3) -1, -2, -1, ...  
            (4) -5, -14, -7, -20, -10, -5, ...  
            (5) -17, -50, -25, -74, -37, -110, -55, -164, -82, -41, -122, -61, -182, -91, -272, -136, -68, -34, -17, ...  

    Args:
        num (int | float): The input real number.

    Returns:
        list: The Collatz Cycle of the input number.
    """
    res = [num]
    while True:
        if num == 1:
            break
        elif num % 2 == 0:
            num /= 2
        else:
            num = 3 * num + 1
        
        if num in res:
            break 
        else:
            res.append(num)
    return res

def continued_fraction(numerator: int | float, denominator : int | float = 1, precision : int = 19) -> list[int]:
    """A continued fraction is an expression obtained through an iterative process of representing a number
        as the sum of its integer part and the reciprocal of another number, then writing this other number
        as the sum of its integer part and another reciprocal, and so on.
        The process ends when either the last number becomes an integer without any fractional part or the precision reachs out.
        
    Example: Consider the rational number (415/93) which is around 4.4624.
        Now we write this number in an iterative process of representing it as:
        [Integer Part] + [Reciprocal of the rest]
        415/93 = 4 + 43/93
               = 4 + 1 / (93/43) 
               = 4 + 1 / (2 + 7/43)
               = 4 + 1 / (2 + (1 / (43/7))) 
               = 4 + 1 / (2 + 1 / (6 + 1/7))
               # Process ends here as the last denominator 7 is an integer and doesn't have fractional part.
        This can be represented by the abbreviated notation:
        415/93 = [4, 2, 6, 7]
        

    Args:
        numerator (int | float): Numerator of the fraction.
        denominator (int | float, optional): Denominator of the fraction. Defaults to 1.
        precision (int, optional): Max number of iterations. Defaults to 19.

    Returns:
        list[int]: List of integer parts in the continued fraction in abbreviated notation.
    """
    frac = numerator / denominator
    fraclist = [int(frac)]
    while precision and (round(frac,9) != fraclist[-1]):
        temp = frac - fraclist[-1]
        frac = 1/(frac - fraclist[-1])
        # mitigating the floating point arithmetic errors for when frac should be an integer
        fraclist.append(int(frac) if frac-int(frac)<0.999999999 else round(frac))
        # mitigating the floating point arithmetic errors for sqrt(2), phi, ...
        if round(temp,10) == round(frac - fraclist[-1],10): 
            frac = temp + fraclist[-1]
        precision -= 1
    return fraclist

def contfrac_to_frac(seq):
    """ Converts the simple continued fraction in "seq" which is:
        (List of integer parts in the continued fraction in abbreviated notation)
        into a fraction, num / den
    """
    n, d, num, den = 0, 1, 1, 0
    for u in seq:
        n, d, num, den = num, den, num*u + n, den*u + d
    return num, den

def egyptian(num: int, den: int = 1) -> list[Fraction]:
    """Egyptian Fractions implemented using Greedy Algorithm.
    
    * Note: 'fractions' module has been used to implement the Fraction instances.
    
    * Definitions:
        - Unit Fraction: Any Fraction with 1 as its numerator and a positive integer for the denominator.
        - Egyptian Fraction: An Egyptian Fraction is the sum of finitely distinct Unit Fractions.
        - Greedy Algorithm: A Greedy Algorithm is any algorithm that follows the problem-solving heuristic
                            of making the locally optimal choice at each stage.
    * Description:
    For any rational number (m/n), where m,n are integers, there exist Unit Fractions (with numerator 1, 
    and integer denominators) such that the sum of these Unit Fractions is equal to (m/n).
    
    * Examples:
    Fraction(2,3) = Fraction(1,2) + Fraction(1,6)
    Fraction(7,12) = Fraction(1,2) + Fraction(1,12)
    Fraction(16,25) = Fraction(1,2) + Fraction(1,8) + Fraction(1,67) + Fraction(1,13400)

    Args:
        num (int): The numerator of Fraction. can use 'str' or 'float' if only one argument is passed. (den=1 by default)
        den (int, optional): The denominator of Fraction. Defaults to 1.

    Returns:
        list[Fraction]: a list of Unit Fractions which their sum is equal to the given Fraction.
        * sum(egyptian(m,n)) = m/n
    """
    neg = 0
    if den == 1:
        curr = Fraction(num)
    else:
        curr = Fraction(num,den)
    if curr < 0:
        neg = 1
        curr *= -1
    if curr.numerator == 1:
        return [Fraction(1,den+1), Fraction(1,den*(den+1))]
    nxt = Fraction(1,ceil(1/curr))
    egpt = [curr, nxt,]
    while curr.numerator != 1:
        curr -= nxt
        nxt = Fraction(1,ceil(1/curr))
        egpt.append(nxt)
    if neg:
        egpt = [-1*i for i in egpt]
    return egpt[1:]

def cnrt(z: complex, n: int) -> dict:
    """Complex nth root theorem:
        Any nonzero complex number has exactly n ∈ N distinct nth roots.
         The roots lie on a circle of radius |z| centered at the origin and spaced out evenly by angles of 2π/n.
         Concretely, if z = complex(x,y) = x + iy = r cis(θ) = r e^{iθ}, then solotions to c^n = z are given by:
         c = z**(1/n) = r**(1/n) cis ((2kπ+θ)/n) = r**(1/n) e^(i((2kπ+θ)/n)) for k ∈ {0,1,2,...,n-1}

    Args:
        z (complex): Any Complex Number
        n (int): The order of number (The degree of roots)

    Returns:
        dict: A dictionary consisting of all of the nth roots for k ∈ {0,1,2,...,n-1}
    """
    real = z.real
    imag = z.imag
    r = (real**2 + imag**2)**(1/2)
    # Finding Angle theta (We need theta such that: "0 = theta <= 2*pi" but using atan we have: "-pi/2 < atan < pi/2")
    # Origin = (0,0):
    if real == 0 and imag == 0:
        theta  = 0
    # Positive X-axis and Quadrant I:
    elif real > 0 and imag >= 0:
        theta = atan(imag / real)
    # Positive Y-axis:
    elif real == 0 and imag > 0:
        theta = pi/2
    # Quadrant II and III:
    elif real < 0 and imag != 0:
        theta = atan(imag / real) + pi
    # Negative X-axis:
    elif real < 0 and imag == 0:
        theta = pi
    # Negative Y-axis:
    elif real == 0 and imag < 0:
        theta = 3*pi/2
    # Quadrant IV:
    else:
        theta = atan(imag / real) + 2*pi 
    nroot = []
    for k in range(n):
        # nroot.append((r**(1/n))*((cos(((2*k*pi) + theta)/n)) + ((-1)**(1/2))*sin(((2*k*pi) + theta)/n)))
        nroot.append((r**(1/n))*(e**(((-1)**(1/2))*(((2*k*pi) + theta)/n))))
    return {i : nroot[i] for i in range(len(nroot))}

if __name__ == "__main__":
    print("Module is being run standalone by the user")
    # print(sieve(1000,7))
    # print(sieve.__annotations__)
    # help(sieve)
    # print(truncate(m.pi,3))
    
else:
    print("Module has been imported externally")
    
    
# if round(frac,10) == int(frac): # mitigating the floating point arithmetic errors
#     break