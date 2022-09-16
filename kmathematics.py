import string
from fractions import Fraction
from math import ceil

def trunc(number : float, digits : int = 0) -> float:
    """truncate will cut a float number at a given decimal without rounding. (round toward zero)

    Args:
        number (float): the number to be cut
        digits (int): number of decimals, Defaults to 0 decimals

    Returns:
        float: truncated version of number with n decimals
    """ 
    try:
        assert isinstance(digits, int)
        pow10 = 10 ** digits
        return number * pow10 // 1 / pow10
    except AssertionError:
        return "Second argument 'n' needs to be an Integer number"


def sieve(n : int | float, p : int = None) -> list[int]:
    """Sieve of Eratosthenes:
        It is an ancient algorithm for finding all prime numbers up to any given limit. 
        It does so by iteratively marking as composite (i.e., not prime) the multiples of each prime,
        starting with the first prime number, 2.

    Args:
        n (int | float): Any number greater than or equal to 2.
        Note: input only 1 argument (n) to output the list of primes up to 'n'.
        
        p = None : (Optional) should be a prime number.
        Note: if 2 arguments are input, it will return the list of composite numbers eliminated in removing the multiples of 'p'.

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


def prime_factorise(n: int) -> list[int]:
    """Prime Factorise (Trial Division): Gives a list of all Prime Factors of 'n'.

    Args:
        n (int): Number to Prime Factorise.

    Returns:
        list[int]: List of all Prime Factors of the given number.
    """
    try:
        assert isinstance(n,int) and n != 0
        prime_factors = []
        if n < 0:
            prime_factors.append(-1)
            n = abs(n)
        while n % 2 == 0: # factor out all occurenses of 2
            prime_factors.append(2)
            n //= 2
        f = 3
        while f * f <= n: # factoring out other odd primes
            if n % f == 0:
                prime_factors.append(f)
                n //= f
            else:
                f += 2
        if n != 1: prime_factors.append(n)
        return prime_factors
    except AssertionError:
        return "AssertionError: trial_division function needs a non-zero integer number as input"
        

def factors(n: int) -> list[int]:
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
        return "AssertionError: factors function needs a non-zero integer number as input"


def perfect(n):
    for i in range(1 , n+1):
        if sum(factors(i)[:-1]) == i:
            yield i
            
def isperfect(n):
    if sum(factors(n)[:-1]) == n:
        return True
    return False

def isprime(n : int) -> bool:
    """Checks whether a number is Prime or Not. in case of Prime number returns 'True', otherwise returns 'False'

    Args:
        n (int): The number to be checked for being prime or not.

    Returns:
        bool: True for Prime numbers, False for non-Prime numbers.
    """
    try:
        assert isinstance(n,int) and n >= 1
        if n == 1:
            return False
        for i in range(2,int((n**(0.5))+1)):
            if n % i == 0:
                return False
        return True
    except AssertionError:
        return "AssertionError: isprime function needs a Natural number: {1,2,3,4,...} as input"


def get_primes(n : int | float):
    """Generates the Prime numbers less than or equal to the given number

    Args:
        n (int | float): Number which will be the upper bound of prime numbers

    Returns:
        _type_: Generator object

    Yields:
        _type_: Prime numbers
    """
    if n == 0 or n == 1:
        return None
    else:
        for num in range(2,int(n)+1):
            if isprime(num):
                yield num


def gcd(*integers : int) -> int:
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
        print("AssertionError: Given numbers must be Integers")
gcf = gcd
""" Greatest Common Factor (Divisor)
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
    

def iscoprime(*integers : any) -> bool:
    """Two numbers are called relatively prime, or coprime, if their greatest common divisor equals 1.

    Returns:
        bool: True if the numbers are Co-Prime, False if they're not.
    """
    try:
        for value in integers:
            assert isinstance(value,int)
        return True if gcd(*integers) == 1 else False
    except AssertionError:
        print("AssertionError: Given numbers must be Integers")
        

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
    

def bezout(a : int, b : int) -> dict[str | int , any]:
    """Bézout's identity — Let a and b be integers or polynomials with greatest common divisor d.
    Then there exist integers or polynomials x and y such that ax + by = d.
    Moreover, the integers or polynomials of the form az + bt are exactly the multiples of d.

    Args:
        a (int): First Integer Number
        b (int): Second Integer Number

    Returns:
        dict[str | int , Any]: A dictionary containing the: ("GCF" , "Bezout_identity", a , b) as 'KEYS' and their corresponding values.
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
            print("Failed")
        else:
            return {"GCF": a, 
                    "Bezout_identity": f"({s[-2]})x({first}) + ({t[-2]})x({second}) = {a}", 
                    first: s[-2], 
                    second: t[-2]}
    except AssertionError:
        print("AssertionError: Given numbers must be Integers")


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
    """Converts the input number "x" from the given base "x_base" to desired base "to_base".

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

def convert(x : int | float | str, x_base : int, to_base : int, precision : int = None) -> str:
    """Converts the input number "x" from the given base "x_base" to desired base "to_base".

    Args:
        x (int | float | str): The input number in any base
        x_base (int): The base of x
        to_base (int): Desired base
        precision (int, optional): Number of digits after floating point. Defaults to number of the input number "x" digits after floating point.

    Returns:
        str: Equivalent number string in desired base
    """
    #from original_base
    if not isinstance(x, str):
        x = str(x)
    integral, point, fractional = x.strip().partition('.')
    num = int(integral + fractional, x_base) * x_base ** -len(fractional)

    #to new_base
    precision = len(fractional) if precision is None else precision
    s = int2base(int(round(num / to_base ** -precision)), to_base)
    if precision:
        return s[:-precision] + '.' + s[-precision:]
    else:
        return s

def deconvert(x, y):
    if not isinstance(x, str):
        x = str(x)
    if not isinstance(y, str):
        y = str(y)
    for i in range(2,17):
        for j in range(2, 17):
            try:
                if convert(x,i,j) == y:
                    print(f"({x})_{i} = ({y})_{j}")
            except ValueError:
                pass
    
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
    ''' Convert the simple continued fraction in `seq`
        into a fraction, num / den
    '''
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