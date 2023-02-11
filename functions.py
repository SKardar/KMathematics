from mmath.exceptions import  MathFunctionError as MTError
from mmath.classes import Math
from mmath.constants import pi,e
from itertools import count
from string import digits,ascii_letters
from math import ceil,sin,cos,atan
from fractions import Fraction
from mmath.exceptions import ArgumentError as AError
from copy import *

def __num(num):
	try:
		return int(num)
	except Exception as e:
		raise AError("Please send num to this function.",cr=2) from e
		
def sum_squares(num):
	__num(num)
	f=0
	for d in str(num):
		f+=int(d)**2
	return f

def is_happy(num):
	__num(num)
	m=str(num)
	items=[]
	while m!=1 and sum_squares(m) not in items:
		items.append(sum_squares(int(m)))
		m=items[-1]
	return items[-1]==1

def factors(num):
	__num(num)
	factorlist=[]
	for i in range(1,num+1):
		if not num%i:
			factorlist.append(i)
	return factorlist

def is_perfect(num):
	__num(num)
	return sum(factors(num))==num
	
def is_prime(num):
	__num(num)
	return len(factors(num))==2
	
def is_emirp(num):
	__num(num)
	return is_prime(int(str(num)[::-1])) and is_prime(num)
	
def is_pythagorean(l1,l2,hyp):
	__num(l1)
	__num(l2)
	__num(hyp)
	return Math.pythagorean(l1,l2)==hyp
	
def is_palindrome(num):
	__num(num)
	pal=str(num)[::-1]==str(num)
	return pal and not len(str(num))==1
	

def collatz_loop(num):
	__num(num)
	loop=[num]
	while 1:
		if num == 1:
			break
		else:
			num=num//2 if not num%2 else num*3+1
		if num in loop:
			break
		res.append(loop)
	return loop
	
def goldbach(num):
	__num(num)
	gbach=[]
	if num%2 or num==2:
		raise AError("num in goldbach conjecture must be even and greater than 2.",cr=1)
	for i in range(num>>1):
		if is_prime(i) and is_prime(num-i):
			gbach.append((i,num-i))
	if not gbach:
		raise MTError("You rejected goldbach conjecture!!!!",cr=3)
	return gbach
	
def legender(n1,n2):
	__num(n1);__num(n2)
	
	nsq1=n1**2
	nsq2=n2**2
	plist=[]
	for i in range(nsq1,nsq2+1):
		if isprime(i):
			plist.append(i)
	if plist:
		return plist
	raise MTError("You rejected legender conjecture!!!!",cr=3)
	
def calc(exp):
	return f"{exp}={eval(exp)}"
	
def infact(num):
	__num(num)
	return sum(i for i in range(num+1) if num%i)
			

def first_primes(num):
	__num(num)
	c=0
	fc=1
	while c < num:
		while True:
			if is_prime(fc):
				yield fc
				c+=1
				fc+=1
				break
			fc+=1
			
def first_primes(num):
	__num(num)
	c=0
	fc=1
	while c < num:
		while True:
			if is_prime(fc):
				yield fc
				c+=1
				fc+=1
				break
			fc+=1			
			
def first_happies(num):
	__num(num)
	c=0
	fc=1
	while c < num:
		while True:
			if is_happy(fc):
				yield fc
				c+=1
				fc+=1
				break
			fc+=1

def first_perfects(num):
	__num(num)
	c=0
	fc=1
	while c < num:
		while True:
			if is_perfect(fc):
				yield fc
				c+=1
				fc+=1
				break
			fc+=1
	
def first_emirps(num):
	__num(num)
	c=0
	fc=1
	while c < num:
		while True:
			if is_emirp(fc):
				yield fc
				c+=1
				fc+=1
				break
			fc+=1


def primes_before(num):
	__num(num)
	for i in range(1,num+1):
		if is_prime(num):
			yield num
			
def happies_before(num):
	__num(num)
	for i in range(1,num+1):
		if is_happy(num):
			yield num

def palindromes_before(num):
	__num(num)
	for i in range(1,num+1):
		if is_palindrome(num):
			yield num

def perfects_before(num):
	__num(num)
	for i in range(1,num+1):
		if is_perfect(num):
			yield num
			
def emirps_before(num):
	__num(num)
	for i in range(1,num+1):
		if is_emirp(num):
			yield num
	
def pythagoreans_before(num):
	__num(num)
	for r in range(1,num+1):
		for g in range(1,num+1):
			for b in range(1,num+1):
				if is_pythagorean(r,g,b):
					yield (r,g,b)
					
def cont_frac(numbers):
	num,den=numbers
	frac=num/den
	fraclist=[int(frac)]
	while frac != fraclist[-1] and round(frac,9) != fraclist[-1]:
		frac2 = frac-fraclist[-1]
		frac=pow(frac2,-1)
		fraclist.append(int(frac) if frac-int(frac) < 0.9999999 else round(frac))
	return fraclist

def arithmetic(*args):
	n=len(args)
	m=0
	for item in args:
		m+=item
	m/=2
	return m
		
def geometric(*args):
	n=len(args)
	m=1
	for item in args:
		m*=item	
	m**=(1/n)
	return m
	
def harmonic(*args):
	n=len(args)
	m=0
	for item in args:
		m+=(1/item)	
	n/=m
	return n
	
def root_square(*args):
	n=len(args)
	m=0
	for item in args:
		m+=(item**2)
	m/=n
	m**=(1/2)
	return m	
	
def gcd(*args):
    cds=set(factors(args[0]))
    for i in args[1:]:
    	cds=cds.union(set(factors(i)))
    return max(cds)
    	
hcf = gcf = gcd


def iscoprime(*args) -> bool:
	return gcd(*args) == 1     

def lcm(*args):
	c=max(args)
	i=1
	while any(c*i%d for d in args):i+=1
	return c*i
	
def rad_to_deg(rads):
	return rads*180/pi
	
def deg_to_rad(degs):
	return degs*pi/180
#########################	
def contfrac_to_frac(seq):
    n, d, num, den = 0, 1, 1, 0
    for u in seq:
        n, d, num, den = num, den, num*u + n, den*u + d
    return num, den   

def digitsum(n):
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

def yielddigitsum(n):
    while len(str(n)) != 1:
        yield n
        n = digitsum(n)
    yield n 
    
def trunc(num, digits=0):
    __num(digits);__num(num)
    pow10 = 10 ** digits
    return number * pow10 // 1 / pow10
    
def prime_factorise(n):
    prime_factors = []
    if n < 0:
        prime_factors.append(-1)
        n = abs(n)
    while n % 2 == 0:
        prime_factors.append(2)
        n //= 2
    f = 3
    while f * f <= n:
         if n % f == 0:
            prime_factors.append(f)
            n //= f
         else:
            f += 2
    if n != 1: 
        prime_factors.append(n)
    return prime_factors
    
   
   
def sieve(n, p=None):
    if n >=2:
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
def bezout(a, b):
     if isinstance(a,int) and isinstance(b,int):
        a , b = max(a,b) , min(a,b)
        first , second = a , b
        s = [1,0]
        t = [0,1]
        while b != 0:
            q = a // b
            s.append(s[-2] - q * s[-1])
            t.append(t[-2] - q * t[-1])
            a , b = b , a%b
        if s[-2] * first + t[-2] * second != a: 
            print("Failed")
        else:
            return {"GCF": a, 
                    "Bezout_identity": f"({s[-2]})x({first}) + ({t[-2]})x({second}) = {a}", 
                    first: s[-2], 
                    second: t[-2]}

def is_factorial(n):
    i = 2
    while n % i == 0:
        n = n // i
        i += 1
        if n == i:
            return i
    return False

digs = digits + ascii_letters

def int2base(x, base):
    if x < 0:
        sign = -1
    elif x == 0:
        return 0
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

def base2base(x, x_base, to_base):
    return int2base(int(str(x), x_base),to_base)

# "base2base" only converts the whole numbers, while "convert" is the general form which converts any floating point number in any base to

def convert(x, x_base, to_base, precision=None):
    if not isinstance(x, str):
        x = str(x)
    integral, point, fractional = x.strip().partition('.')
    num = int(integral + fractional, x_base) * x_base ** -len(fractional)

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
def cnrt(z, n):
    
    real = z.real
    imag = z.imag
    r = (real**2 + imag**2)**(1/2)
    
    if real == 0 and imag == 0:
        theta  = 0
    
    elif real > 0 and imag >= 0:
        theta = atan(imag / real)
    
    elif real == 0 and imag > 0:
        theta = pi/2
    
    elif real < 0 and imag != 0:
        theta = atan(imag / real) + pi
    
    elif real < 0 and imag == 0:
        theta = pi
    
    elif real == 0 and imag < 0:
        theta = 3*pi/2
 
    else:
        theta = atan(imag / real) + 2*pi 
    nroot = []
    for k in range(n):
        nroot.append((r**(1/n))*(e**(((-1)**(1/2))*(((2*k*pi) + theta)/n))))
    return {i : nroot[i] for i in range(len(nroot))}
    
def egyptian(num: int, den: int = 1) -> list[Fraction]:
    neg = 0
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
    
    
