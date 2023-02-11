from dataclasses import dataclass,field
from typing import Tuple,Union
from array import array
from abc import ABC,ABCMeta,abstractmethod
from abc import abstractclassmethod
from math import prod,gcd,lcm
from random import randint
from functools import partial,update_wrapper
from copy import copy,deepcopy
from mmath.exceptions import *

class __saver(ABC,metaclass=ABCMeta):
	"""
		Id saver for classes. You Cannot make Instance from this class.
	"""
	@abstractmethod
	def __init__(self):pass
	
	@abstractclassmethod
	def __init_subclass__(cls):raise NotImplementedError("Cannot subclass __saver!")
	identifiers=array("h",[])
	clone_maker=lambda a:deepcopy(a)
	clone_make_flow=partial(clone_maker)
	update_wrapper(clone_make_flow,clone_maker)
	class Base(object):pass
class Math(__saver.Base):
	"""
		Class for mathematical operations.
	
	"""
	def multiply(self,a,b):
		return a*b
		
	def divide(self,a,b,round=False):
		return a//b if round else a/b
		
	def add(self,a,b):
		return a+b
		
	def subtract(self,a,b,absolute=False):
		return abs(a-b) if absolute else a-b
		
	def power(self,a,b,c=1):
		return (a**b)%c
		
	def root(self,a,b=2):
		return a**(1/b)
		
	def triangle(self,a):
		return sum(range(a+1))
		
	def pythagorean(self,a,b):
		return (a**2+b**2)**(1/2)
		
	def factorial(self,a):
		return math.prod(range(1,a))
	
	def log(self,a,base):
		return [i for i in range(a) if self.power(base,i)==a][0]
	
@dataclass(order=True)
class Number(__saver.Base):
    """
    	Class for numbers like builtin <int> class.
    
    """
    m=0
    __value:int

    def __init__(self):
        Number.m+=1
    def __add__(self,other):
        return self.__value+other.__value
    def __sub__(self,other):
        return self.__value-other.__value
    def __mul__(self,other):
        return self.__value*other.__value
    def __truediv__(self,other):
        return self.__value/other.__value
    def __floordiv__(self,other):
        return self.__value//other.__value
    def __pow__(self,other):
        return pow(self.__value,other.__value)
    def __ne__(self,other):
        return self.__value != other.__value or self.__class__ == other.__class__     
    def __iadd__(self,other):
        return Number(self.__value+other.__value)
    def __isub__(self,other):
        return Number(self.value-other.__value)
    def __imul__(self,other):
        return Number(self.__value*other.__value)
    def __itruediv__(self,other):
        return Number(self.__value/other.__value)
    def __ifloordiv__(self,other):
        return Number(self.__value//other.__value)
    def __ipow__(self,other):
        return Number(self.__value**other.__value)
    def __bool__(self):
        return self.__value!=0
    def as_integer(self):
        return int(self.__value)
        
class Fraction(__saver.Base):
  	def __init__(self,num,den,*,simplify=True):
  		self.numerator=num
  		self.denominator=den
  		if simplify:
  			self.simplify()
  			
  	def __repr__(self):
  		return f"<Fraction ({self.numerator},{self.denominator})>"
  		
  	def __add__(self,other):
  		den=lcm(self.denominator,other.denominator)
  		m1,m2=den//self.denominator,den//other.denominator
  		num=self.numerator*m1+other.numerator*m2
  		return self.__class__(num,den)
  		
  	def __sub__(self,other):
  		den=lcm(self.denominator,other.denominator)
  		m1,m2=den//self.denominator,den//other.denominator
  		num=self.numerator*m1-other.numerator*m2
  		return self.__class__(num,den)  
  		
  	def __mul__(self,other):
  		num=self.numerator*other.numerator
  		den=self.denominator*other.denominator
  		return self.__class__(num,den).simplify()
  				
  	def __truediv__(self,other):
  		t=other.reverse()
  		return self*t
  		
  	def simplify(self):
  		simp=gcf(self.numerator,self.denominator)
  		self.numerator//=simp
  		self.denominator//=simp
  		return self
  		
  	def reverse(self):
  		return self.__class__(self.denominator,self.numerator)
 
class MultiplicationTable(__saver.Base):
 	def __init__(self,point):
 		self.point=point
 		self.quadline=range(1,self.point+1)
 		self.tripline=range(1,self.point+1)		
 		
 	def __repr__(self):
 		text=""
 		for x in self.quadline:
 			for y in self.tripline:
 				text+=str(x*y)+"\t"
 			text+="\n" 			
 		return text
 	
 	def __getitem__(self,p:Tuple):
 		x,y=p
 		if x in self.quadline and y in self.tripline:
 			return x*y
 		raise TableError("Unexpected dimensions.",cr=2)
 		
 	def __setitem__(self,tpl,val):
 		if prod(tpl)==val and tpl[0] in self.quadline and tpl[-1] in self.tripline:
 			return super().__init__()
 		raise TableError("Unequal ndim separator!",cr=2)
 		 			
@dataclass(order=True)
class Point(__saver.Base):
	"""
		Class For points in Cartesian coordinate system diagram.
	
	"""
	___id: int =field(init=False,compare=False,repr=False)
	
	x: int = 0
	y: int = 0
	def __post_init__(self):
		id=self.__set_id()
		while not id:
			id=self.__set_id()
		self.t_pointer : Tuple[int,int]=(self.x,self.y)
	def __set_id(self):
		id=randint(-32768,32768)
		ids=__saver.identifiers
		if id not in ids:
			ids.append(id)
			self.___id=id
			return id
		return False
	def __sub__(self,other):
		x_dis=abs(self.x-other.x)
		y_dis=abs(self.y-other.y)
		dist=Math.pythagorean(self,x_dis,y_dis)
		return Number(dist)

class Vector(__saver.Base):
		"""
			Class for Vectors in Cartesian Coordinate System.
		"""
		def __init__(self,head:Union[Point,Number],tail:Union[Point,Number],*,ht=False):
			if ht:
				self.head=head
				self.tail=tail
				x_line:Number=Number(abs(self.head.x-self.tail.x))
				y_line:Number=Number(abs(self.tail.y-self.head.y))
			else:
				self.x_line=head
				self.y_line=tail
			
		def __repr__(self):
			return f"{self.__class__.__name__}(x={self.x_line},y={self.y_line})"
		def __eq__(self,other):
			return self.y_line==other.y_line and self.x_line == other.x_line and type(self) == type(other)
		def __reversed__(self):
			return self.__class__(x_line=self.x_line*-1,y_line=self.y_line*-1)
			
		def __add__(self,other):
			new_v=Vector(x_line=self.x_line+other.x_line,y_line=self.y_line+other.y_line)
			return new_v
		
		def __sub__(self,other):
			return self+reversed(other)
			
		def __mul__(self,other):
			return Number((self.x_line*other.x_line)+(self.y_line*other.y_line))
