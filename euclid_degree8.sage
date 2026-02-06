#Code for computing the Euclidean division in some degree 8 number fields, for ISSAC submission #4295

F15.<z15>=CyclotomicField(15)
F20.<z20>=CyclotomicField(20)
F24.<z24>=CyclotomicField(24)
F8.<i,s3,s5>=NumberField([x^2+1,x^2-3,x^2-5])

E=[(vector([1,-1,-1,-1,-1,-1,-1,1])/2),
	vector([1,1,0,0,0,0,0,0]),
	vector([-1,1,0,0,0,0,0,0]),
	vector([0,-1,1,0,0,0,0,0]),
	vector([0,0,-1,1,0,0,0,0]),
	vector([0,0,0,-1,1,0,0,0]),
	vector([0,0,0,0,-1,1,0,0]),
	vector([0,0,0,0,0,-1,1,0])] #classical base of E8 root system
	
Einvmat=matrix([ #the inverse of the transpose of matrix(E)
	[   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  2],
	[ 1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,5/2],
	[-1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,1/2  ,7/2],
	[   0  ,  0  ,  1  ,  1  ,  1  ,  1  ,  1  ,  5],
	[   0  ,  0  ,  0  ,  1  ,  1  ,  1  ,  1  ,  4],
	[   0  ,  0  ,  0  ,  0  ,  1  ,  1  ,  1  ,  3],
	[   0  ,  0  ,  0  ,  0  ,  0  ,  1  ,  1  ,  2],
	[   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1  ,  1]
])	

	
def closest_vector_E8(x): #returns vector in E8 closest to x following Conway and Sloane (1982)
	def f(x): #rounds every coordinate of x, returning the floor in case of a half-integer 
		l=len(x)
		y=vector(ZZ,l)
		for i in range(l):
			u=round(x[i])
			if u-x[i]==0.5:
				y[i]=u-1
			else:
				y[i]=u
		return y
		
	def g(x): #rounds every coordinate of x, but rounds the worst one badly
		l=len(x)
		y=f(x)
		difflist=[abs(y[i]-x[i]) for i in range(l)]
		i=difflist.index(max(difflist))
		if y[i]>x[i]:
			y[i]-=1
		else:
			y[i]+=1
		return y
	l=len(x)
	u=f(x)
	v=g(x)
	if sum(u)%2==0:
		y=u
	else:
		y=v
	s=f([x[i]-1/2 for i in range(l)])
	t=g([x[i]-1/2 for i in range(l)])
	if sum(s)%2==0:
		z=vector([s[i]+1/2 for i in range(l)])
	else:
		z=vector([t[i]+1/2 for i in range(l)])
	if norm(y-x)>norm(z-x):
		return z
	else:
		return y

def basis_E8(F): #returns integral basis B of F for F=Q(zeta_n) for n=15,20,24 or F=Q(sqrt(-1),sqrt(3),sqrt(5)) such that the trace matrix of B is the E8 Cartan matrix, and base change matrix M from original basis
	if F==F15:
		zn=F.gens()[0]
		B=[zn^7, -zn^7 + zn^6 + zn^3 - zn^2, zn^5, -zn^5 - zn^3, zn^3 + zn, -zn^7 - zn^6 + zn^5 - zn^4 - 2*zn + 1, zn^7 + zn^6 + zn^4 + zn^2 + zn, zn^7 - zn^6 - zn^5 + zn^4 - zn^3 - 2]
		M=matrix([
			[2,0,-2,0,2,0,-1,1],
			[2,0,-2,0,3,0,-1,0],
			[3,1,-3,-1,4,1,-2,0],
			[5,1,-4,-1,6,0,-3,0],
			[4,1,-3,0,5,0,-3,0],
			[3,0,-2,0,4,0,-2,0],
			[2,0,-1,0,3,0,-1,0],
			[1,0,-1,0,2,0,-1,0]
		])
	elif F==F20:
		zn=F.gens()[0]
		B=[zn^7, -zn^6 - zn^5 + zn^3, zn^6, zn^5, -zn^5 - zn^4 - zn^3, zn^4 + zn^3 + zn^2, -zn^7 + zn^5 - zn^3 - zn^2 + zn, -zn^6 + zn^4 - zn^2 - zn + 1]
		M=matrix([
			[ 1 , 1 , 0 , 0 , 0 , 0 , 0 , 1],
			[ 2 , 1 , 0 , 1 ,-1 , 0 , 0 , 0],
			[ 3 , 1 , 0 , 1 ,-1 , 0 , 1 , 0],
			[ 4 , 1 , 1 , 1 ,-2 , 1 , 0 , 0],
			[ 3 , 1 , 1 , 0 ,-1 , 0 , 0 , 0],
			[ 2 , 1 , 1 , 0 , 0 , 0 , 0 , 0],
			[ 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0],
			[ 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0]
		])
	elif F==F24:
		zn=F.gens()[0]
		B=[zn^7, zn^4, zn^6, zn^5, -zn^7 - zn^6 - zn^5 + zn^3, zn^7 + zn^6 - zn^4 - 2*zn^3 - zn^2, -zn^6 - zn^5 + zn^3 + 2*zn^2 + zn, -zn^7 + zn^5 - zn^2 - zn + 1]
		M=matrix([
			[ 1 , 1 ,-1 , 1 , 0 , 0 , 0 , 1],
			[ 1 , 2 ,-1 , 0 , 1 , 0 , 0 , 0],
			[ 1 , 2 ,-1 , 1 , 0 , 0 , 1 , 0],
			[ 1 , 4 ,-2 , 1 , 0 , 1 , 0 , 0],
			[ 1 , 3 ,-2 , 1 , 0 , 0 , 0 , 0],
			[ 1 , 2 ,-1 , 0 , 0 , 0 , 0 , 0],
			[ 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0],
			[ 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0]		
		])
	elif F==F8: #only works when F=Q(sqrt(-1),sqrt(3),sqrt(5)) with generators given in this order
		i,s3,s5=F.gens()
		B=[(-1/4*s5 + 1/4)*i + (1/4*s5 - 1/4)*s3, ((1/4*s5 - 1/4)*s3 - 1)*i - 1/4*s5 + 1/4, (-1/4*s5 + 1/4)*i + (-1/4*s5 + 1/4)*s3, 1/2*s3*i + 1/2, -1/2*i - 1/2*s3 + 1/2*s5 - 1/2, ((-1/4*s5 - 1/4)*s3 + 1/4*s5 + 3/4)*i + (1/4*s5 + 3/4)*s3 - 3/4*s5 - 3/4, ((1/4*s5 - 1/4)*s3)*i + 1/4*s5 - 1/4, ((-1/4*s5 - 1/4)*s3 + 1/4*s5 + 1/4)*i + (-1/4*s5 - 1/4)*s3 + 1/4*s5 + 5/4]
		M=matrix([
			[ 1 , -1  ,  -7  ,  -17   ,-111  , -193  ,-1847 , -1553],
			[ 1 , -3  ,  -3  ,  -55   , -39  , -843  , -747 ,-10175],
			[ 1 , -3  ,  -7  ,  -51   ,-143  , -611  ,-2839 , -4883],
			[ 2 , -4  , -10  ,  -80   ,-182  ,-1084  ,-3202 ,-11336],
			[ 1 , -3  , -11  ,  -55   ,-151  , -811  ,-1955 , -9951],
			[ 1 , -3  ,  -7  ,  -51   , -79  , -611  , -855 , -4883],
			[ 1 , -1  ,  -3  ,  -37   ,  -7  , -521  ,  245 , -6061],
			[ 1 , -1  ,   1  ,  -17   ,  33  , -193  ,  353 , -1553]		
		])
	else: 
		raise ValueError("Field not supported")
	return B,M
		
	
def vector_in_E8basis(M,x): #takes the element x of F, computes its coordinates in B, and computes the corresponding vector with these coordinates in the E8 basis
	V=F.absolute_vector_space()[2]	
	w= M*V(x)
	return sum([w[i]*E[i] for i in range(8)])
	
def number_field_elt(B,v): #takes an element given by its vector of coordinates in the E8 basis and returns the corresponding number field element
	w=Einvmat*(v)
	return sum([w[i]*B[i] for i in range(len(v))])
		
def euclidean_div(F,a,b): #computes the euclidean division of a by b in F
	B,M=basis_E8(F)
	x=vector_in_E8basis(M,a/b)
	q=number_field_elt(B,closest_vector_E8(x))
	return q,a-b*q
		

