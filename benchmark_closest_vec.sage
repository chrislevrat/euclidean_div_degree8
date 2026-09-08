F=F15


def reflexion(x,y): #s_y(x), s_y is the reflexion about the hyperplane orthogonal to y, when ||y||²=2
	if x==0:
		return vector([0,0,0,0,0,0,0,0])
	else:
		return x-(x*y)*y

E8=E
for e in E:
	for f in E:
		if reflexion(e,f) not in E8:
			E8.append(reflexion(e,f))
		
		
maxroot=vector([0,0,0,0,0,0,1,1])
pt=[
vector([0,0,0,0,0,0,0,0]),
vector([0,1/2,0,0,0,0,0,0]),
vector([0,0,0,0,0,0,1/2,1/2]),
vector([0,0,0,0,0,0,0,1]),
vector([0,0,0,0,0,1/3,1/3,2/3])
]

maxlist=[]


def rand_algint(B,u):
	L=[randint(-u,u)*b for b in B]
	return sum(L)	
		

def norm_polynomial(F): #for cyclotomic fields
	B,M=basis_E8(F)
	ConjList=[]
	R=PolynomialRing(QQ,8,'t')
	t=R.gens()
	for phi in F.automorphisms():
		ConjList.append(sum([t[i]*phi(B[i]) for i in range(8)]))
	t=var('t0 t1 t2 t3 t4 t5 t6 t7')
	return t0.parent(product(ConjList)) 


t=var('t0 t1 t2 t3 t4 t5 t6 t7')
t=vector([t0,t1,t2,t3,t4,t5,t6,t7])

NF=norm_polynomial(F)

#for i in range(240):
#	r=E8[i]
#	consr=[sum([t[i]*(E[i]*reflexion(E[j],r)) for i in range(8)]) for j in range(8)] #contraintes pour être 
#	consr.append(1-sum([t[i]*(E[i]*reflexion(maxroot,r)) for i in range(8)]))
#	M=0
#	for k in range(5):
#		ur=minimize_constrained(1/NF,consr,list(reflexion(pt[k],r)))
#		m=NF.subs(t0=ur[0],t1=ur[1],t2=ur[2],t3=ur[3],t4=ur[4],t5=ur[5],t6=ur[6],t7=ur[7])
#		M=max(m,M)
#	maxlist.append(M)
#	if i%10==0:
#		print(i)
		
		
		
		
# pour F=F15=Q(zeta15) E8[183] donne uex=(3/2,3/2,5/2,3,5/2,3/2,1,1/2) pour lequel on a absolute_norm=61/256


#u1=vector([0,0,0,0,0,0,0,1])
#u2=(1/6)*vector([1,1,1,1,1,1,1,5])

#VoronoiVertices=[u1,u2]
#for r in E8:
#	for ui in VoronoiVertices:
#		vi=reflexion(ui,r)
#		if vi not in VoronoiVertices:
#			VoronoiVertices.append(vi)
			
		

# pour le corps F8:
L8=F8.absolute_field('c')
B8,M8=basis_E8(F8)
B8=[L8(b) for b in B8]
ConjList=[]
R=PolynomialRing(QQ,8,'t')
t=R.gens()
for phi in L8.automorphisms():
	ConjList.append(sum([t[i]*phi(B8[i]) for i in range(8)]))
t=var('t0 t1 t2 t3 t4 t5 t6 t7')
NF8= t0.parent(R(product(ConjList))) 
t=vector([t0,t1,t2,t3,t4,t5,t6,t7])
#pour F8 il suffit ci-dessous de prendre PEF[0] : on y trouve déjà le maximum


PE=load("VoronoiPolyhedronInE8Coordinates.sobj") #plus facile de charger cet objet: la région de Voronoi autour de 0 exprimée en coordonnées dans la base E
PEF=PE.facets()
PEFacesEqs=[]
for PEFi in PEF:
     U=PEFi.as_polyhedron()
     Ui=U.inequalities()
     Ue=U.equations()
     L=[]
     for u in Ui:
         L.append(u.A()*t+u.b()) #les inégalités exprimées en fonction de t
     for u in Ue:
         L.append(u.A()*t+u.b()) #les égalités sont données comme deux inégalités
         L.append(-(u.A()*t+u.b()))
     PEFacesEqs.append(L)
     
     
ListMaxOpt=[]     
     
for i in range(len(PEF)):
	U=PEF[i].as_polyhedron()
	pt=list(U.center())
	ur=minimize_constrained(1/NF,PEFacesEqs[i],pt)
	uval=NF.subs(t0=ur[0],t1=ur[1],t2=ur[2],t3=ur[3],t4=ur[4],t5=ur[5],t6=ur[6],t7=ur[7])
	ListMaxOpt.append(uval)
	if i%10==0:
		print(i)

#fonction d'arrondi des coordonnées dans la base [1,zeta,zeta^2,...]

def stupid_round(F,x):
	zn=F.gens()[0]
	v=F.absolute_vector_space()[2](x)
	return sum([round(v[i])*zn^i for i in range(8)])

def eucl_min(F):
	if F==F15:
		return 1/16
	if F==F20:
		return 1/5
	if F==F24:
		return 1/4

#compare l'arrondi stupide à la méthode par recherche de plus proche vecteur dans E8
#sur nb éléments tirés au hasard avec |coeffs|<=u, compte le nombre de fois que chaque méthode est en-dessous du minimum euclidien et le nombre de fois que la méthode par arrondi donne le résultat le plus petit des deux
def benchmark_compare(F,nb,u):
	count_E8_less_minimum = 0
	count_stupid_less_minimum = 0
	count_E8_less_stupid = 0
	count_E8_more_stupid = 0
	count_E8_equal_stupid = 0
	count_opt_best = 0
	m = eucl_min(F)
	B,M=basis_E8(F)
	for j in range(nb):
		x=rand_algint(B,u)/randint(1,u)
		y_stupid=stupid_round(F,x)
		y_E8=number_field_elt(B,closest_vector_E8(vector_in_E8basis(F,M,x)))
		_,y_opt,_=closest_number_field_elt(F,x)
		n_opt=y_opt.absolute_norm()
		n_stupid=(x-y_stupid).absolute_norm()
		n_E8=(x-y_E8).absolute_norm()
		if ((n_opt<n_E8) & (n_opt<n_stupid)):
			count_opt_best+=1
		if n_E8<=m:
			count_E8_less_minimum += 1
		if n_stupid<=m:
			count_stupid_less_minimum += 1
		if n_E8 < n_stupid:
			count_E8_less_stupid += 1
		elif n_E8 > n_stupid:
			count_E8_more_stupid += 1
		elif n_E8 == n_stupid:
			count_E8_equal_stupid += 1
	return count_E8_less_minimum, count_stupid_less_minimum, count_E8_less_stupid, count_E8_equal_stupid, count_E8_more_stupid


		



