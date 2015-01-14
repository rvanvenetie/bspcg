import numpy as np
import sys
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
import cProfile
import plottri


'''
Assemble matrix with inproduct between al basis functions in the continous linear piecewise sheisse
'''
def AssembleM(N,T,G):
	EM = np.mat([[2,1,1],
							 [1,2,1],
							 [1,1,2]]) / 24.0
							 
	n = N.shape[0] #Aantal hoekpunten
	m = T.shape[0] #Aantal driehoeken
	
	M = np.zeros((n,n))
	for k in range(m):
		Tk = T[k,:]									 #Indices van hoekpunten die bij Tk horen
		C = N[Tk,:].T								 #Coordinaten van deze hoekpunten
		B = C[:,[1,2]] - C[:,[0, 0]] #Bepaal transformatiematrix van D(F_k^-1)
		#Now we can calculate all the nonzero inproducts on triangle k.
		M[np.ix_(Tk,Tk)] = M[np.ix_(Tk,Tk)] + abs(np.linalg.det(B)) * EM 
	return M

'''
Assemble matrix that holds all the inproducts between gradients of basis functions
'''
#@profile
def AssembleA(N,T,G):
	EA1 = np.mat([[-1,-1],[1,0],[0,1]])
	EA2 = np.mat([[-1,1,0],[-1,0,1]])
	n = N.shape[0] #Amount of vertices
	m = T.shape[0] #Amount of triangles
	A = np.zeros((n,n)) #Create matrix that is going to hold all the inproducts
	for k in range(m):
		Tk = T[k,:] #Vertices current triangle
		C = N[Tk,:].T #Get coordinates of the vertices
		B = np.mat(C[:,[1,2]] - C[:,[0, 0]]) #Assemble transformation matrix of D(F_k^-1)
		BI = np.linalg.inv(B) #Inverse of B (TODO hard code inverse?)
		#
		Ak = 0.5 * abs(np.linalg.det(B)) * EA1 * BI * BI.T * EA2 #See theory
		#print Tk
		#print Ak
		#D = np.mat(C[:,[2,0,1]] - C[:,[1,2,0]])
		#Ak = D.T * D / (2 * abs(np.linalg.det(B)))

		#Now we can calculate all the nonzero inproducts on triangle k.
		A[np.ix_(Tk,Tk)] = A[np.ix_(Tk,Tk)] + Ak
	return A
	

def Norms(A,M,v):
	v = np.mat(v).T
	return (np.sqrt(v.T*(A+M)*v), np.sqrt(v.T*M*v))
	
'''
Apply FEM on f
'''
#@profile
def LinFEM(N,T,G,f,g = 1):
	A = AssembleA(N,T,G)
	M = AssembleM(N,T,G)
	n = N.shape[0] #Amount of vertices

	GInd = np.where(G == 0)[0] #Find the indices of basis functions that are not on the border
	GIx = np.ix_(GInd,GInd) #Indices of all 
	#Remove basis functions on edge
	Ag = A[GIx]			 #Stiffness matrix
	Mg = g * M[GIx] #Mass matrix times reaction term g
	Mf	= M[GInd,:]  #Matrix that holds the inner products for the calculation of (f,\psi_i)
	
	#Bereken de waarde van f in de punten van onze partitie
	F = map(f,N)
	F = np.mat(F).T


	#Finally solve for uh!
	sol = np.linalg.solve(Ag + Mg,Mf*F)
	#print Ag
	print Mf*F
	#print sol
	#print "HOERENZOON"
	uh = np.zeros(n) #Values of the solution at vertices
	uh[GInd] = sol
	
	return uh, sol, A, M
	
	
	
'''
TODO: Replace all appends with predefined matrices
'''
def Refine(N,T,G):
	n = N.shape[0] #Amount of vertices
	m = T.shape[0] #Amount of triangles
	refined = np.zeros((n,n)) #Indicates if the edge has a new refined vertice yet
	edges = [(0,1),(0,2),(1,2)] #All edges in a triangle
	TT = np.zeros((m*4,3),dtype=np.int)
	TIndex = 0 #Current index in TT
	for k in range(m):
		Tk = T[k,:] #Vertices
		C = N[Tk,:] #Get coordinates of the vertices
		new_indices = np.zeros(3,dtype=np.int) #Store indices of new vertices in triangle
		for j in range(3): #Find indices of vertices of sub-triangles
			a,b = edges[j] #One of the three sides of triangle Tk
			if refined[Tk[a],Tk[b]] == 0: #We have not yet created new vert on this edge
				#print "Creating midpoint on", Tk[a],Tk[b]
				new_vert =	(C[a,:] + C[b,:])/ 2.0
				N = np.append(N,[new_vert],axis=0)
				G = np.append(G,1) #Assume this vertice is on the boundary
				new_indices[j] = n
				refined[Tk[a],Tk[b]] = refined[Tk[a],Tk[b]] = n
				n = n +1
			else:
				new_indices[j] = refined[Tk[a],Tk[b]]
				G[new_indices[j]] = 0 #Not on the boundary
				#print "Already have a midpoint on", Tk[a],Tk[b]
		#New_indices now holds the vertices of the sub-triangles
		TT[k*4:4*(k+1),:] = [[Tk[0],new_indices[0],new_indices[1]],
										 [Tk[1],new_indices[0],new_indices[2]],
										 [Tk[2],new_indices[1],new_indices[2]],
										 [new_indices[0],new_indices[1],new_indices[2]]]
	return N,TT,G
	
def Voorbeeld():
	'''Exacte oplossing'''
	fig = plt.figure()
	x= np.linspace(0, 1, 50)
	y = np.linspace(0, 1, 50)
	X,Y = np.meshgrid(x, y) # grid of point
	u  = np.sin(np.pi * X) * np.sin(np.pi* Y)
	ax = fig.gca(projection='3d')
	ax.plot_surface(X,Y,u,rstride=2, cstride=2,		linewidth = 0.1, cmap=cm.RdBu)

	plt.show()
	'''Unit SQUARE'''
	N = np.array([[0,0],
								[0,1],
								[1,1],
								[1,0]])
	T = np.array([[0,1,2],
								[0,2,3]])
	G = np.array([1,1,1,1])
	j = 5
	fig = plt.figure()
	normsa = np.zeros(j)
	norms2 = np.zeros(j)
	for i in range(j):
		N,T,G = Refine(N,T,G)
			 
		u  = np.sin(np.pi * N[:,0]) * np.sin(np.pi* N[:,1])
		#f = np.mat(u * (2 * np.pi**2 + 1)).T
		f = lambda pt: np.sin(np.pi * pt[0]) * np.sin(np.pi * pt[1]) *(2 * np.pi**2 + 1)
		g = 1
		#f = lambda pt: g * -(2 * pt[0] - 1)**2 * (2 * pt[1] - 1)**2 
		
		uh, A, M = LinFEM(N,T,G,f,g)
		normsa[i],norms2[i] = Norms(A,M,u-uh)
		
		
		#fig = plt.figure()
		#ax = fig.add_subplot(111)
		ax = fig.add_subplot(2,j,j+i+1)
		ax.triplot(N[:,0],N[:,1],T)
		ax.set_title('Triangulatie')
		#fig = plt.figure()
		#ax = fig.add_subplot(111,projection='3d')
		ax = fig.add_subplot(2,j,i+1, projection='3d')
		ax.plot_trisurf(N[:,0],N[:,1],uh,		 cmap=cm.RdBu)
		ax.set_title('$u_h$')
		#ax = fig.add_subplot(j,3,i*3+3, projection='3d')
		#ax.plot_trisurf(N[:,0],N[:,1],u-uh)
		#ax.set_title('$u - u_h$')
		'''
		plt.figure()
		ax = fig.add_subplot(2,j,j+i+1)
		plt.triplot(N[:,0],N[:,1],T)
		plt.title('Triangulatie')
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.plot_trisurf(N[:,0],N[:,1],uh)
		ax.set_title('$u_h$')
		'''


	ax = plt.figure().add_subplot(111)
	#plt.yscale('log')
	plt.plot(normsa)
	plt.plot(norms2)
	#print normsa.shape
	plt.legend([r"$||\Pi u - u_V^\Pi||_A$",r"$||\Pi u - u_V^\Pi ||_2$"])
	plt.xticks(range(j))
	ax.set_xticklabels(range(1,j+1))
	plt.show()		
	
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#Voorbeeld()
fn = sys.argv[1] if len( sys.argv) > 1 else "lshaped.m"
(distributed, nverts, ntris, x, y, b, t, P, p) = plottri.readmesh( fn)
N = np.array([x, y]).T
T = t.astype('int')
G = b.astype('int')
	
for i in range(10):
    print "Refine"
    N,T,G = Refine(N,T,G)
    plottri.writemesh( "faka%d.m" % i, N, T, G)
    print G.shape

f = lambda pt: 1
g = 0
uh, sol, A, M = LinFEM(N,T,G,f,g)
for i in range( len( sol)):
    print "%.5f" % sol[i]
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.triplot(N[:,0],N[:,1],T)
ax.set_title('Triangulatie')
ax = fig.add_subplot(2,1,2, projection='3d')
ax.plot_trisurf(N[:,0],N[:,1],uh,		 cmap=cm.RdBu)
ax.set_title('$u_h$')
plt.show()
