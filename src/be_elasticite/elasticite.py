import numpy as np 
import matplotlib.pyplot as plt 
import scipy.io 
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D



#%% PARAMÈTRES DÉTERMINANTS LE CALCUL
################################################

maille = 3  # Choix du maillage 1 (gros),2 (plus fin) ,3 (encore plus fin) ou 4 (trou)

# Proprietes mecaniques
E = 100 # Young's modulus  
v = 0.3  # Poisson's coefficient

# Forces (densité surfacique ici)
fvx = 100          
fvy = 0 
fv = np.array([fvx,fvy])



#%% LECTURE ET CHARGEMENT DU MAILAGE
#######################################################

if maille == 1:
    mat = scipy.io.loadmat('./maillage_1.mat')
elif maille ==2:
    mat = scipy.io.loadmat('./maillage_2.mat')
elif maille ==3:
    mat = scipy.io.loadmat('./maillage_3.mat')
elif maille ==4:
    mat = scipy.io.loadmat('./maillage_trou.mat')
else:
    raise ValueError('Maille 1,2,3 ou 4 ')
    
# paramètres déduits
#===============================================

# coordonnes des noeuds  
X = mat['p'].T  

# table de connectivite
C = mat['t'].T[:,:-1] -1   

# nombre de noeuds du maillage
n_nodes = X.shape[0] 

# nombre d'elements du maillage
n_elems = C.shape[0] 

# Nombre de noeuds par element
n_nodes_elem = 3 # ce sont des triangles P1

# Nombre de ddl par element
ndofe= 2*n_nodes_elem # champ vectoriel en 2D

# Nombre de ddl dans tout le domaine  
ndof = 2*n_nodes # champ vectoriel en 2D


# dimensions de la plaque 
xmin =  np.min(X[:,0])
xmax =  np.max(X[:,0])
ymin =  np.min(X[:,1])
ymax =  np.max(X[:,1])
# vérifier dimensions de la plaque [0,2]x[0,1]

# affichage du maillage
Initial_triangulation =  mtri.Triangulation(X[:,0],X[:,1],C)
plt.figure(0)
plt.triplot(Initial_triangulation, color = 'black')
plt.title("maillage")
plt.show()


#%% CALCUL DES FUNS EF P1 LAGRANGE
############################################
 
def fun_tri_P1_lag(x,y,x_nodes,y_nodes):
    # x,y -> point d'évaluation
    # x_nodes -> tableau 1D numpy avec abscisses des noeuds
    # y_nodes -> tableau 1D numpy avec ordonnées des noeuds
    
    Te = 1/2*np.abs((x_nodes[1]-x_nodes[0])*(y_nodes[2]-y_nodes[0])
    - (y_nodes[1]-y_nodes[0])*(x_nodes[2]-x_nodes[0]))
    
    N1 = 1/(2*Te)*(
        x_nodes[1]*y_nodes[2] 
        - x_nodes[2]*y_nodes[1]
        + (y_nodes[1]-y_nodes[2])*x
        + (x_nodes[2]-x_nodes[1])*y
        )
    N2 =  1/(2*Te)*(
        x_nodes[2]*y_nodes[0] 
        - y_nodes[2]*x_nodes[0]
        + (y_nodes[2]-y_nodes[0])*x
        + (x_nodes[0]-x_nodes[2])*y
        )
    N3 =  1/(2*Te)*(
        x_nodes[0]*y_nodes[1] 
        - y_nodes[0]*x_nodes[1]
        + (y_nodes[0]-y_nodes[1])*x
        + (x_nodes[1]-x_nodes[0])*y
        )
    N = np.array([N1,N2,N3])
    dNdx = np.array([
        (y_nodes[1]-y_nodes[2])/(2*Te),
        (y_nodes[2]-y_nodes[0])/(2*Te),
        (y_nodes[0]-y_nodes[1])/(2*Te)
        ])
    
    dNdy = np.array([
        (x_nodes[2]-x_nodes[1])/(2*Te),
        (x_nodes[0]-x_nodes[2])/(2*Te),
        (x_nodes[1]-x_nodes[0])/(2*Te)
        ])
    return [N, dNdx, dNdy]
    # N = [N1(x,y),N2(..),N3(..)]-> tableau 1D numpy avec valeur des 3 funs au point d'éval
    # dNdx -> idem avec dérivée par rapport à x des funs
    # dNdy -> idem avec dérivée par rapport à y des funs
    
def GetBe(dNdx,dNdy):
    # dNdx=[dN1dx(x,y),dN2dx(x,y),dN3dx(x,y)] -> tableau 1D numpy avec valeur à un certain point d'éval
    
    Be = np.zeros((3,6))
    Be[0,:3] = dNdx
    Be[1,3:] = dNdy
    Be[2] = np.hstack((dNdy,dNdx))
    
    return Be
    # Be : matrice élémentaire lien déformation - déplacement (evaluee à un certain point)


print(GetBe(np.array([1,2,3]), np.array([11,12,13])))

def GetNe(N):
    Ne_matrix = np.zeros((2,6))
    Ne_matrix[0,:3] = N
    Ne_matrix[1,3:] = N
    
    return Ne_matrix
    # Ne_matrix : matrice élémentaire lien déplacement - ddls 
    
print(GetNe(np.array([1,2,3])))

#%%
#--------------------------------------------------------------------------
#
# 1. CALCUL DES DEPLACEMENTS
#
#--------------------------------------------------------------------------
    

#%% CALCUL DE LA MATRICE DE RIGIDITE
#######################################################
    

# Loi de comportement : Contraintes - Deformations
# Hypothese de contraintes planes 
# convention : [sigmaxx,sigmayy,sigmaxy] = H [epsxx, epsyy, 2*epsxy] 
H = np.array([[E/(1-v**2), v*E/(1-v**2), 0],
     [v*E/(1-v**2), E/(1-v**2), 0],
     [0, 0, E/(2*(1+v))]])

    
# initialisation de la matrice de rigiditié 
K = np.zeros((ndof, ndof))

# Boucle sur les éléments

for e in range(n_elems):
    
    # coordonnes des noeuds de l'element e   
    x_nodes = np.zeros(3)
    y_nodes = np.zeros(3)
    
    C_e = C[e]
    for i in range(3):
        x_nodes[i], y_nodes[i] = X[C_e[i]]

    # calcul de la matrice de rigidité elementaire
    
    # A REMPLIR
    # POUR CELA, AVANT, REMPLIR LES FUNS fun_tri_P1_lag ET GetBe CI-DESSUS
    N, dNdx, dNdy = fun_tri_P1_lag(x_nodes[0], y_nodes[0],x_nodes,y_nodes)
    Be = GetBe(dNdx, dNdy)
    
    Te = 1/2*np.abs((x_nodes[1]-x_nodes[0])*(y_nodes[2]-y_nodes[0])
    - (y_nodes[1]-y_nodes[0])*(x_nodes[2]-x_nodes[0]))
    
    Ke = Be.T @ H @ Be
    Ke = np.multiply(Te, Ke)
    
    # Assemblage des contributions élémentaires
    for i in range(n_nodes_elem):
        for j in range(n_nodes_elem):
            K[C_e[i], C_e[j]] += Ke[i,j]
            K[C_e[i]+n_nodes, C_e[j]] += Ke[i+n_nodes_elem,j]
            K[C_e[i], C_e[j]+n_nodes] += Ke[i,j+n_nodes_elem]
            K[C_e[i]+n_nodes, C_e[j]+n_nodes] += Ke[i+n_nodes_elem,j+n_nodes_elem]
            


#%% CALCUL DU SECOND MEMBRE
#######################################################
           
# initialisation du second membre
F = np.zeros(ndof)

# Boucle sur les éléments

for e in range(n_elems):
    
    # coordonnes des noeuds de l'element e   
    x_nodes = np.zeros(3)
    y_nodes = np.zeros(3)
    
    C_e = C[e]
    for i in range(3):
        x_nodes[i], y_nodes[i] = X[C_e[i]]
    
    # POUR CELA, AVANT, REMPLIR LA FUN GetNe CI-DESSUS    
    N1, _, _ = fun_tri_P1_lag(x_nodes[0], y_nodes[0],x_nodes,y_nodes)
    N2, _, _ = fun_tri_P1_lag(x_nodes[1], y_nodes[1],x_nodes,y_nodes)
    N3, _, _ = fun_tri_P1_lag(x_nodes[2], y_nodes[2],x_nodes,y_nodes)
    Ne1 = GetNe(N1)
    Ne2 = GetNe(N2)
    Ne3 = GetNe(N3)
    
    Te = 1/2*np.abs((x_nodes[1]-x_nodes[0])*(y_nodes[2]-y_nodes[0])
    - (y_nodes[1]-y_nodes[0])*(x_nodes[2]-x_nodes[0]))

    Fe = Te/3 * (Ne1.T @ fv + Ne2.T @ fv + Ne3.T @ fv)
    
    for i in range(n_nodes_elem):
        F[C_e[i]] += Fe[i]
        F[C_e[i] + n_nodes] += Fe[i+n_nodes_elem]

#%% IMPOSITION DES CL DE DIRICHLET
#######################################################
        
        
for n in range(n_nodes): # boucle sur les noeuds
    

    if X[n,0] == 0: # noeuds bloqués
    # On bloque 
        K[n,:] = 0
        K[:,n] = 0
        F[n] = 0
        K[n,n] = 1 # Sinon pas inversible !
        
        K[n+n_nodes] = 0
        K[:,n+n_nodes] = 0
        F[n+n_nodes] = 0
        K[n+n_nodes, n+n_nodes] = 1 # Sinon pas inversible

    
#%% RESOLUTION DU SYSTÈME LINÉAIRE ET VISUALISATION
#######################################################
        
# résolution
U = np.linalg.solve(K,F)

#Calcul des coordonnes des noeuds apres deformation
x = X[:,0] + U[:n_nodes]
y = X[:,1] + U[n_nodes:]

# affichage du maillage
Deformed_triangulation  =  mtri.Triangulation(x,y,C)
plt.figure(1)
plt.triplot(Initial_triangulation, color = 'black')
plt.triplot(Deformed_triangulation, color = 'blue')
plt.show()
 

#%%
#--------------------------------------------------------------------------
#
# 2. POST-TRAITEMENT : CALCUL DES CONTRAINTES DE VON MISES
#
#--------------------------------------------------------------------------
    

# initialisations

T = np.zeros(n_nodes)       # Surface
SVM = np.zeros(n_nodes)     # contrainte de Von Mises


#boucle sur les elements 

for e in range(n_elems):
    
    # coordonnes des noeuds de l'element e   
    x_nodes = X[C[e,:],0]
    y_nodes = X[C[e,:],1] 
    
    x1=x_nodes[0]; x2=x_nodes[1]; x3=x_nodes[2] 
    y1=y_nodes[0]; y2=y_nodes[1]; y3=y_nodes[2] 
    
    # surface de l'élément e
    Te = 0.5*np.abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
    
    # Calcul des contraintes de Von Mises
    [Nm,dNmdx,dNmdy] = fun_tri_P1_lag(x1,y1,x_nodes,y_nodes)
    Be = GetBe(dNmdx,dNmdy)
    Ue = np.r_[U[C[e,:]],U[C[e,:]+n_nodes]]
    sigma = H.dot(Be.dot(Ue))
    sxx = sigma[0] 
    syy = sigma[1]  
    sxy = sigma[2] 
    svm = (sxx**2+syy**2+3*sxy**2-sxx*syy)**0.5 
    
    # Assemblage : surface et sigma
    for i in range(n_nodes_elem):
        T[C[e,i]] +=  Te
        SVM[C[e,i]] +=  Te*svm

# Calcul de la contrainte de Von Mises aux noeuds
SVM = SVM/T

# Visualisation des containtes 
t = mtri.Triangulation(x,y,C)
plt.figure(2)
plt.triplot(t)
plt.tricontourf(t,SVM)
plt.colorbar()
plt.title('Contraintes de Von Mises')
 