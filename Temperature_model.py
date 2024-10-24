# %%
import numpy as np
import csv
from mec3075f_p2_functions import *
import math as math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 


#%%
def getparameters(): 
    EMPL_ID = input("Enter peoplsoft code:\n") #entering people soft id

    data = [] 

    with open("allocation.csv") as csvfile: #opening file
        reader = csv.reader(csvfile)  #reading file
        for row in reader:    #for loop checks every row in the file, if starts with the correct ID then the data is appended to a list
            if row[0].startswith(EMPL_ID):
                data.append(row)
            
    return data                            
#%% 
parameters = getparameters()

#  getting list of parameters, accessing list made in get paramters but now i am creating  new a list
# setting a1.b1.a3.b2.b3 to the correct values so that they can be accessed later on 

a1 = float(parameters[0][1])
b1 = float(parameters[0][2])
a2 = float(parameters[0][3])
b2 = float(parameters[0][4])
a3 = float(parameters[0][5])
b3 = float(parameters[0][6])


parameters_list = [a1,b1,a2,b2,a3,b3]  


print(parameters_list)
# %%
m_ = int(input("enter the number of nodes along the width")) #input for mesh size, user can choose level of accuracy 
n_ = int(input("enter the number of nodes along the length"))

corners = np.array([[0, 0], [0, 0.044], [0.048, 0.060], [0.048, 0.044]])
node_list_, el_list_ = get_node_element_arrays(corners, m_, n_) #gives the coordinates of points of the nodes from the left to right 

print(node_list_)
print(el_list_)

# %%
order = int(input("what is the order of the system?")) 

def find_ecm_(elements, node_list_, el_list_, order):  #defining function ecm 
    k = np.zeros(shape=(4, 4))

    wvector = []
    phivector = []

    if order == 1:
        wvector = [2]
        phivector = [0]
    elif order == 2:
        wvector = [1, 1]
        phivector = [-1/math.sqrt(3), 1/math.sqrt(3)]
    elif order == 3:
        wvector = [5/9, 8/9, 5/9]
        phivector = [-math.sqrt(3/5), 0, math.sqrt(3/5)]

    elenod = el_list_[elements]
    coordinates = node_list_[elenod - 1]
    
    k = np.zeros(shape=(4,4))

    for i in range(order):
        for j in range(order):
            D = [[(1/4)*-(1-phivector[j]), (1/4)*(1-phivector[j]),(1/4)*(1+phivector[j]), -(1/4)*(1+phivector[j])],
                 [-(1/4)*(1-phivector[i]),-(1/4)*(1+phivector[i]),(1/4)*(1+phivector[i]),(1/4)*(1-phivector[i])]]
            J = np.matmul(D, coordinates)  #calculating jacobian 
            Determinantj = np.linalg.det(J) #calculating determinant of jacobian 
            Jinverse = np.linalg.inv(J) #calculating inverse of jacobian 
            B = np.matmul(Jinverse, D) #calculating B matrix 
            k = k+(wvector[i]*wvector[j]*np.transpose(B).dot(B)*Determinantj) #calculating K matrix 

    return k

GCM = np.zeros(shape=(5, 5))
GCM = assemble_gcm(node_list_, el_list_, find_ecm_)

print("GCM matrix is:", GCM)

#%% get constraints

def make_constrained_system(GCM, m_,n_, a, b): 
    product = m_*n_ 
    lefthandside_boundry = range(0, product, m_) #defining the changes on either side of the boundries respectiviley
    righthandside_boundry = range(m_ -1, product, m_)#defining the indicis from l&R boundaries
    hc=[] 
    delta_temperature = (b-a) / (n_ - 1) #defining the changes on either side of the boundries respectiviley
    
    tk_1 = np.zeros(shape =(product,1)) #temperature vector with constrants 
    
    for i in lefthandside_boundry:
        tk_1[i] = 100 
    for i in righthandside_boundry:
        tk_1[i] = a + (delta_temperature*((i//m_)-1))
        
    hc1 = -GCM @ tk_1 #heat matrix with constraints 
    
    
    hc = hc1.copy()
    
    for i in lefthandside_boundry: 
        hc[i] = 100 
    for i in righthandside_boundry:
        hc[i] =  a + (delta_temperature*((i//m_)-1))
        
    Kc = GCM.copy() #altered GCM matrix now due to constraints
    for i in range(product):
        if tk_1[i] !=0:
            Kc[i, :] = 0
            Kc[i, i] = 1
            
    return hc, Kc       

H1, GCM1  = make_constrained_system(GCM, m_, n_, a1,b1) #matrixes outputted must be symmetric for cholesky 

print("hc_1 is:\n",H1)

# for each set of values
# cholesky function

def cholesky(A):
    #compute matrix size of A
    len_A = len(A) 
    L = np.zeros(shape=(len_A, len_A)) #empty matrix L to be defined 
    if np.any(np.linalg.eigvals(A) <= 0):# Check if A is positive definite (if any eigenvalue is less than or equal to 0)
        A = A + 1e-8 * np.eye(len_A) 

    for x in range(len_A): # this is now the cholesky decomposition 
        u = A[x,x] - np.dot( L[x, 0:x] , L[x, 0:x] ) 
        L[x,x] = math.sqrt(u)
        for i in range(x+1, len_A):  
            L[i,x] = ( A[i,x] - np.dot(L[i,0:x], L[x,0:x]))/L[x,x]
    U = np.transpose(L) # U is calculated by transposing L matrix 

    return (L,U)

L, U = cholesky(GCM1)

# foward and back substitution 
def foward_substitution(L, b): 
    n = L.shape[0]
    y = np.linalg.solve(L, b)
    return y

def backward_substitution(U, y):
    n = U.shape[0]
    x = np.linalg.solve(U, y)
    return x

T1 = backward_substitution(U,foward_substitution(L, H1))


#plotting

xcord=[]
ycord=[]
xcord, ycord = node_list_.T
zcord = T1


fig = plt.figure(figsize=(11,15))
ax = plt.axes(projection='3d')

 
colourMap = plt.cm.jet
ax.scatter(xcord, ycord, zcord, c=T1, cmap=colourMap)
ax.scatter(xcord, ycord)

 

plt.title("Temperature Distribution of Systen")
plt.xlabel("m nodes")
plt.ylabel("n nodes")

ax.set_zlabel("T (degrees Celsius)")
plt.show()







