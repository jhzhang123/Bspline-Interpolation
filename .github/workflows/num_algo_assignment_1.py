#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate

# In[8]:

input1 = np.array([(0,0),(1,0),(1,1),(0,1),(0,2)])
input2 = np.array([(0,0),(1,0),(1,1),(0,1),(0,2),(1,3),(15,16)])
input3 = np.array([(0,5),(4,9),(4,4),(4,6.5),(5,6.7),(4,5.7),(5,8),(5,3),(0,1),(7,4),(0,0),(1,1),(2,2)])
input4 = np.array([(0,0),(1,2),(4,8),(3,6),(2,4),(1,2),(2,3.2),(3,3.8),(4.3,3.7),(3,2.3),(2,2),(1,2),(2,1.8),(3,0.8),(4,0),(5,0.8)])
input5 = np.array([(0,5),(4.2,9.2),(4,9),(4,1.2),(4.1,4.9),(4.5,6),(4.4,5.9),(4,3.8),(4.5,7),(4.5,6),(4.5,2),                    (0,-1),(0.1,-0.9),(4.5,2),(8,4),(7.2,3.2),(4.5,2),                    (2,0),(2.1,0.1),(4,1),(2,-2),(2,-0.2),(4.5,-2),                    (4.5,1),(6.2,1),(5,-1.8),(4.5,0.6),(7,-2),(7.5,-1)])
input6 = np.array([(1,7.5), (1.75, 7.5), (0,5),(4,9), (2,0), (4,8), (8,8), (8,6), (4,6), (3,7), (5,7), (6,8), (3,0)])

input7 = []

u = np.linspace(0,1,num=10)

for u in u:
    x = np.exp(np.cos(6.2*u+51/30)+0.1)*np.cos(12.4*u)
    y = np.exp(np.cos(6.2*u+51/30)+0.1)*np.sin(12.4*u)
    input7.append((x,y))
print(input7)

#input_file = open('input.txt','r').readlines()
#plist = list()
#for line in input_file:
#    plist.append([int(i) for i in line.split()])

plist = np.array(input7)

#print(plist)
#plist = input5
#n : No. of data points - 1
#n_d : No. of De Boor (control point)
#k : Degree
n = len(plist) - 1

k = 3

# (k+1) + n_d = knots = n + 1 + k * 2
n_d = n+1+2*k-k-1

#print(n_d,k,n)

#Step 1 : Parameterization

D = np.zeros((n+1))
#print('D: ',D)
#print(plist[1][1])
for i in range(n):
    D[i] = math.hypot(plist[i+1][0] - plist[i][0], plist[i+1][1] - plist[i][1])
#print(D)
#print(np.sum(D))
t = np.linspace(0,1,n+1,endpoint=True)

#apply chord length parameterizaion
for i in range(1,n):
    t[i] = t[i-1] + D[i-1]/np.sum(D)
#print(t)   

u = np.append([0]*k,t)
u = np.append(u,[1]*k)

#print('t vector: ',t)
#print('u vector: ',u)


# In[4]:


N = np.zeros([n+1,n+3])
#print(N.shape)
def cubic_basis(u,t,i):
    if ( (u[i+1]-u[i])*(u[i+2]-u[i])*(u[i+3]-u[i]) != 0 and t>=u[i] and t<=u[i+1] ):
#         if ((u[i+1]-u[i])*(u[i+2]-u[i])*(u[i+3]-u[i]) == 0):
#             return 0
        #print("interval 1")
        return ((t-u[i])**3) / ((u[i+1]-u[i])*(u[i+2]-u[i])*(u[i+3]-u[i]))
    
    elif ( (u[i+2]-u[i+1])*(u[i+3]-u[i])*(u[i+2]-u[i])*(u[i+3]-u[i+1])*(u[i+4]-u[i+1]) != 0 and t>=u[i+1] and t<=u[i+2]):
#         if ((u[i+2]-u[i+1])*(u[i+3]-u[i])*(u[i+2]-u[i])*(u[i+3]-u[i+1])*(u[i+4]-u[i+1]) == 0):
#             return 0
        #print("interval 2")
        return ((t-u[i])**2)*(u[i+2]-t)/((u[i+2]-u[i+1])*(u[i+3]-u[i])*(u[i+2]-u[i])) +                (u[i+3]-t)*(t-u[i])*(t-u[i+1])/((u[i+2]-u[i+1])*(u[i+3]-u[i+1])*(u[i+3]-u[i])) +                (u[i+4]-t)*((t-u[i+1])**2)/((u[i+2]-u[i+1])*(u[i+4]-u[i+1])*(u[i+3]-u[i+1]))
    
    elif ( ((u[i+3]-u[i+2])*(u[i+3]-u[i+1])*(u[i+3]-u[i])*(u[i+4]-u[i+1])*(u[i+4]-u[i+2])) != 0 and t>=u[i+2] and t<=u[i+3]):
#         if (((u[i+3]-u[i+2])*(u[i+3]-u[i+1])*(u[i+3]-u[i])*(u[i+4]-u[i+1])*(u[i+4]-u[i+2])) == 0):
#             return 0
        #print("interval 3")
        return (t-u[i])*(u[i+3]-t)**2/((u[i+3]-u[i+2])*(u[i+3]-u[i+1])*(u[i+3]-u[i])) +                (u[i+4]-t)*(u[i+3]-t)*(t-u[i+1])/((u[i+3]-u[i+2])*(u[i+4]-u[i+1])*(u[i+3]-u[i+1])) +                (u[i+4]-t)**2*(t-u[i+2])/((u[i+3]-u[i+2])*(u[i+4]-u[i+2])*(u[i+4]-u[i+1]))
    
    elif ( (u[i+4]-u[i+3])*(u[i+4]-u[i+2])*(u[i+4]-u[i+1]) !=0 and t>=u[i+3] and t<=u[i+4]):
        #print("interval 4")
        return (u[i+4]-t)**3/((u[i+4]-u[i+3])*(u[i+4]-u[i+2])*(u[i+4]-u[i+1]))
    
    else:
        return 0
    
for i in range(1,n):
    #print(i)
    #N_k(tk), N_k+1(tk), N_k+2(tk)
    N[i][i] = cubic_basis(u,t[i],i)
    N[i][i+1] = cubic_basis(u,t[i],i+1)
    N[i][i+2] = cubic_basis(u,t[i],i+2)
    

print(N)
O = N
#add endpoint conditions
#when t0 = u3, N[0][0] = (N0(t0)'') = 6(u4-u)/[(u4-u3)(u4-u2)(u4-u1)]

N[0][0] = 6*(u[4]-t[0])/((u[4]-u[3])*(u[4]-u[2])*(u[4]-u[1]))
N[0][1] = (6*u[3]-2*u[1]-4*u[4])/((u[4]-u[3])*(u[4]-u[2])*(u[4]-u[1])) +           (6*u[3]-2*u[2]-2*u[5]-2*u[4])/((u[4]-u[3])*(u[5]-u[2])*(u[4]-u[2])) +           (6*u[3]-4*u[5]-2*u[3])/((u[4]-u[3])*(u[5]-u[3])*(u[5]-u[2]))
N[0][2] = (4*u[2]+2*u[4]-6*u[3])/((u[4]-u[3])*(u[5]-u[2])*(u[4]-u[2])) +           (2*u[2]+2*u[3]+2*u[5]-6*u[3])/((u[4]-u[3])*(u[5]-u[3])*(u[5]-u[2])) +           (4*u[3]+2*u[6]-6*u[2])/((u[4]-u[3])*(u[6]-u[3])*(u[5]-u[3]))
#when tn= u(n+3), N[n][n]
#print(n)
N[n][n] = (6*t[n]-2*u[n]-4*u[n+3])/((u[n+3]-u[n+2])*(u[n+3]-u[n+1])*(u[n+3]-u[n])) +            (6*t[n]-2*u[n+1]-2*u[n+3]-2*u[n+4])/((u[n+3]-u[n+2])*(u[n+4]-u[n+1])*(u[n+3]-u[n+1])) +           (6*t[n]-4*u[n+4]-2*u[n+2])/((u[n+3]-u[n+2])*(u[n+4]-u[n+2])*(u[n+4]-u[n+1]))

N[n][n+1] = (4*u[n+1]+2*u[n+3]-6*t[n])/((u[n+3]-u[n+2])*(u[n+4]-u[n+1])*(u[n+3]-u[n+1])) +           (2*u[n+1]+2*u[n+2]+2*u[n+4]-6*t[n])/((u[n+3]-u[n+2])*(u[n+4]-u[n+2])*(u[n+4]-u[n+1])) +           (2*u[n+5]+4*u[n+2]-6*t[n])/((u[n+3]-u[n+2])*(u[n+5]-u[n+2])*(u[n+4]-u[n+2]))

N[n][n+2] = 6*(t[n]-u[n+2])/((u[n+3]-u[n+2])*(u[n+4]-u[n+2])*(u[n+5]-u[n+2]))


#P0 interpolate D0 ==> N0(t0) = 1
N = np.insert(N,0,np.zeros(n+3),axis = 0)
N = np.insert(N,n+2,np.zeros(n+3),axis = 0)
N[0][0] = 1
N[n+2][n+2] = 1

#print(N)


# In[5]:


P = np.zeros((n+3,2))
#print(P)
plist_r = np.insert(plist,1,[0,0],axis=0)
plist_r = np.insert(plist_r,n+1,[0,0],axis=0)
#plist_r = plist_r.astype(np.float)
#print(plist_r)

#print(plist)

P = np.linalg.solve(N, plist_r)
#np.around(P,decimals=2)
print(P)


# In[6]:


def r(tau):
    sum = 0;
    for i in range(n+3):
        #print(cubic_basis(u,tau,i))
        sum += P[i]*cubic_basis(u,tau,i)
    return sum

u3=np.linspace(0,1,100,endpoint=True)

B_spline_dot = []
for j in range(len(u3)):
    B_spline_dot.append(r(u3[j]))

B_spline_dot = np.array(B_spline_dot)

#print(B_spline_dot)
x_i,y_i=plist.T
plt.plot(B_spline_dot[:,0],B_spline_dot[:,1],'b',linewidth=2.0,label='B-spline curve')
plt.scatter(x_i,y_i)

x,y=P.T
plt.plot(x,y,'k--',label='Control polygon',marker='o',markerfacecolor='red')

plt.legend(loc='best')
plt.gca().set_aspect('equal')
plt.axis([min(x)-0.5, max(x)+0.5, min(y)-0.5, max(y)+0.5])
plt.rcParams["figure.figsize"] = [9,9]

plt.title('Cubic B-spline curve evaluation')
plt.show()


# In[7]:


output = [k,'\n',n_d,'\n',u,'\n',P]
#print(output)

output_file = open('output.txt','w')
for line in output:
    output_file.write(str(line))
output_file.close()


# In[285]:


# Back up solution to draw spline using scipy package. 

# #Control Points
# x,y=P.T

# #Interpolated Points
# x_i,y_i=plist.T

# u3=np.linspace(0,1,700,endpoint=True)
# #print(u3)
# out = interpolate.splev(u3,[u,[x,y],3])
# #print(out)

# plt.scatter(x_i,y_i)
# plt.plot(x,y,'k--',label='Control polygon',marker='o',markerfacecolor='red')
# plt.plot(out[0],out[1],'b',linewidth=2.0,label='B-spline curve')
# plt.legend(loc='best')
# plt.gca().set_aspect('equal')
# plt.axis([min(x)-0.5, max(x)+0.5, min(y)-0.5, max(y)+0.5])
# plt.rcParams["figure.figsize"] = [9,9]

# plt.title('Cubic B-spline curve evaluation')
# plt.show()


# In[ ]:





# In[ ]:





