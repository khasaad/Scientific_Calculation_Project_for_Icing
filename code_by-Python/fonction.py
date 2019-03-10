import scipy.sparse as sp
import numpy as np
import math


def matrice(h,gama,lamda,N,Ns):
    dim=N+1
    diags=np.zeros((4,dim))

    #diagonale principale
    diags[2,0]=1.
    diags[2,1:Ns]=-2./(h*h)
    diags[2,Ns-1:]=1.

    #diagonale +1
    diags[3,2:Ns]=1./(h*h)
    #diagonal -1
    diags[1,:Ns-1]=1./(h*h)
    diags[1,Ns-2]=-gama

    #diagonale -2
    diags[0,Ns-3]=lamda


    # Construction de la matrice creuse A
    A = sp.spdiags(diags,[-2,-1,0,1], dim, dim, format="csr")
    return h*h* A


def f(x,a):
    return a

def solution(x,s,a):
    return a*x*(x-s)/2

def secf(h,Tg,x,Ts,N,Ns,a):
    b=np.zeros((N+1))
    b[0]=Tg
    for i in np.arange(1,Ns-1):
        b[i]=f(x[i],a)
    b[Ns-1]=0
    b[Ns:]=Ts

    return h*h*b

def exact(x,s,N,Ns,Ts,a):
    b=np.zeros((N+1))
    for i in np.arange(Ns):
        b[i]=solution(x[i],s,a)
    b[Ns:]=Ts
    return  b



def matriceo1(h,delta,N,Ns):
    dim=N+1
    diags=np.zeros((3,dim))

    #diagonale principale
    diags[1,0]=1.
    diags[1,1:Ns-1]=-2./(h*h)
    diags[1,Ns-1]=-1.
    diags[1,Ns:]=1.

    #diagonale +1
    diags[2,2:Ns]=1./(h*h)

    #diagonal -1
    diags[0,:Ns-2]=1./(h*h)
    diags[0,Ns-2]=delta/(delta+h)


    # Construction de la matrice creuse A
    A = sp.spdiags(diags,[-1,0,1], dim, dim, format="csr")
    return h*h* A

def matriceCD(h,beta,mu,N,Ns):
    dim=N+1
    diags=np.zeros((4,dim))

    #diagonale principale
    diags[1,:Ns+1]=1.
    diags[1,Ns+1:N]=-2./(h*h)
    diags[1,N]=1.

    #diagonale +1
    diags[2,Ns+1:]=1./(h*h)
    diags[2,Ns+1]=beta

    #diagonal -1
    diags[0,Ns:N-1]=1./(h*h)
    diags[0,N-1]=0.


    #diagonale +2
    diags[3,Ns+2]=mu


    # Construction de la matrice creuse A
    A = sp.spdiags(diags,[-1,0,1,2], dim, dim, format="csr")
    return h*h*A

def secfCD(h,Td,x,Ts,N,Ns,a):
    b=np.zeros((N+1))
    b[:Ns-1]=Ts
    b[Ns]=0
    for i in np.arange(Ns,N):
        b[i]=f(x[i],a)
    b[N]=Td
    return h*h*b

def solutionCd(x,s,a):
    return a*(x*x -(s+1)*x +s)/2

def exactCD(x,s,N,Ns,Ts,a):
    b=np.zeros((N+1))
    b[:Ns]=Ts
    for i in np.arange(Ns,N+1):
        b[i]=solutionCd(x[i],s,a)
    return  b

def matriceprincipale(h,gama,lamda,beta,mu,N,Ns):
    dim=N+1
    diags=np.zeros((5,dim))

    #diagonale principale
    diags[2,0]=1.
    diags[2,1:Ns]=-2./(h*h)
    diags[2,Ns-1:]=1.
    diags[2,Ns:Ns+1]=1.
    diags[2,Ns+1:N]=-2./(h*h)
    diags[2,N]=1.

    #diagonale +1
    diags[3,2:Ns]=1./(h*h)
    diags[3,Ns+1:]=1./(h*h)
    diags[3,Ns+1]=beta
    #diagonal -1
    diags[1,:Ns-1]=1./(h*h)
    diags[1,Ns-2]=-gama
    diags[1,Ns:N-1]=1./(h*h)
    diags[1,N-1]=0.

    #diagonale -2
    diags[0,Ns-3]=lamda

    #diagonale +2
    diags[4,Ns+2]=mu

    # Construction de la matrice creuse A
    A = sp.spdiags(diags,[-2,-1,0,1,2], dim, dim, format="csr")
    return h*h*A

def secfprincipale(h,Td,Tg,x,Ts,N,Ns,a):
    b=np.zeros((N+1))
    b[0]=Tg
    for i in np.arange(1,Ns-1):
        b[i]=f(x[i],a)
    for i in np.arange(Ns-1,N):
        b[i]=f(x[i],a)
    b[Ns-1]=0.
    b[Ns]=0
    b[N]=Td
    return h*h*b

def exactprincipale(x,s,N,Ns,Ts,a):
    b=np.zeros((N+1))
    for i in np.arange(Ns):
        b[i]=solution(x[i],s,a)
    for i in np.arange(Ns,N):
        b[i]=solutionCd(x[i],s,a)
    b[N]=solutionCd(x[N],s,a)
    return  b

def qd(x,s,t):
    return np.sin( np.pi*(x-s)/(1-s) )*np.exp( -np.pi*np.pi*t/((1-s)*(1-s)) )

def qg(x,s,t):
  
    return np.sin( np.pi*x/s)*np.exp( -(np.pi*np.pi/(s*s))*t)

def solexacte(x,t,Ns,s):
    b=np.zeros((np.size(x)))
    for i in np.arange(np.size(x)):
        if(i<Ns):
            b[i]=qg(x[i],s,t)
        else:
            b[i]=qd(x[i],s,t)
    return b        
