#coding: utf-8

import numpy as np
import scipy.sparse as sp
import fonction as k
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.integrate import odeint
import sys
from mpl_toolkits.mplot3d import Axes3D
import math
# reload (sys)
# sys.setdefaultencoding('utf8')


N=int(input('(Nombre intervalles) N = '))# nombre d'intervalle
l=1.
h=l/(N)
j=1
x=np.linspace(0,l,N+1)
dt=float(input('(pas de temps) dt = '))
# s=(x[5]+x[6])/2#1/np.sqrt(2)
s = 0.8
Ns=int(s/h)+1
tmin=float(input('(temps initial) tmin = '))
tmax=float(input('(temps final) tmax = '))
Tg=10.
Td=1.
Ts=0.
a=1.
F=np.zeros(N+1)
F[0]=-Tg
F[1:Ns-1]=-Tg+Tg*(x[1:Ns-1]/s)
F[Ns:N]=(x[Ns:N] -s)*Td/(1-s)
F[Ns-1]=0
F[Ns]=0
F[N]=Td
while tmin<=tmax :
    Ns=int(s/h)+1
    delta=s-h*(Ns-1)
    alpha=2*h*h +3*h*delta +delta*delta
    lamda=(delta*h+delta*delta)/alpha
    gama=(4*h*delta+2*delta*delta)/alpha
    #print("Ns= ",Ns," delta=",delta," lamda= ",lamda," gama=",gama)


    #B=k.matrice(h,gama,lamda,N,Ns)
    #A=k.matriceo1(h,delta,N,Ns)
    #F=k.secf(h,Tg,x,Ts,N,Ns,a)
    #sol=k.exact(x,s,N,Ns,Ts,a)
    alphacd=-5*delta*h +delta*delta +6*h*h
    mu=(h-delta)/(3*h-delta)
    beta=(-2*delta*delta+8*delta*h- 6*h*h)/alphacd
    #print("beta= ",alphacd," sigma=  ",mu," ro=  ",beta)
    #C=k. matriceCD(h,beta,mu,N,Ns)
    #sol=k.exactCD(x,s,N,Ns,Ts,a)
    #F=k.secfCD(h,Td,x,Ts,N,Ns,a)

    A=k.matriceprincipale(h,gama,lamda,beta,mu,N,Ns)
    #F=k.secfprincipale(h,Td,Tg,x,Ts,N,Ns,a)
    #sol=k.solexacte(x,t,Ns,s)#exactprincipale(x,s,N,Ns,Ts,a)
    #T=np.linalg.solve(A.todense(),F)
    #F=k.solexacte(x,0.,Ns,s)

    #I=np.eye(N+1)
    I=sp.csr_matrix(np.identity(N+1))
    I[Ns-1,Ns-1]=0
    I[Ns,Ns]=0

    #T=np.zeros((N+1,N+1))
    #print('second_menmbre_initiale ',' ')
    #print(F/dt)
    #print("\n")
    F[Ns-1]=0
    F[Ns]=0
    #T=sp.linalg.spsolve((I/(dt))-A/(h*h),F/dt)
    T=np.linalg.solve((I.todense()/(dt))-A.todense()/(h*h),F/dt)
    #print(T[Ns+1]," ",T[Ns+2])
    #print(T)
    dTcd = ((h+delta)/(2*h*h-h*delta) - 2*(delta-h)/(2*h*h-h*delta) )*T[Ns+1] + (-delta/(3*h*h-h*delta) + 2*(delta-h)/(3*h*h-h*delta))*T[Ns+2]
    dTcg = ((delta-2*h)/(h*h+h*delta) - 2*delta/(h*h +h*delta))*T[Ns-2] + ((h-delta)/(2*h*h + h*delta) + 2*delta*1/(2*h*h+h*delta))*T[Ns-3]
    #print("dTcd= ",dTcd)
    #print("dTcg= ",dTcg)
    #s1=s +dt*(dTcg-dTcd)
    #if(s1-s>=0.5*h):
    #    dt=h/2
    si=s +dt*(dTcg-dTcd)
    print("s= ",s)
    if si-s>=0.7*h:
        while si-s>=0.7*h:
            dt=dt/2
            si=s +dt*(dTcg-dTcd)
            print(dt," ",si)
    else:
        s= s +dt*(dTcg-dTcd)
    print("s= ",s)
    if j == 1:
        line, = plt.plot(x, T)
    else:
        line.set_ydata(T)

    plt.pause(0.0001) # pause avec duree en secondes
    tmin=tmin+dt
    j=j+1
    F=T
    tmin=tmin+dt

#print(F)
#print("\n")
#print(T)
#np.savetxt('matrice.txt',(I/(dt))-A.todense()/(h*h),fmt='%f',delimiter=' ')
#np.savetxt('second_membre_fin.txt', F,fmt='%f',delimiter=' ')
#np.savetxt('solution_exacte.txt', sol,fmt='%f',delimiter=' ')
#np.savetxt('solution_approche.txt', T,fmt='%f',delimiter=' ')



print(tmin," ",tmax)


#erreur=np.abs(F-A.dot(sol))
#erreurg=np.abs(T[1:Ns-1]-sol[1:Ns-1])
#erreurd=np.abs(T[Ns+1:N]-sol[Ns+1:N])
#print(T[Ns-2]," ",T[Ns])
#print (' max erreur gauche ',np.max(erreurg) )
#print (' max erreur droite ',np.max(erreurd) )


#plt.plot(x,sol,'b',label=" solution_exacte ")
#plt.plot(x,T,'r',label="solutuion_approché ")


#erreur en O(dx²)


#erreur en O(dt)
#plt.grid(True)
#ertg=np.array([ 1.286e-11, 2.573e-11, 5.146e-11,1.028e-10,1.285e-10])
#ertd=np.array([1.189e-09, 2.382e-09,4.753e-09,9.509e-09, 1.188e-08])
#NPt=np.array([ 1e-10,2e-10,4e-10, 8e-10, 10e-10])
#plt.title("erreur en O(dt) Nx=300 et nombre itération=10*dt")
#plt.plot(np.log(NPt),np.log(ertg),'b',linewidth=0.8,label=" erreur_temps_coté_gauche ")
#plt.plot(np.log(NPt),np.log(ertd),'r',linewidth=0.8,label=" erreur_temps_coté_droit ")
#plt.ylabel('log(|T_ex_fina-T_app_final|)')
#plt.xlabel('dt')


plt.legend()
plt.show()
