# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 2023

@author: Tim Sparing
"""

import math as m
import matplotlib as plt

#Inputs

Mdes = 2.8
gamma = 1.4
thetamin = 0.4
nolines = 5
Steps = 10
Throat = 1

#Prandtl-Meyer function
def PrandtlMeyer(Mdes=Mdes, gamma=gamma):
    return ((gamma+1)/(gamma-1))**0.5*m.degrees(m.atan((((gamma-1)/(gamma+1))*(Mdes**2-1))**0.5))-m.degrees(m.atan((Mdes**2-1)**0.5))

PrandtlMeyer = PrandtlMeyer()

#Find the maximum turning angle
def thetamax(Mdes=Mdes):
    return PrandtlMeyer/2

thetamax = thetamax()

#Find the respectice turning angles
def theta(nolines=nolines, thetamin=thetamin):
    TurningAngles = []
    for l in range(nolines):
        TurningAngles.append(thetamin+((l*(thetamax-thetamin))/(nolines-1)))
    return TurningAngles
        
theta = theta()

#Calculate the Riemann Invariants for each chacarteristic line
def CharacteristicRiemann(nolines=nolines):
    CRPlus = []
    CRMinus = []
    for l in range(nolines):
        CRPlus.append(0)
        CRMinus.append(2*(theta[l]))
    return [CRPlus, CRMinus]

CharacteristicRiemann = CharacteristicRiemann()

#Calculate the Riemann Invariants for each point
def RiemannInvariants(nolines=nolines):
    points = nolines
    for l in range(nolines+1):
        points = points + l
    Riemanns = []
    for l in range(points):
        Riemanns.append(l)
    pos = 0
    RPlus = []
    RMinus = []
    for l in range(nolines):
        for i in range(nolines-l+1):
            RPlus.append(CharacteristicRiemann[1][l])
        for i in range(nolines):
            if i >= pos:
                RMinus.append(CharacteristicRiemann[1][i])
                if i == nolines-1:
                    RMinus.append("N/A")
        pos = pos + 1
    return [RPlus, RMinus]

RiemannInvariants = RiemannInvariants()

#Calculate Nu and the turning angles from the Riemann Invariants
def ThetaAndNu(nolines=nolines):
    points = nolines
    for l in range(nolines+1):
        points = points + l
    Theta = []
    Nu = []
    for l in range(points):
        if RiemannInvariants[1][l] != "N/A":
            Theta.append((RiemannInvariants[1][l]-RiemannInvariants[0][l])/2)
        else:
            Theta.append((RiemannInvariants[1][l-1]-RiemannInvariants[0][l-1])/2)
        Nu.append(RiemannInvariants[0][l]+Theta[l])
    return [Theta, Nu]

ThetaAndNu = ThetaAndNu()
            
#Use the Newton Rhapson method to get the Mach number form the PrandtMeyer function
def PrandtlMeyerNR(Nu, Steps=Steps):
    M = 1.1
    for l in range(Steps):
        M = M - (((gamma+1)/(gamma-1))**0.5*(m.atan((((gamma-1)/(gamma+1))*(M**2-1))**0.5))-(m.atan((M**2-1)**0.5))-m.radians(Nu))/((((gamma+1)/(gamma-1))**0.5*((gamma-1)/(gamma+1))**0.5*M)/((((gamma-1)/(gamma+1))*(M**2-1)+1)*(M**2-1)**0.5)-1/(M*(M**2-1)**0.5))
    return M

#Calculate Mu
def Mu(nolines=nolines, Steps=Steps):
    points = nolines
    mu = []
    for l in range(nolines+1):
        points = points + l
    for l in range(points):
        mu.append(m.degrees(m.asin(1/(PrandtlMeyerNR(ThetaAndNu[1][l], Steps)))))
    return mu

Mu = Mu()
        
def Coordinates(nolines=nolines, Throat=Throat, Steps=Steps):
    points = nolines
    for l in range(nolines+1):
        points = points + l
    X = [None]*points
    Y = [None]*points
    sym = 0
    leap = nolines
    for l in range(nolines):
        if leap > 0:
            Y[sym] = 0
            X[sym] = -1*(Throat/(m.tan(m.radians(0.5*(theta[l]-m.degrees(m.asin(1/(PrandtlMeyerNR(theta[l], Steps))))-Mu[sym])))))
            sym = sym + leap +1
            leap = leap - 1
            if sym > 0:
                X[sym-1] = "Boundary"
                Y[sym-1] = "Boundary"
        X[-1] = "Boundary"
        Y[-1] = "Boundary"
    p = 0
    for l in range(points):
        if X[l] == None and Y[l] == None:
            if l < nolines:
                X[l] = (-1*X[l-1]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l])))+Y[l-1]-1)/(m.tan(m.radians(0.5*(theta[l]-m.degrees(m.asin(1/(PrandtlMeyerNR(theta[l], Steps))))+ThetaAndNu[0][l]-Mu[l])))-m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l]))))
                Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(0.5*(Mu[l-1]+ThetaAndNu[0][l-1]+Mu[l]+ThetaAndNu[0][l])))
            else:
                sub = 0
                count = nolines+1
                for i in reversed(range(nolines+1)):
                    if l > count + i:
                        count = count+i
                        sub = sub+1
                    else:
                        sub = nolines-sub
                        X[l] = ((X[l-sub]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-sub]-Mu[l-sub]+ThetaAndNu[0][l]-Mu[l]))))-(X[l-1]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l]))))+Y[l-1]-Y[l-sub])/(m.tan(m.radians(0.5*(ThetaAndNu[0][l-sub]-Mu[l-sub]+ThetaAndNu[0][l]-Mu[l])))-m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l]))))
                        Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(0.5*(Mu[l-1]+ThetaAndNu[0][l-1]+Mu[l]+ThetaAndNu[0][l])))
                        count = -1*nolines
        elif X[l] == "Boundary" and Y[l] == "Boundary":
            if l == nolines:
                X[l] = (Y[l-1]-X[l-1]*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))-Throat)/(m.tan(m.radians(0.5*(theta[nolines-1]+ThetaAndNu[0][l])))-m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1])))
                Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))
                p = l
            else:
                print(p)
                print(X[p], Y[p])
                print(X[l-1], Y[l-1])
                print(Mu[p], ThetaAndNu[0][p])
                X[l] = (X[p]*m.tan(m.radians(0.5*(ThetaAndNu[0][p]+ThetaAndNu[0][l])))-X[l-1]*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))+Y[l-1]-Y[p])/(m.tan(m.radians(0.5*(ThetaAndNu[0][p]+ThetaAndNu[0][l])))-m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1])))
                Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))
                p = l
    return [X, Y]        

Coordinates = Coordinates()

plt.pyplot.scatter(Coordinates[0], Coordinates[1])
    

