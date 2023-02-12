# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 2023

@author: Tim Sparing
"""

import math as m
import matplotlib as plt
import numpy as np
import pandas as pd

# Inputs

Mdes = 10
gamma = 5/3
thetamin = 0.4
nolines = 10
Steps = 5
Throat = 1

# Prints progress for debug


def Percentage(current, total):
    if round((current/total)*100) != round(((current-1)/total)*100):
        print(f"{round((current/total)*100)}% complete")

# Prandtl-Meyer function


def PrandtlMeyer(Mdes=Mdes, gamma=gamma):
    return ((gamma+1)/(gamma-1))**0.5*m.degrees(m.atan((((gamma-1)/(gamma+1))*(Mdes**2-1))**0.5))-m.degrees(m.atan((Mdes**2-1)**0.5))


PrandtlMeyer = PrandtlMeyer()

# Find the maximum turning angle


def thetamax(Mdes=Mdes):
    return PrandtlMeyer/2


thetamax = thetamax()

# Find the respective turning angles


def theta(nolines=nolines, thetamin=thetamin):
    TurningAngles = []
    for l in range(nolines):
        TurningAngles.append(thetamin+((l*(thetamax-thetamin))/(nolines-1)))
    return TurningAngles


theta = theta()

# Calculate the Riemann Invariants for each chacarteristic line


def CharacteristicRiemann(nolines=nolines):
    CRPlus = []
    CRMinus = []
    for l in range(nolines):
        CRPlus.append(0)
        CRMinus.append(2*(theta[l]))
    return [CRPlus, CRMinus]


CharacteristicRiemann = CharacteristicRiemann()

# Calculate the Riemann Invariants for each point


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

# Calculate Nu and the turning angles from the Riemann Invariants


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
            Theta.append(
                (RiemannInvariants[1][l-1]-RiemannInvariants[0][l-1])/2)
        Nu.append(RiemannInvariants[0][l]+Theta[l])
    return [Theta, Nu]


ThetaAndNu = ThetaAndNu()

# Use the Newton Rhapson method to get the Mach number form the PrandtMeyer function


def PrandtlMeyerNR(Nu, Steps=Steps):
    M = 1.1
    for l in range(Steps):
        M = M - (((gamma+1)/(gamma-1))**0.5*(m.atan((((gamma-1)/(gamma+1))*(M**2-1))**0.5))-(m.atan((M**2-1)**0.5))-m.radians(Nu)) / \
            ((((gamma+1)/(gamma-1))**0.5*((gamma-1)/(gamma+1))**0.5*M) /
             ((((gamma-1)/(gamma+1))*(M**2-1)+1)*(M**2-1)**0.5)-1/(M*(M**2-1)**0.5))
    return M

# Calculate Mu


def Mu(nolines=nolines, Steps=Steps):
    points = nolines
    mu = []
    for l in range(nolines+1):
        points = points + l
    for l in range(points):
        mu.append(
            m.degrees(m.asin(1/(PrandtlMeyerNR(ThetaAndNu[1][l], Steps)))))
    return mu


Mu = Mu()

# Calculates the coordinates in the X-Y plane of the characteristic line intersections and boundary


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
            X[sym] = -1*(Throat/(m.tan(m.radians(0.5*(theta[l] -
                         m.degrees(m.asin(1/(PrandtlMeyerNR(theta[l], Steps))))-Mu[sym])))))
            sym = sym + leap + 1
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
                X[l] = (-1*X[l-1]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l])))+Y[l-1]-1)/(m.tan(m.radians(0.5*(theta[l]-m.degrees(
                    m.asin(1/(PrandtlMeyerNR(theta[l], Steps))))+ThetaAndNu[0][l]-Mu[l])))-m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l]))))
                Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(0.5 *
                                                            (Mu[l-1]+ThetaAndNu[0][l-1]+Mu[l]+ThetaAndNu[0][l])))
            else:
                stack = 0
                count = nolines+1
                comp = 0
                for i in reversed(range(nolines+1)):
                    if comp == 1:
                        comp = 1
                    elif l > count + i:
                        count = count+i
                        stack = stack+1
                    else:
                        sub = nolines - stack
                        # print(X[l-sub], ThetaAndNu[0][l-sub], Mu[l-sub], ThetaAndNu[0]
                        #      [l], Mu[l], ThetaAndNu[0][l-1], Mu[l-1], Y[l-1], Y[l-sub])
                        X[l] = ((X[l-sub]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-sub]-Mu[l-sub]+ThetaAndNu[0][l]-Mu[l]))))-(X[l-1]*m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l])))) +
                                Y[l-1]-Y[l-sub])/(m.tan(m.radians(0.5*(ThetaAndNu[0][l-sub]-Mu[l-sub]+ThetaAndNu[0][l]-Mu[l])))-m.tan(m.radians(0.5*(ThetaAndNu[0][l-1]+Mu[l-1]+ThetaAndNu[0][l]+Mu[l]))))
                        Y[l] = Y[l-1]+(X[l]-X[l-1])*m.tan(m.radians(0.5 *
                                                                    (Mu[l-1]+ThetaAndNu[0][l-1]+Mu[l]+ThetaAndNu[0][l])))
                        count = nolines*nolines
                        comp = 1
        elif X[l] == "Boundary" and Y[l] == "Boundary":
            if l == nolines:
                X[l] = (Y[l-1]-X[l-1]*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))-Throat)/(m.tan(m.radians(
                    0.5*(theta[nolines-1]+ThetaAndNu[0][l])))-m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1])))
                Y[l] = Y[l-1]+(X[l]-X[l-1]) * \
                    m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))
                p = l
            else:
                X[l] = (X[p]*m.tan(m.radians(0.5*(ThetaAndNu[0][p]+ThetaAndNu[0][l])))-X[l-1]*m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1])) +
                        Y[l-1]-Y[p])/(m.tan(m.radians(0.5*(ThetaAndNu[0][p]+ThetaAndNu[0][l])))-m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1])))
                Y[l] = Y[l-1]+(X[l]-X[l-1]) * \
                    m.tan(m.radians(ThetaAndNu[0][l-1]+Mu[l-1]))
                p = l
        Percentage(l, points)
    return [X, Y]


Coordinates = Coordinates()

# Show only boundary coordinates


def Boundary(nolines=nolines):
    pos = -1
    BP = []
    X = [0]
    Y = [Throat]
    for l in reversed(range(nolines+2)):
        if l > 1:
            pos = pos+l
            BP.append(pos)
    for l in BP:
        X.append(Coordinates[0][l])
        Y.append(Coordinates[1][l])
    return [X, Y]


Boundary = Boundary()

# Fit curve to datapoints


def Curve():
    c = np.poly1d(np.polyfit(Boundary[0], Boundary[1], 6))
    plt.pyplot.plot(Boundary[0], c(Boundary[0]))


# Curve()


# Returns all the values (mainly for degbugging)

def All(nolines=nolines, Steps=Steps):
    Mach = []
    cMach = []
    cMu = []
    cx = []
    cy = []
    points = nolines
    for l in range(nolines):
        cx.append(0)
        cy.append(1)
        cMach.append(PrandtlMeyerNR(theta[l], Steps))
        cMu.append(m.degrees(m.asin(1/cMach[l])))
    for l in range(nolines+1):
        points = points + l
    for l in range(points):
        Mach.append(PrandtlMeyerNR(ThetaAndNu[1][l], Steps))
    CPoints = pd.DataFrame({f'R\u207A': RiemannInvariants[0], f'R\u207B': RiemannInvariants[1], f'\u03B8': ThetaAndNu[0],
                           f'\u03BD': ThetaAndNu[1], 'M': Mach, f'\u03BC': Mu, 'x': Coordinates[0], 'y': Coordinates[1]})
    CLines = pd.DataFrame({f'R\u207A': CharacteristicRiemann[0], f'R\u207B': CharacteristicRiemann[1],
                          f'\u03B8': theta, f'\u03BD': theta, 'M': cMach, f'\u03BC': cMu, 'x': cx, 'y': cy})
    CPoints.to_csv("CPoints.csv", sep='\t', encoding='utf-16')
    CLines.to_csv("CLines.csv", sep='\t', encoding='utf-16')


# All()

plt.pyplot.scatter(Boundary[0], Boundary[1])
