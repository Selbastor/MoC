# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 2023

@author: Tim
"""

import math as m
import matplotlib.pyplot as plt
import pandas as pd

# Inputs

Mdes = 2.8
gamma = 5 / 3
thetamin = 0.4
nolines = 100
Steps = 4
Throat = 1
Perc = 100


# Prints progress for debug


def Percentage(current, total):
    if round((current / total) * 100) != round(((current - 1) / total) * 100):
        print(f"{round((current/total)*100)}% complete")


# Prandtl-Meyer function


def PrandtlMeyer(Mdes=Mdes, gamma=gamma):
    return ((gamma + 1) / (gamma - 1)) ** 0.5 * m.degrees(
        m.atan((((gamma - 1) / (gamma + 1)) * (Mdes**2 - 1)) ** 0.5)
    ) - m.degrees(m.atan((Mdes**2 - 1) ** 0.5))


# Find the maximum turning angle


def thetamax(PM, Mdes=Mdes):
    return PM / 2


# Find the respective turning angles


def Ctheta(TMax, nolines=nolines, thetamin=thetamin):
    TurningAngles = []
    for l in range(nolines):
        TurningAngles.append(
            thetamin + ((l * (TMax - thetamin)) / (nolines - 1))
        )
    return TurningAngles


# Calculate the Riemann Invariants for each chacarteristic line


def CharacteristicRiemann(theta, nolines=nolines):
    CRPlus = []
    CRMinus = []
    for l in range(nolines):
        CRPlus.append(0)
        CRMinus.append(2 * (theta[l]))
    return [CRPlus, CRMinus]


# Calculate the Riemann Invariants for each point


def RiemannInvariants(CR, nolines=nolines):
    points = nolines
    for l in range(nolines + 1):
        points = points + l
    Riemanns = []
    for l in range(points):
        Riemanns.append(l)
    pos = 0
    RPlus = []
    RMinus = []
    for l in range(nolines):
        for i in range(nolines - l + 1):
            RPlus.append(CR[1][l])
        for i in range(nolines):
            if i >= pos:
                RMinus.append(CR[1][i])
                if i == nolines - 1:
                    RMinus.append("N/A")
        pos = pos + 1
    return [RPlus, RMinus]


# Calculate Nu and the turning angles from the Riemann Invariants


def GetThetaAndNu(RI, nolines=nolines):
    points = nolines
    for l in range(nolines + 1):
        points = points + l
    Theta = []
    Nu = []
    for l in range(points):
        if RI[1][l] != "N/A":
            Theta.append((RI[1][l] - RI[0][l]) / 2)
        else:
            Theta.append((RI[1][l - 1] - RI[0][l - 1]) / 2)
        Nu.append(RI[0][l] + Theta[l])
    return [Theta, Nu]


# Use the Newton Rhapson method to get the Mach number form the PrandtlMeyer function


def PrandtlMeyerNR(Nu, Steps=Steps):
    M = 1.1
    for l in range(Steps):
        M = M - (
            ((gamma + 1) / (gamma - 1)) ** 0.5
            * (m.atan((((gamma - 1) / (gamma + 1)) * (M**2 - 1)) ** 0.5))
            - (m.atan((M**2 - 1) ** 0.5))
            - m.radians(Nu)
        ) / (
            (
                ((gamma + 1) / (gamma - 1)) ** 0.5
                * ((gamma - 1) / (gamma + 1)) ** 0.5
                * M
            )
            / (
                (((gamma - 1) / (gamma + 1)) * (M**2 - 1) + 1)
                * (M**2 - 1) ** 0.5
            )
            - 1 / (M * (M**2 - 1) ** 0.5)
        )
    return M


# Calculate Mu


def GetMu(ThetaAndNu, nolines=nolines, Steps=Steps):
    points = nolines
    mu = []
    for l in range(nolines + 1):
        points = points + l
    for l in range(points):
        mu.append(
            m.degrees(m.asin(1 / (PrandtlMeyerNR(ThetaAndNu[1][l], Steps))))
        )
    return mu


# Calculates the coordinates in the X-Y plane of the characteristic line intersections and boundary


def Coordinates(
    theta, Mu, ThetaAndNu, nolines=nolines, Throat=Throat, Steps=Steps
):
    points = nolines
    for l in range(nolines + 1):
        points = points + l
    X = [None] * points
    Y = [None] * points
    sym = 0
    leap = nolines
    reflec = []
    for l in range(nolines):
        if leap > 0:
            Y[sym] = 0
            if sym == 0:
                X[sym] = -Throat / (
                    m.tan(
                        m.radians(
                            0.5
                            * (
                                theta[l]
                                - m.degrees(
                                    m.asin(
                                        1 / (PrandtlMeyerNR(theta[l], Steps))
                                    )
                                )
                                - Mu[sym]
                            )
                        )
                    )
                )
            if sym != 0:
                reflec.append(sym)
            sym = sym + leap + 1
            leap = leap - 1
            if sym > 0:
                X[sym - 1] = "Boundary"
                Y[sym - 1] = "Boundary"
        X[-1] = "Boundary"
        Y[-1] = "Boundary"
    p = 0
    for l in range(points):
        if X[l] == None and Y[l] == None or l in reflec:
            if l < nolines:
                X[l] = (
                    -1
                    * X[l - 1]
                    * m.tan(
                        m.radians(
                            0.5
                            * (
                                ThetaAndNu[0][l - 1]
                                + Mu[l - 1]
                                + ThetaAndNu[0][l]
                                + Mu[l]
                            )
                        )
                    )
                    + Y[l - 1]
                    - 1
                ) / (
                    m.tan(
                        m.radians(
                            0.5
                            * (
                                theta[l]
                                - m.degrees(
                                    m.asin(
                                        1 / (PrandtlMeyerNR(theta[l], Steps))
                                    )
                                )
                                + ThetaAndNu[0][l]
                                - Mu[l]
                            )
                        )
                    )
                    - m.tan(
                        m.radians(
                            0.5
                            * (
                                ThetaAndNu[0][l - 1]
                                + Mu[l - 1]
                                + ThetaAndNu[0][l]
                                + Mu[l]
                            )
                        )
                    )
                )
                Y[l] = Y[l - 1] + (X[l] - X[l - 1]) * m.tan(
                    m.radians(
                        0.5
                        * (
                            Mu[l - 1]
                            + ThetaAndNu[0][l - 1]
                            + Mu[l]
                            + ThetaAndNu[0][l]
                        )
                    )
                )
            else:
                stack = 0
                count = nolines + 1
                comp = 0
                for i in reversed(range(nolines + 1)):
                    if comp == 1:
                        comp = 1
                    elif l >= count + i:
                        count = count + i
                        stack = stack + 1
                    elif l in reflec:
                        sub = nolines - stack
                        X[l] = X[l - sub] - (
                            Y[l - sub]
                            / (
                                m.tan(
                                    m.radians(
                                        0.5
                                        * (
                                            ThetaAndNu[0][l - sub]
                                            - Mu[l - sub]
                                            + ThetaAndNu[0][l]
                                            - Mu[l]
                                        )
                                    )
                                )
                            )
                        )
                        comp = 1
                    else:
                        sub = nolines - stack
                        X[l] = (
                            (
                                X[l - sub]
                                * m.tan(
                                    m.radians(
                                        0.5
                                        * (
                                            ThetaAndNu[0][l - sub]
                                            - Mu[l - sub]
                                            + ThetaAndNu[0][l]
                                            - Mu[l]
                                        )
                                    )
                                )
                            )
                            - (
                                X[l - 1]
                                * m.tan(
                                    m.radians(
                                        0.5
                                        * (
                                            ThetaAndNu[0][l - 1]
                                            + Mu[l - 1]
                                            + ThetaAndNu[0][l]
                                            + Mu[l]
                                        )
                                    )
                                )
                            )
                            + Y[l - 1]
                            - Y[l - sub]
                        ) / (
                            m.tan(
                                m.radians(
                                    0.5
                                    * (
                                        ThetaAndNu[0][l - sub]
                                        - Mu[l - sub]
                                        + ThetaAndNu[0][l]
                                        - Mu[l]
                                    )
                                )
                            )
                            - m.tan(
                                m.radians(
                                    0.5
                                    * (
                                        ThetaAndNu[0][l - 1]
                                        + Mu[l - 1]
                                        + ThetaAndNu[0][l]
                                        + Mu[l]
                                    )
                                )
                            )
                        )
                        Y[l] = Y[l - 1] + (X[l] - X[l - 1]) * m.tan(
                            m.radians(
                                0.5
                                * (
                                    Mu[l - 1]
                                    + ThetaAndNu[0][l - 1]
                                    + Mu[l]
                                    + ThetaAndNu[0][l]
                                )
                            )
                        )
                        comp = 1
        elif X[l] == "Boundary" and Y[l] == "Boundary":
            if l == nolines:
                X[l] = (
                    Y[l - 1]
                    - X[l - 1]
                    * m.tan(m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1]))
                    - Throat
                ) / (
                    m.tan(
                        m.radians(0.5 * (theta[nolines - 1] + ThetaAndNu[0][l]))
                    )
                    - m.tan(m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1]))
                )
                Y[l] = Y[l - 1] + (X[l] - X[l - 1]) * m.tan(
                    m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1])
                )
                p = l
            else:
                X[l] = (
                    X[p]
                    * m.tan(
                        m.radians(0.5 * (ThetaAndNu[0][p] + ThetaAndNu[0][l]))
                    )
                    - X[l - 1]
                    * m.tan(m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1]))
                    + Y[l - 1]
                    - Y[p]
                ) / (
                    m.tan(
                        m.radians(0.5 * (ThetaAndNu[0][p] + ThetaAndNu[0][l]))
                    )
                    - m.tan(m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1]))
                )
                Y[l] = Y[l - 1] + (X[l] - X[l - 1]) * m.tan(
                    m.radians(ThetaAndNu[0][l - 1] + Mu[l - 1])
                )
                p = l
        Percentage(l, points)
    return [X, Y]


# Show only boundary coordinates


def BoundaryPoints(Coords, nolines=nolines, Throat=Throat):
    pos = -1
    BP = []
    X = [0]
    Y = [Throat]
    for l in reversed(range(nolines + 2)):
        if l > 1:
            pos = pos + l
            BP.append(pos)
    for l in BP:
        X.append(Coords[0][l])
        Y.append(Coords[1][l])
    return [X, Y]


# Returns all the values (mainly for degbugging)


def All(theta, RI, ThetaAndNu, Mu, Coords, CR, nolines=nolines, Steps=Steps):
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
        cMu.append(m.degrees(m.asin(1 / cMach[l])))
    for l in range(nolines + 1):
        points = points + l
    for l in range(points):
        Mach.append(PrandtlMeyerNR(ThetaAndNu[1][l], Steps))
    CPoints = pd.DataFrame(
        {
            f"R\u207A": RI[0],
            f"R\u207B": RI[1],
            f"\u03B8": ThetaAndNu[0],
            f"\u03BD": ThetaAndNu[1],
            "M": Mach,
            f"\u03BC": Mu,
            "x": Coords[0],
            "y": Coords[1],
        }
    )
    CLines = pd.DataFrame(
        {
            f"R\u207A": CR[0],
            f"R\u207B": CR[1],
            f"\u03B8": theta,
            f"\u03BD": theta,
            "M": cMach,
            f"\u03BC": cMu,
            "x": cx,
            "y": cy,
        }
    )
    CPoints.to_csv("CPoints.csv", sep="\t", encoding="utf-16")
    CLines.to_csv("CLines.csv", sep="\t", encoding="utf-16")


# Shortens the Nozzle


def Cutoff(BP, Perc=Perc):
    Max = BP[0][-1]
    CutX = []
    CutY = []
    point = 0
    for l in BP[0]:
        if l <= (Perc / 100) * Max:
            CutX.append(l)
            CutY.append(BP[1][point])
            point = point + 1
    return [CutX, CutY]


# Show characteristic lines


def ShowCLines(BP, color, size, Throat=Throat):
    LX = []
    LY = []
    point = 0
    for l in range(nolines + 1):
        if l != 0:
            LX.append(BP[0][l])
            LY.append(-1 * BP[1][l])
    for l in LX:
        lX = [0]
        lY = [Throat]
        lX.append(l)
        lY.append(LY[point])
        rY = [i * -1 for i in lY]
        plt.plot(lX, lY, c=color, lw=size)
        plt.plot(lX, rY, c=color, lw=size)
        point = point + 1


# 1D approximation


def oneD(Mdes, gamma):
    Ae = (1 / Mdes) * (
        ((2 / (gamma + 1)) * (1 + (((gamma - 1) / 2) * (Mdes**2))))
        ** ((gamma + 1) / (2 * (gamma - 1)))
    )
    return Ae


def Rao(Mdes, gamma):
    thetaE = []
    thetaN = []
    ye = []
    # ye = P*xe+Q+(S*xe+T)**0.5
    xe = []
    P = (
        ye * m.tan(m.radians(thetaN))
        + ye * m.tan(m.radians(thetaE))
        - 2 * xe * m.tan(m.radians(thetaE)) * m.tan(m.radians(thetaN))
    ) / (2 * ye - xe * m.tan(m.radians(thetaN)) - xe * m.tan(m.radians(thetaE)))
    S = (((ye - P * xe) ** 2) * (m.tan(m.radians(thetaN)) - P)) / (
        xe * m.tan(m.radians(thetaN)) - ye
    )
    Q = (S) / (2 * (m.tan(m.radians(thetaN)) - P))
    T = Q**2


def Conical():
    Epsilon = []
    Rt = []
    R1 = []
    thetaE = []
    X = []
    Y = []
    L = (
        Rt * ((Epsilon**0.5) - 1) + R1 * (1 / (m.cos(m.radians(thetaE)) - 1))
    ) / m.tan(m.radians(thetaE))
    return L
