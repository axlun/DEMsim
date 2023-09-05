import numpy as np
from math import pi
from scipy.optimize import fsolve
from scipy.optimize import fmin


#Small python class for bookkeeping the spheres conviniently
class Sphere():
    def __init__(self, R, mat):
        self.R = R
        self.E = mat[0]
        self.E0 = mat[0]/(1-mat[1]**2)
        self.params = len(mat)
        self.mat = mat
        self.idealPlastic=False
        self.Hmax = 3
        if self.params ==3:  #Ideal plastic material
            self.idealPlastic=True
            self.sY = mat[2]
        elif self.params ==4: # strain hardening material according to Eq. (6)
            self.sY = (mat[2]**mat[3]/mat[0])**(1./(mat[3]-1)) #Eq. (7)
            if mat[3] < 5:
                self.Hmax = 2.8
        else:  #Three parameter strain hardening material Eq. (8)
            self.sY = mat[2]
            if mat[3] >= 10 and mat[4] < 5:
                self.Hmax = 2.8

        #Determining the parameter log(a/R0)_FP, in the code denoted xHmax
        def  eq_xHmax(x):
            return mat[0]/(1-mat[1]**2)/sRep(x,self,1)*x -160
        self.xHmax =  np.log(fsolve(eq_xHmax,1)) 

def Lambda(a, s,R):
    return s.E0/sRep(a,s,R)*a/R       #Eq. (4)     

def sRep(a, s, R):
    #Representative stress Eq. (11)
    if s.idealPlastic or 0.2*a/R*s.E < s.sY:
        return s.sY
    else:
        if s.params == 4:
            return s.mat[2]*(0.2*a/R)**(1./s.mat[3])
        elif s.params == 5:
            return s.sY*(1+s.mat[3]*(0.2*a/R)**(1./s.mat[4])) #Eq. (12)
    return 0 #Should never happen
            
def HbarLargeDef(a, R,s,aLD):
    L = Lambda(a,s,R)
    #Calculating the normalized hardness according to Eq. (13)
    if L < 3:
        H = 4.0/(3.0*pi)*L
    elif 4.0/(3.0*pi)*L < s.Hmax - 0.0965*(np.log(L)-np.log(160))**2:
        H = 4.0/(3.0*pi)*L
    elif L < 160:
        H = s.Hmax - 0.0965*(np.log(L)-np.log(160))**2
    else:
        H = s.Hmax
        
    #Compensating for large deformations
    if a/R > aLD:
        # print(a/R)
        # print(aLD)

        q = np.log(a/R) - np.log(aLD)        #Eq. (27)
        f = 0.0460789*q**2.2139              #Eq. (28)
        phi = q/(s.xHmax - np.log(aLD))      #Eq. (30)
        if phi<1 and phi>0:
            return phi*H*(1-f) + (1-phi)*H   #Eq. (29)
        return H*(1-f)                       #Eq. (28)
    return H                                 #If a < a_LD

def Hbar(a, s,R):
    L = Lambda(a,s,R)
    #Calculating the normalized hardness according to Eq. (13)
    if L < 3:
        H = 4.0/(3.0*pi)*L
    elif 4.0/(3.0*pi)*L < s.Hmax - 0.0965*(np.log(L)-np.log(160))**2:
        H = 4.0/(3.0*pi)*L
    elif L < 160:
        H = s.Hmax - 0.0965*(np.log(L)-np.log(160))**2
    else:
        H = s.Hmax
    return H

def generateContactData(s1,s2,hMax,dF, largeDef=True):
    def getContactRadiusAndCurvature(F, s1, s2, x0,largeDef):
        #Function that deterimines the contact radius and the contact curvature
        def calc_a(var, *data):
            a, R0 = var
            s1, s2 = data

            #Defiing effective radii according to Eq. (9)
            #Observe that in the code that R0 is the curvature and not the radii
            #of the assumed spherical contact surfacce
            R1 = 1./(1./s1.R - R0)   
            R2 = 1./(1./s2.R + R0)
            if largeDef ==True:
                #Calculating (a/R0)_LD according to Eq. (25)
                def onsetOfHDrop(x):
                    #Parameters according to Eq. (26)
                    A = 0.041
                    B = 2*(0.1-A)/pi
                    C = np.tan((0.072-A)/B)

                    # print('x: ' + str(x))
                    if x < -1:
                        x = -1
                    if x>0:
                        return A + B*np.arctan(C*x)
                    return B*C*(x+x**2) + A*(1- x**2)
                aLD2 = onsetOfHDrop(+s2.R*R0)
                aLD1 = onsetOfHDrop(-s1.R*R0)
                H1 = HbarLargeDef(a,R1,s1,aLD1)*sRep(a,s1,R1) 
                H2 = HbarLargeDef(a,R2,s2,aLD2)*sRep(a,s2,R2)
            else:
                H1 = Hbar(a,s1,R1)*sRep(a,s1,R1)
                H2 = Hbar(a,s2,R2)*sRep(a,s2,R2)

            #returning a function that should be minimized in order to solve
            #Eq. (10)
            # print(F)
            # print((1 - H2/H1)**2 + (H1*pi*a**2/F - 1)**2)
            return (1 - H2/H1)**2 + (H1*pi*a**2/F - 1)**2

        data = s1,s2
        #Solving Eq. 10
        x = fmin(calc_a, x0, args=data, ftol=1E-12, disp=False,
                 maxfun=1E6,maxiter=1E6)
        a,R = x
        # print('a: ' + str(a))
        R1 = 1./(1./s1.R - R)
        R2 = 1./(1./s2.R + R)
        # print('R0: ' +str(R) + '\nR1: ' + str(R1) + '\nR2: ' + str(R2))
        return a, R    #Contact radius and curvature are calculated

    def computeIndentationDepth(F, s1,s2, x0, largeDef=False):

        def c2(L,s): #Calculates c2 according to Eq. (16) - Eq. (23)
            #Function for determining c2 in the fully plastic region
            def c2Max(s):  
                if s.idealPlastic:
                    return 1.43
                if s.params==4:
                    #2 parameter plasticity model, use Eq. (22) directly for
                    #c2Max
                    return 1.43*np.exp(-0.97/s.mat[3])

                k = s1.mat[3]
                m = s1.mat[4]
                c2M = 1.43*np.exp(-0.97/m)
                A = 1.43 - c2M
                b = 0.27
                return A*np.exp(-b/A*np.log(1+k)) + c2M #Eq. (23)

            #Calculating c2 according to Eq. (16)
            c2M = c2Max(s)
            L0 = 4
            Lmax = 1527
            if L < L0:
                return 0.5
            elif L < Lmax:
                LP = np.exp((0.5 - 0.235)/0.191)
                if L < LP:
                    return 0.235 + 0.191*np.log(L)
                x = (np.log(L) - np.log(LP))
                xm = np.log(Lmax) - np.log(LP)
                x = x/xm               #Eq. (17)
                ym = c2M - 0.5
                A1 = 0.191*xm          #Eq. (18)
                A2 = (3*ym - 2*A1)     #Eq. (19)
                A3 = -(2*ym - A1)      #Eq. (20)
                return A1*x + A2*x**2 +A3*x**3 + 0.5
            return c2M
        a, R0 = getContactRadiusAndCurvature(F, s1, s2, x0, largeDef)
        R1 = 1./s1.R - R0
        R2 = 1./s2.R + R0
        L1 = Lambda(a,s1,1./R1)
        L2 = Lambda(a,s2,1./R2)
        #Calculating the indentation depths, Eq. (14)
        h1 = a**2/(2*c2(L1,s1))*R1
        h2 = a**2/(2*c2(L2,s2))*R2
        return h1 + h2,a,R0

    #empty vectors for the contact data
    F = []
    h = []
    a = []
    R = []

    R0 = (1./s2.R - 1./s1.R)/2
    F.append(0)
    h.append(0)
    a.append(0)
    R.append(R0)
    x0 = [0,R0]    #Starting guess for the contact radius and curvature
    Force = 0
    while h[-1]<hMax:
        #iterating until hmax is obtained
        Force = Force+dF
        IndDepth, cRadius, curv = computeIndentationDepth(Force, s1,s2, x0,
                                                          largeDef)
        x0 = [cRadius, curv]
        F.append(Force)
        h.append(IndDepth)
        a.append(cRadius)
        R.append(curv)

    #Converting to numpy arrays
    F = np.array(F)
    h = np.array(h)
    a = np.array(a)
    R = np.array(R)
    return h,F,a,1/R

#==============================================================================
# *** *** Functions for fitting the results according to the appendix *** ***
#==============================================================================
def baseFunc(h, par):
    A1,a1,A2,a2 = par
    return A1*h**a1 + A2*h**a2

#Function to be minimized for the fitting, least squares method
def baseFuncOpt(par, *data):
    f,h = data
    return np.sum((f-baseFunc(h,par))**2)
    
