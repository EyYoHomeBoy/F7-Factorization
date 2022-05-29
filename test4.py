import math
from time import time
import numpy as np
from math import gcd
from sympy import legendre_symbol

def TrialDiv(n, P):
    '''
    Input: En int n>0, en liste af primfaktorer
    Output: En liste indeholdende eksponenten af hver primfaktor modulo 2
    '''
    gamma = []
    for p in P:
        if n == 1:
            gamma.append(0)
        else:
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            gamma.append(k%2)
    return n,gamma

def intSqrt(N):
    '''
    Input: En int N>0
    Output: heltalsværdien af sqrt(N)
    '''
    l = 0
    h = N + 1
    while l + 1 < h:
        mid = (h + l) // 2
        if mid**2 > N:
            h = mid
        else:
            l = mid
    return l
 
def AQpair(k,N,B):
    '''
    Input: Ints k,N,B > 0
    Out: To lister bestående af værdierne Q_n og A_n genereret ud fra sqrt(k*N)
    for n <= B
    '''
    g = intSqrt(k*N)
    Lq = []
    LP = []
    index = [0]
    LQ=[k*N,1]
    LA = [0,1]
    P = 0
    LP.append(P)
    for i in range(1,B):
        index.append(index[i-1]+1)
        q = (g+P)//LQ[i]
        Lq.append(q)
        Pformer = P
        P = q*LQ[i] - Pformer
        LP.append(P)
        LQ.append( LQ[i-1] - q * (P-Pformer) )
        LA.append(( q * LA[i] + LA[i-1])%N )
    return LQ[2:], LA[2:]

def primesInRange(n):
    '''
    Input: En int n>0
    Output: En liste af alle primtal op til n
    '''
    prime = [True for i in range(n+1)]
    primes = []
    p = 2
    while (p**2 <= n):
        if (prime[p] == True):
            for i in range(p**2, n+1, p):
                prime[i] = False
        p += 1
    for p in range(2, n+1):
        if prime[p]:
            primes.append(p)
    return primes

def Factorbase(k,B,N):
    '''
    Input: Ints k,B,N > 0
    Output: En liste bestående af de primtal p, hvorom det gælder, at (k*N/p)=1 eller (k*N/p)= 0
    '''
    P = primesInRange(B)[1 : ]
    aList = []
    aList.append(2)
    for p in P:
        if legendre_symbol(k*N,p) == 1:
            aList.append(p)
        if legendre_symbol(k*N,p) == 0:
            aList.append(p)
    return aList

def Qfactor(Qs, As ,base, UB, B1= intSqrt(2**128+1)/500,B2 = intSqrt(2**128+1)/2*10**7):
    '''
    Input: Lister Qs, As, base, som indeholder hhv. værdierne for Q_n ,A_{n-1}
    samt primtal i en faktorbase. Yderligere tages som input en int UB svarende
    til den øvre grænse U_B beskrevet i afsnit 2.1 og floats B1, B2 svarende
    til de øvre grænser der indgår i Early Stopping.
    Output: Et np.array med shape = (T,S+1), hvor T er antal Q_n'er som
    faktoriseres og S er antallet af primtal der bliver faktoriseret over.
    Yderligere to np.arrays med shape = (1,T) indeholdende værdierne for de
    Q_n'er som bliver faktoriseret og de tilhørende A_{n-1}-værdier.
    '''
    goodQ = []
    goodA = []
    mat = []
    baseExt  = []
    for i in range(len(Qs)):
        gamma = []
        if i%2 == 0:
            gamma.append(1)
        else:
            gamma.append(0)
        n , gamma = TrialDiv(Qs[i],base[:15])[0], gamma+TrialDiv(Qs[i],base[:15])[1]
        if n < B1:
            s = TrialDiv(n,base[15:80])
            n = s[0]
            gamma += s[1]
            if n < B2:
                s = TrialDiv(n,base[80:])
                n , gamma = s[0] , gamma + s[1]     
        if n < UB:
            if n not in set(baseExt) and n > 1:
                baseExt.append(n)
                for row in mat:
                    row.append(0)
            for p in baseExt:
                k = 0
                if n == 1:
                    gamma.append(k)
                else:
                    while n%p == 0:
                        n //= p
                        k += 1
                    gamma.append(k%2)
            mat.append(gamma)
            goodQ.append(Qs[i])
            goodA.append(As[i])
        if i % 100000 == 0:
            print(i)
    return np.array(mat), np.array(goodQ), np.array(goodA)    

def Qsqrt(L,N):
    '''
    Input: En liste af ints L, sådan at produktet er et kvadrattal,
    og en int N
    Output: Kvadratroden af produktet af elementerne i L modulo N 
    '''
    Q = L[0]
    SR = 1
    for i in range(len(L)-1):
        G = gcd(int(Q),int(L[i+1]))
        Q = ((Q // G)*(L[i+1] // G))
        SR *= G
        SR = SR % N
    return int(SR * intSqrt(Q)) % N

def isPivot(V,p):
    '''
    Input: Et array V indeholdende 1 og 0, en int p>0
    Output: True vis V[p] er det yderst højreliggende 1-tal, False ellers.
    '''
    if p == len(V):
        return True
    if V[p] == 1:
        for i in range(p+1,len(V)):
            if V[i] == 1:
                return False
        return True

def GaussElim(M):
    '''    
    Input: Et np.array M
    Output: np.arrays svarende til M' og H' som beskrevet i afsnit 2.1
    '''
    m = len(M[1]) - 1
    n = len(M)
    j = m
    i = 0
    row = n
    H = np.identity(n)
    while j >= 0:
        c = False
        while c == False and i < n:
            if M[i][j]== 1 and isPivot(M[i],j):
                c  = True
                row = i
            i += 1
        if row < n-1:
            for l in range(row + 1,n):
                if M[l][j] == 1 and isPivot(M[l],j):
                    M[l] = np.mod(M[l]+M[row],2)
                    H[l] = np.mod(H[l]+H[row],2)
        j -= 1
        i = 0
        row = n
    return M,H

def FindFactor(M,H,Qvect,Avect,N):
    '''
    Input: Et np.array M med shape = (T,S+1), Et np.array H med shape = (T,T),
    np.arrays Qvect, Avect af længde T indeholdende hhv. værdier for Q_n og
    A_{n-1} og en int N
    Output: To ikke-trivielle faktorer af N ellers printes:
    "Ingen ikke-trivielle faktorer fundet"
    '''
    n = len(M)
    j = len(M[1]) - 1
    for i in range(n):
        if np.array_equiv(M[i],np.zeros((1,j+1))):
            L =[]
            As =[]
            for k in range(n):
                if H[i][k] == 1:
                    L.append(Qvect[k])
                    As.append(Avect[k])
            Q = Qsqrt(L,N)
            A = 1
            for a in As:
                A = A * a%N
            if A != Q and gcd(N, A - Q) != 1 and A - Q != N:
                return gcd(N, A-Q), gcd(N,A+Q) 
    return "Ingen ikke-trivielle faktorer fundet"

total_time = time()
t = time()
N = 2**128+1
k = 257
UB = 52000
var = AQpair(k, N,1150000)
Qs = var[0]
As = var[1]
F = Factorbase(k,UB,N)
base = F[:800]
print(time()-t)
t = time()
D = Qfactor(Qs, As, base, UB)
print("Eksponent er lavet")
print(time()-t)
t = time()
Gamma = D[0]
Q = D[1]
A = D[2]
K = GaussElim(Gamma)
M = K[0]
H = K[1]
print(FindFactor(M,H,Q,A,N))
print(time()-t)
print(time()-total_time)
