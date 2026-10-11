
# Test AGHT, OJ, IOJ on the RD problem for Multi-UR-AG, LowMS, RYDE

from sage.all import *
from pulp import * 
import math

def Classical_AGHT(q,m,n,k,r,w): # AGHT - 2018
    WF1 = w*math.log((n-k-1)*m, 2) + (r*math.ceil((k+1)*m/n) - m)*math.log(q,2) 
    WF2 = w*math.log((n-k-1), 2) + math.log(m**2, q) + (r*(k+1))*math.log(q,2) 
    return math.floor(min(WF1, WF2))

def Classical_OJ_Basis(q,m,n,k,r,w): # Basis Enumeration
    N = math.ceil(k*r/(m-r)); WF = w*math.log(N*m, q) + ((m-r)*(r-1))*math.log(q,2)
    return round(WF, 3)

def Classical_OJ_Coordinates(q,m,n,k,r,w): # Coordinates Enumeration
    N = math.ceil((m*(r-1) + k+1-r)/(m-1))
    WF = w*math.log(m*N, 2) + (r-1)*(k+1+N-r)*math.log(q,2)
    return round(WF, 3)

def Classical_OJ(q,m,n,k,r,w): #  Our OJ  Attack
    WF1 = Classical_OJ_Basis(q,m,n,k,r,w); WF2 = Classical_OJ_Coordinates(q,m,n,k,r,w)
    #print(WF1, WF2)
    return math.floor(min(WF1, WF2))

def Classical_IOJ(q,m,n,k,r,w): # IOJ  Attack; Appendix A
    WF1 = w*math.log((k+1)*r, q) + ((m-r)*(r-1))*math.log(q,2)
    if 0 < k+1 <= r: # Case 1
        WF1 = w*math.log(m**(r-1),2)
        return math.floor(min(WF1, WF2))
    elif r < k+1 <= n: # Case 2
        NN = m*(r-1) + k+2-r; WF2 = w*math.log(NN, 2) + (r-1)*(k+2-r)*math.log(q,2)
        return math.floor(min(WF1, WF2))


# Multi-UR-AG: Attack-w
(q,m,n,k,r,w)= (2,97,2*24,24,8,2.81) # 128, AGHT & OJ & IOJ = 225 & 201 & 152
#(q,m,n,k,r,w)= (2,127,2*35,35,9,2.81) # 192, AGHT & OJ & IOJ = 352 & 316 & 252

# Multi-UR-AG: Attack-w1
#(q,m,n,k,r,w)= (2,97,2*24+14,24,8,2.81) # 128, AGHT & OJ & IOJ = 227 & 201 & 152
#(q,m,n,k,r,w)= (2,127,2*35+15,35,9,2.81) # 192, AGHT & OJ & IOJ = 353 & 316 & 252

# LowMS: [n, k, r]
#(q,m,n,k,r,w)= (2,61,2*25,25,7,2.81) # 128, AGHT & OJ & IOJ = 192 & 180 & 144
#(q,m,n,k,r,w)= (2,101,2*37,37,8,2.81) # 192, AGHT & OJ & IOJ = 331 & 293 & 243
#(q,m,n,k,r,w)= (2,101,2*44,44,9,2.81) # 256, AGHT & OJ & IOJ = 400 & 387 & 323

# RYDE-v 2.0.0: [n, k, r],  NIST PQC 2023 [10] &  DCC 2025 [11]
#(q,m,n,k,r,w)= (2,31,33,15,10,2) # 128, AGHT & OJ & IOJ = 147 & 160 & 79
#(q,m,n,k,r,w)= (2,37,41,18,13,2) # 192, AGHT & OJ & IOJ = 216 & 245 & 101
#(q,m,n,k,r,w)= (2,43,47,18,17,2) # 256, AGHT & OJ & IOJ = 283 & 323 & 66

# RYDE-v 2.1.0: [n, k, r], NIST PQC 2025 [12]
#(q, m, n, k, r, w)= (2, 53, 53, 45, 4, 2) # 128, AGHT & OJ & IOJ = 148 & 153 & 144
#(q, m, n, k, r, w)= (2, 61, 61, 51, 5, 2) # 192, AGHT & OJ & IOJ = 217 & 224 & 208
#(q, m, n, k, r, w)= (2, 61,57, 44, 7, 2) # 256, AGHT & OJ & IOJ = 301 & 287 & 251

AGHT = Classical_AGHT(q,m,n,k,r,w)
OJ = Classical_OJ(q,m,n,k,r,w)
IOJ = Classical_IOJ(q,m,n,k,r,w)

print(AGHT, "&", OJ, "&", IOJ) 
