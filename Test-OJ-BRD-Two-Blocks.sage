# Test AGHT, PRR, OJ on the BRD problem for BRQC, BOuroboros-R
# Two blocks

from sage.all import *
from pulp import * 
from math import *
from tabulate import tabulate

def Classical_AGHT(q,m,n,k,r,n1,n2,r1,r2,w): # AGHT - 2025
    WF1 = w*math.log((n-k-1)*m, 2) + (r*math.ceil((k+1)*m/n) - m)*math.log(q,2) 
    WF2 = w*math.log((n-k-1), 2) + math.log(m**2, q) + (r*(k+1))*math.log(q,2)
    return round(min(WF1, WF2), 3)

def Classical_PRR1(q, m, n, k, r, n1, n2, r1, r2, w): # PRR - 2025
    prob1 = LpProblem('PRR1:', sense=LpMinimize) # sense: LpMinimize or LpMaximize
    t1 = LpVariable('t1', r1, m, LpInteger); t2 = LpVariable('t2', r2, m, LpInteger)
    prob1 += (m-t1)*r1 + (m-t2)*r2 
    prob1 += t1 + t2 <= m;  prob1 += t1*n1 + t2*n2 <= m * (n-k-1) 
    status1 = prob1.solve(PULP_CBC_CMD(msg=False))  
    WF = (value(prob1.objective) - m) * math.log(q,2) + w*math.log((n-k-1)*m,2)
    return round(WF,3)

def Classical_PRR2(q, m, n, k, r, n1, n2, r1, r2, w): # PRR - 2025
    prob2 = LpProblem('PRR2:', sense=LpMinimize)
    t1 = LpVariable('t1', r1, n1, LpInteger); t2 = LpVariable('t2', r2, n2, LpInteger)
    prob2 += (n1-t1)*r1 + (n2-t2)*r2 
    prob2 += t1 + t2 <= n-k-1  
    status2 = prob2.solve(PULP_CBC_CMD(msg=False)) 
    WF = value(prob2.objective) * math.log(q,2) + w*math.log(n-k-1,2) + math.log(m**2,2)
    return round(WF,3)

def Classical_PRR(q, m, n, k, r, n1, n2, r1, r2, w): # PRR - 2025
    WF1 = Classical_PRR1(q, m, n, k, r, n1, n2, r1, r2, w)
    WF2 = Classical_PRR2(q, m, n, k, r, n1, n2, r1, r2, w)
    return round(min(WF1, WF2), 3)

def Classical_OJ_Basis(q, m, n, k, r, n1, n2, r1, r2, w): # Basis Enumeration
    N = math.ceil(k*r/(m-r)); WF = w*math.log(N*m, 2) + ((m-r)*(r-1))*math.log(q,2)
    return round(WF, 3)

def Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, r1, r2, w): # Coordinates Enumeration 
    gamma = min(r1, r2); N = math.ceil((m*(r-1) + (n1-r1))/m); T = (r1-1)*(n1-r1)
    if 0 < k+1 <= r1: # Case 1
       WF1 = w*math.log(m*N,2); WF2 = (N-r1-r2+k+1)*gamma * math.log(q,2)
       return round(WF1 + WF2, 3)
    elif r1 < k+1 <= n1: # Case 2
       WF1 = w*math.log(m*N,2); WF2 = ((r1-1)*(k+1-r1) + (N-r2)*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1 < k+1 <= n1+r2: # Case 3
       WF1 = w*math.log(m*N,2); WF2 = (T + (N-n1-r2+k+1)*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)
    else: # Case 4: n1+r2 < k+1 <= n1+n2
       WF1 = w*math.log(m*N,2); WF2 = (T + (k+1-n1-r2)*r2 + N*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)

def Classical_OJ(q, m, n, k, r, n1, n2, r1, r2, w):
    WF1 = Classical_OJ_Basis(q, m, n, k, r, n1, n2, r1, r2, w)
    WF2 = Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, r1, r2, w)
    return round(min(WF1, WF2), 3)


# BRQC: Attack-2n
#(q,m,n,k,r,n1,n2,r1,r2,w)= (2,83,79*2,79,8,79,79,4,4,2.81) 
# 128   296 & 292 & 271
#(q,m,n,k,r,n1,n2,r1,r2,w)= (2,127,113*2,113,10,113,113,5,5,2.81)  
# 192   561 & 556 & 490
(q,m,n,k,r,n1,n2,r1,r2,w)= (2,139,137*2,137,10,137,137,5,5,2.81)  
# 256   610 & 605 & 587


# BOuroboros-R: Attack-2n
#(q,m,n,k,r,n1,n2,r1,r2,w)= (2,59,73*2,73,8,73,73,4,4,2.81) 
# 128   214 & 214 & 256
(q,m,n,k,r,n1,n2,r1,r2,w)= (2,89,101*2,101,9,101,101,4,5,2.81)  
# 192   352 & 312 & 342 
#(q,m,n,k,r,n1,n2,r1,r2,w)= (2,97,109*2,109,10,109,109,5,5,2.81)  
# 256   430 & 430 & 479


AGHT = Classical_AGHT(q, m, n, k, r, n1, n2, r1, r2, w)
PRR = Classical_PRR(q, m, n, k, r, n1, n2, r1, r2, w)
OJ = Classical_OJ(q, m, n, k, r, n1, n2, r1, r2, w)

headers = ["Attacks", "Classical Bit Complexity"]
data = [["AGHT attack ", AGHT], ["PRR attack  ", PRR], ["OJ attack  ", OJ]]
table =  tabulate(data, headers=headers, tablefmt="grid", numalign="center", stralign="center", showindex=True)
print(table) 


print(floor(AGHT), "&", floor(PRR), "&", floor(OJ))
