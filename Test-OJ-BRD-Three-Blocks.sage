# Test AGHT, PRR, OJ on the BRD problem for BRQC, BOuroboros-R
# Three blocks

from sage.all import *
from pulp import * 
import math
from tabulate import tabulate

def Classical_AGHT(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w): # AGHT - 2025
    WF1 = w*math.log((n-k-1)*m,2) + (r*math.ceil((k+1)*m/n) - m) * math.log(q,2) 
    WF2 = w*math.log((n-k-1),2) + math.log(m**2,2) + (r*(k+1)) * math.log(q,2) 
    return round(min(WF1, WF2), 3)

def Classical_PRR1(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w):  # PRR-1 2025
    prob = LpProblem('PRR1:', sense=LpMinimize)  # sense: LpMinimize or LpMaximize
    t1 = LpVariable('t1', r1, m, LpInteger); t2 = LpVariable('t2', r2, m, LpInteger)
    t3 = LpVariable('t3', r3, m, LpInteger)
    prob += (m-t1)*r1 + (m-t2)*r2 + (m-t3)*r3 
    prob += t1 + t2 + t3 <= m;   prob += t1*n1 + t2*n2 + t3*n3 <= m*(n-k-1)  
    status = prob.solve(PULP_CBC_CMD(msg=False))  
    #print(value(prob.objective), {v.name: v.varValue for v in prob.variables()})
    WF = (value(prob.objective) - m) * math.log(q,2)  + w*math.log((n-k-1)*m, 2)
    #print(WF, {v.name: v.varValue for v in prob.variables()})
    return round(WF,3)

def Classical_PRR2(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w): # PRR-2 2025
    prob = LpProblem('PRR2:', sense=LpMinimize)
    t1 = LpVariable('t1', r1, n1, LpInteger); t2 = LpVariable('t2', r2, n2, LpInteger)
    t3 = LpVariable('t3', r3, n3, LpInteger)
    prob += (n1-t1)*r1 + (n2-t2)*r2 + (n3-t3)*r3  
    prob += t1 + t2 + t3 <= n-k-1  
    status = prob.solve(PULP_CBC_CMD(msg=False)) 
    #print(value(prob.objective), {v.name: v.varValue for v in prob.variables()})
    WF = value(prob.objective) * math.log(q,2) + w*math.log(n-k-1,2) + math.log(m**2,2)
    #print(WF, {v.name: v.varValue for v in prob.variables()})
    return round(WF,3)

def Classical_PRR(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w):
    WF1 = Classical_PRR1(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
    WF2 = Classical_PRR2(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
    return round(min(WF1, WF2), 3)

def Classical_OJ_Basis(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w): # Basis Enumeration
    N = math.ceil(k*r/(m-r)); WF = w*math.log(N*m, q) + ((m-r)*(r-1)) * math.log(q,2) 
    return round(WF, 3)

def Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w): # Coordinates Enumeration
    gamma = min(r1, r2, r3); N = math.ceil((m*(r-1) + (n1-r1))/m); T = (r1-1)*(n1-r1)
    if 0 < k+1 <= r1: # Case 1
       WF1 = w*math.log(m*N,2); WF2 = (N-r1-r2-r3+k+1)*gamma * math.log(q,2)
       return round(WF1 + WF2, 3)
    elif r1 < k+1 <= n1: # Case 2
       WF1 = w*math.log(m*N,2); WF2 = ((r1-1)*(k+1-r1) + (N-r2-r3)*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1 < k+1 <= n1+r2: # Case 3
       WF1 = w*math.log(m*N,2); WF2 = (T + (N-n1-r2-r3+k+1)*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1+r2 < k+1 <= n1+n2: # Case 4
       WF1 = w*math.log(m*N,2); WF2 = T + (k+1-n1-r2)*r2 + (N-r3)*gamma + math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1+n2 < k+1 <= n1+n2+r3: # Case 5
       WF1 = w*math.log(m*N,2); WF2 = (T + (n2-r2)*r2 + (N-n1-n2-r3+k+1)*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)
    else: # Case 6: n1+n2+r3 < k+1 <= n1+n2+n3
       WF1 = w*math.log(m*N,2); WF2 = (T + (n2-r2)*r2 + (k+1-n1-n2-r3)*r3 + N*gamma) * math.log(q,2)
       return round(WF1 + WF2, 3)

def Classical_OJ(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w):
    WF1 = Classical_OJ_Basis(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
    WF2 = Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
    return round(min(WF1, WF2), 3)


# BRQC: Attack-3n
# (q,m,n,k,r,n1,n2,n3,r1,r2,r3,w) = (2,83,3*79,79,12,79,79,79,4,4,4,2.81) 
# 128, 303 & 353 & 272
# (q,m,n,k,r,n1,n2,n3,r1,r2,r3,w) = (2,127,3*113,113,15,113,113,113,5,5,5,2.81) 
# 192, 559 & 605 & 492
(q,m,n,k,r,n1,n2,n3,r1,r2,r3,w) = (2,139,3*137,137,19,137,137,137,6,7,6,2.81)
# 256, 796 & 864 & 728

# BOurobors: : Attack-3n
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3,w) = (2,59,3*73,73,13,73,73,73,4,4,5,2.81) 
# 128   237 & 327 & 258 
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3,w) = (2,89,3*101,101,18,101,101,101,6,6,6,2.81) 
# 192   490 & 646 & 553 
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3,w)= (2,97,3*109,109,19,109,109,109,6,6,7,2.81) 
# 256   570 & 695 & 593 


AGHT = Classical_AGHT(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
PRR = Classical_PRR(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)
OJ = Classical_OJ(q, m, n, k, r, n1, n2, n3, r1, r2, r3, w)

headers = ["Attacks", "Classical Bit Complexity"]
data = [["AGHT attack ", AGHT], ["PRR attack  ", PRR], ["OJ attack  ", OJ]]
table =  tabulate(data, headers=headers, tablefmt="grid", numalign="center", stralign="center", showindex=True)
print(table)

print(floor(AGHT), "&", floor(PRR), "&", floor(OJ))
