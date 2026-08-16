# Test AGHT, NH-AGHT, OJ for NH-Multi-RQC-AG and NH-Multi-UR-AG

from sage.all import *
from pulp import * 
import math
from tabulate import tabulate

def Classical_AGHT(q, m, n, k, r, n1, n2, n3, r1, r2, w): # AGHT 
    WF1 = w*math.log((n-k-1)*m,2) + (r*math.ceil((k+1)*m/n) - m) * math.log(q,2) 
    WF2 = w*math.log((n-k-1),2) + math.log(m**2,2) + (r*(k+1)) * math.log(q,2) 
    return round(min(WF1, WF2), 3)

def Classical_AGHT_NHRD(q, m, n, k, r, n1, n2, n3, r1, r2, w):  # NH-AGHT
    prob = LpProblem('NH-AGHT:', sense=LpMinimize)  # sense: LpMinimize or LpMaximize
    t1 = LpVariable('t1', r1, m, LpInteger); t2 = LpVariable('t2', r2, m, LpInteger)
    prob += (m-t1)*(r1 + r2) - t2*r2 # Objective
    prob += t1 + t2 <= m;   prob += t1*(2*n1+n2) + t2*n2 <= m*(n1+n2-1)  # constraints
    status = prob.solve(PULP_CBC_CMD(msg=False))  # solving
    diction = {v.name: v.varValue for v in prob.variables()}
    t1 = diction['t1']; t2 = diction['t2']
    WF = (value(prob.objective) - m) * math.log(q,2) + w*math.log(m*(n1+n2-1)-n2*(t1+t2), 2)
    print(WF, {v.name: v.varValue for v in prob.variables()})
    return round(WF,3)

def Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, n3, r1, r2, w): # Coordinates Enumeration
    N = math.ceil((m*(r-1) + (k+1-r1))/(m-1))
    if 0 < k+1 <= n1+r2: # Case 1
       WF1 = w*math.log(m*N,2); WF2 = (N+k+1-r1)*(r1-1) + math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1+r2 < k+1 <= n1+n2: # Case 2
       WF1 = w*math.log(m*N,2); WF2 = (N+k+1-r1)*(r1-1) + r2*(k+1-r2) + math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n1+n2  < k+1 <= n1+n2+n3: # Case 3
       WF1 = w*math.log(m*N,2); WF2 = (N+k+1-r1)*(r1-1) + r2*(n2-r2) + math.log(q,2)
       return round(WF1 + WF2, 3)

def Classical_OJ_Coordinates_Permuted(q, m, n, k, r, n1, n2, n3, r1, r2, w): # Coordinates Enumeration
    N = math.ceil((m*(r-1) + (k+1-r))/m)
    if 0 < k+1 <= r: # Case 1
       WF1 = w*math.log(m*N,2); WF2 = (N-r+k+1)*(r1-1) + math.log(q,2)
       return round(WF1 + WF2, 3)
    elif r < k+1 <= n2: # Case 2
       WF1 = w*math.log(m*N,2); WF2 = (r-1)*(k+1-r1) + N*(r1-1) + math.log(q,2)
       return round(WF1 + WF2, 3)
    elif n2  < k+1 <= n1+n2+n3: # Case 3
       WF1 = w*math.log(m*N,2); WF2 = (r-1)*(n2-r) + (k+1-n2)*(r1-1) + N*(r1-1) + math.log(q,2)
       return round(WF1 + WF2, 3)

# NH-Multi-RQC-AG
(q,m,n,k,r,n1,n2,n3,r1,r2,w)= (2,61,3*50,50,12,50,50,50,7,5,2.81) # 128  (e1,e2,e3) 
#(q,m,n,k,r,n1,n2,n3,r1,r2,w)= (2,79,3*95,95,13,95,95,95,8,5,2.81) # 192  (e1,e2,e3)

# NH-Multi-UR-AG
#(q,m,n,k,r,n1,n2,n3,r1,r2,w)= (2,73,2*22+13,22,12,22,13,22,8,4,2.81) # 128  (e1,e2,e3) 
# (q,m,n,k,r,n1,n2,n3,r1,r2,w)= (2,97,2*30+14,30,13,30,14,30,9,4,2.81) # 192  (e1,e2,e3)


AGHT = Classical_AGHT(q, m, n, k, r, n1, n2, n3, r1, r2, w)
AGHT_NHRD = Classical_AGHT_NHRD(q, m, n, k, r, n1, n2, n3, r1, r2, w)
OJ_NHRD = Classical_OJ_Coordinates(q, m, n, k, r, n1, n2, n3, r1, r2, w)
OJ_NHRD_Permuted = Classical_OJ_Coordinates_Permuted(q, m, n, k, r, n1, n2, n3, r1, r2, w)

headers = ["Attacks", "Classical Bit Complexity"]
data = [["AGHT attack ", AGHT], ["AGHT_NHRD attack ", AGHT_NHRD], ["OJ_NHRD attack  ", OJ_NHRD], ["OJ_NHRD_Permuted attack  ", OJ_NHRD_Permuted]]
table =  tabulate(data, headers=headers, tablefmt="grid", numalign="center", stralign="center", showindex=True)
print(table)
