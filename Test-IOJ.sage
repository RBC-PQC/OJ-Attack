from sage.all import *
import random

def random_small_vector_generator(n,t):
    S = matrix(Fqm.base_ring(), m, t, 0)
    while S.rank() != t:
        S = random_matrix(Fqm.base_ring(), m, t)
    C = matrix(Fqm.base_ring(), m, n, 0)
    while C.rank() != t:
        C = S * random_matrix(Fqm.base_ring(), t, n) 
    A =  C.transpose()   
    return vector(Fqm, [A[i] for i in range(n)])  # A[i] : the first i row of A 

def vector_matrix(small_vector):    # To  m * length  matrix
    length = len(list(small_vector))
    return matrix(Fqm.base_ring(),length, m, [vector(small_vector[j]) for j in range(length)])

def vector_space(small_vector):    # row_space; support
    return vector_matrix(small_vector).row_space()

def vector_to_basis_list(small_vector):    # basis
    return [Fqm(x) for x in vector_space(small_vector).basis_matrix().rows()]

def Rank_Submatrix(error):
    r = len(vector_to_basis_list(error))
    S = vector_matrix(vector_to_basis_list(error)).transpose() # m * r  matrix
    SC =  vector_matrix(error).transpose()  # m * n  matrix
    C = S.solve_right(SC) # r * n  matrix
    C_Row_Space_Basis_Matrix = C.row_space().basis_matrix()
    sub_C_Row_Space_Basis_Matrix = C_Row_Space_Basis_Matrix.matrix_from_columns(list(range(r)))
    return sub_C_Row_Space_Basis_Matrix.rank()

def test(totalltests):
    succ = 0
    failure = 0
    G = random_matrix(Fqm, k, n)
    while G.rank() != k:
        G = random_matrix(Fqm, k, n)
    for npair in range(totalltests):
        e = random_small_vector_generator(n, r)
        y = Message*G + e
        P = identity_matrix(Fqm.base_ring(), n)
        while Rank_Submatrix(e) != r:
            SG = SymmetricGroup(n)
            P = SG.random_element().matrix()
            e  = e*P
        # Construct instance （yP, GP, eP）such that yp = Message*GP + ep
        y = y*P;  G = G*P;  G_y = block_matrix(2, 1, [matrix(y), G])
        Systematic_G_y = G_y.echelon_form()
        Identity_Matrix  = identity_matrix(Fqm.base_ring(), k+1, k+1)
        
        R = Systematic_G_y.matrix_from_columns(Columns_Index_1)  
        S = vector_matrix(vector_to_basis_list(e)).transpose()  # m * r  matrix
        SC =  vector_matrix(e).transpose()   # n * m  matrix
        C = S.solve_right(SC)   # r * n  matrix

        A1 = C.matrix_from_rows_and_columns(Rows_Index, Columns_Index)  
        A2 = C.matrix_from_rows_and_columns(Rows_Index, list(range(k+1, n))) 
        j = random.choice(list(range(n-k-1)))
        aj = A2.matrix_from_columns([j]); rj = R.matrix_from_columns([j]) 
        A1_aj_matrix = block_matrix(1, 2, [A1, aj])     # [A1  aj] matrix
        list_rj_one = rj.list() + list(vector(Fqm, 1, [-1]))  # [rj  -1]
        Tj = vector_matrix(list_rj_one)    # (k+2) * m matrix
        try:
            if (S * A1_aj_matrix * Tj == zero_matrix(Fqm.base_ring(), m, m)): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f" % (succ,totalltests,succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" % (failure,totalltests,failure/totalltests))


(q,m,n,k,r) = (2, 31, 30, 15, 9)  

Fqm = GF(q**m)

G = random_matrix(Fqm, k, n)
while G.rank() != k:
    G = random_matrix(Fqm, k, n)

Message = random_vector(Fqm, k)  


Rows_Index = [i for i in range(0, r)]
Columns_Index = [i for i in range(0, k+1)]
Columns_Index_1 = [i for i in range(k+1, n)]


test(100)
