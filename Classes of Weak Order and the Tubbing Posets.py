#!/usr/bin/env python
# coding: utf-8



load("Unit_Interval_graphs.py")
QSym = QuasiSymmetricFunctions(QQ)
YQS = QSym.YQS()
Fund = QSym.F()
QS = QSym.QS()
QQ = sage.rings.rational_field.RationalField()
schur = SymmetricFunctions(QQ).schur()



def is_cut_vertex_set(g,u,v,L):
    """
    Given a graph g, vertices u and v on graph g and a subset of vertices L of g that 
    do not contain u and v. We see if g minus the vertices in L still has a path from u to v,
    i.e. checking that L is a uv-cut set in g. 
    """
    h = g.copy()
    h.delete_vertices(L)
    if len(h.shortest_path(u,v))==0:
        return True
    return False




def are_related(pi,sigma,g):
    """
    Pi and sigma are permutations of {1,2,3,...,n}.
    We determine if the two are related according to g. 
    This means that they are distance 1 in the weak order and 
    at the two places pi and sigma differ, the numbers to the 
    right make a cut set. 
    """
    n = len(pi)
    equal = 0
    a = 0
    for i in range(n):
        if pi[i]==sigma[i]:
            equal += 1
        else:
            a = i
    if not equal == n-2:
        return False
    if pi[a-1]==sigma[a-1]:
        return False
    return is_cut_vertex_set(g,pi[a-1],pi[a],pi[a+1:])






def print_equivalence_classes_Bruhat(g):
    """
    Given a graph g (it should be a unit interval graph)
    we build the equivalence classes from the identity up. 
    More or less we go layer by layer of the weak order starting with the identity, 
    and consider what is related between each layer. 
    We return PI, the list of all permutations of 1,2,...,n
    and C, which is a tuples n! long. A matching number in the tuple implies
    the associated permutations are in the same equivalence class. 
    """
    C = []
    PI = []
    n = g.order()
    for l in range(n*(n-1)/2+1):#focusing on permutations of length l
        pi = [i+1 for i in range(n)]
        while (not pi ==False):
            if Permutation(pi).length()==l:
                #print("pi:",pi)
                found = False
                for j in range(len(PI)):
                    sigma = PI[j]
                    if found == False and are_related(pi,sigma,g):
                        C += [C[j]]
                        found = True
                if found == False:
                    if C == []:
                        C += [0]
                    else:
                        C += [max(C)+1]
                PI += [pi]
            pi = Permutation(pi).next()
    return [C,PI]





def gather_classes(g):
    """
    We reorder the information from print_equivalence_classes_Bruhat(g)
    so that permutations in the same class are gathered. 
    """
    [C,PI] = print_equivalence_classes_Bruhat(g)
    CC = []
    for i in range(len(C)):
        if C[i]>=len(CC):
            CC += [[PI[i]]]
        else:
            CC[C[i]]+=[PI[i]]
    return CC





def is_max(pi):
    """
    Return True if pi is the maximum element of its class. 
    """
    n = len(pi)
    for i in range(n-1):
        if pi[i]<pi[i+1] and is_cut_vertex_set(g,pi[i],pi[i+1],pi[i+2:]):
            return False
    return True
def is_min(pi):
    """
    Return True if pi is the minimum element of its class. 
    """
    n = len(pi)
    for i in range(n-1):
        if pi[i]>pi[i+1] and is_cut_vertex_set(g,pi[i],pi[i+1],pi[i+2:]):
            return False
    return True
def mins_of_classes(g):
    """
    We will find permutations that are minimum elements of the classes. 
    """
    mins = []
    CC = gather_classes(g)
    for C in CC:
        for c in C:
            if is_min(c)==True:
                mins += [c]
    return mins
def maxs_of_classes(g):
    """
    We will find permutations that are maximum elements of the classes. 
    """
    maxs = []
    CC = gather_classes(g)
    for C in CC:
        for c in C:
            if is_max(c)==True:
                maxs += [c]
    return maxs





def print_classes_nicely(g):
    """
    We will print the information about about the equivaance classes
    given a unit interval graph g. 
    """
    CC = gather_classes(g)
    for C in CC:
        print(C)
        for c in C:
            if is_max(c):
                print("max:", c)
        for c in C:
            if is_min(c):
                print("min:", c)
        print()





def tubing_poset_mins(g):
    """
    We form the tubbing poset associated to graph g. 
    The vertices are all written as the minimum elements of the associated class. 
    """
    V = mins_of_classes(g)
    #print(V)
    E = []
    n = len(V)
    elm_labs = {}
    for i in range(n):
        elm_labs.update({i:str(V[i])})
        for j in range(i+1,n):
            if Permutation(V[i]).permutohedron_lequal(Permutation(V[j])):
                #print(V[i], "<", V[j])
                E += [[i,j]]
    #print(E)
    VV = [i for i in range(len(V))]
    return Poset((VV,E)).relabel(elm_labs)
def tubing_poset_maxs(g):
    """
    We form the tubbing poset associated to graph g. 
    The vertices are all written as the maximum elements of the associated class. 
    """
    V = maxs_of_classes(g)
    #print(V)
    E = []
    n = len(V)
    elm_labs = {}
    for i in range(n):
        elm_labs.update({i:str(V[i])})
        for j in range(i+1,n):
            if Permutation(V[i]).permutohedron_lequal(Permutation(V[j])):
                #print(V[i], "<", V[j])
                E += [[i,j]]
    #print(E)
    VV = [i for i in range(len(V))]
    return Poset((VV,E)).relabel(elm_labs)






def print_all_classes_all_graphs(n):
    """
    For all unit interval graphs we will print out what the classes are
    noting specifically the min and max of each class. 
    """
    for g in all_unit_interval_graphs_connected(n):
        g.show(figsize=2)
        print_classes_nicely(g)
        print("_______________________________________________")
    return 0





def string_to_list(string):
    """
    Since we are feeding in the chains as a sequence of permutations that are stored as strings, 
    we will need to be able to turn that string back into a list. This method does that. 
    """
    L = string[1:-1].split(",")
    return [int(a) for a in L]


def differences(A,B):
    """
    Supposing A and B have the same length we return all indices (they start at 1)
    where A and B are different. 
    """
    return [i+1 for i in range(len(A)) if not A[i]==B[i]]
        

def cycle_sequence(C):
    """
    We input C, a list of strings where each string is a permutation. 
    We turn that string into a list. 
    We then look at adjacent permutations and determine
    where they take on different values. Because of our structure on the poset
    This is esentially the cycle that brings you from one permutation to the next. 
    We return all the differences as a list. 
    """
    C_lists = [string_to_list(c) for c in C]
    R = []
    for i in range(len(C)-1):
        R += [differences(C_lists[i],C_lists[i+1])]
    return R



def composition_of_chain_from_minweakdec(cycles):
    """ 
    Given a list of cycles, we return a composition associated to 
    where the smallest elements in the cycles have weak decreases. 
    """
    comp = []
    c = 0
    for i in range(len(cycles)-1):
        c += 1
        if cycles[i][0]>=cycles[i+1][0]:
            comp = comp + [c]
            c = 0
    c += 1
    comp = comp + [c]
    return comp

def composition_of_chain_from_maxweakdec(cycles):
    """
    Given a list of cycles, we return a composition associated to 
    where the largest elements in the cycles have weak decreases. 
    """
    comp = []
    c = 0
    for i in range(len(cycles)-1):
        c += 1
        if cycles[i][-1]>=cycles[i+1][-1]:
            comp = comp + [c]
            c = 0
    c += 1
    comp = comp + [c]
    return comp


def fundamental_quasisymmetric_of_tubing_maxreps_from_minweakdec(g):
    """
    Using the compositions from composition_of_chain_from_minweakdec(cycles), 
    we construct a quasisymmetric function in the fundamental basis. 
    """
    chains_as_list_of_str = tubing_poset_maxs(g).maximal_chains()
    f = 0
    for C in chains_as_list_of_str:
        comp = composition_of_chain_from_minweakdec(cycle_sequence(C))#change this to get different versions
        f = f + Fund(comp)
    return f


#this is the one we use in our paper (that we did in calculations, but I noticed a typo in the definition of Des). 
#vertices are maximums of equivalence classes
#descents come from looking at the top of each cycle and seeing if we have a weak decrease
def fundamental_quasisymmetric_of_tubing_maxreps_from_maxweakdec(g):
    """
    Using the compositions from composition_of_chain_from_maxweakdec(cycles), 
    we construct a quasisymmetric function in the fundamental basis.  
    """
    chains_as_list_of_str = tubing_poset_maxs(g).maximal_chains()
    f = 0
    for C in chains_as_list_of_str:
        comp = composition_of_chain_from_maxweakdec(cycle_sequence(C))#change this to get different versions
        f = f + Fund(comp)
    return f
        







def print_F_per_degree(f,n):
    """
    This code  will nicely separate your quasisymmetric function by degree in the fundamental basis.
    Note that n is the homogeneous degree of f, and this is likely to be n*(n-1)/2 where
    n is the number of vertices of the graph g. 
    """
    g = Fund(f)
    per_degree = [0 for i in range(n+1)]
    total_coeff = [0 for i in range(n+1)]
    for i in range(n+1):
        for A in Compositions(i):
            c = g.coefficient(A)
            per_degree[i] = per_degree[i] + c*Fund(A)
            total_coeff[i] = total_coeff[i] + c
    for i in range(n+1):
        print("Degree", i,": (coeff total =",total_coeff[i], ") ", per_degree[i])
def print_YQS_per_degree(f,n):
    """
    This code  will nicely separate your quasisymmetric function by degree in the YQS basis
    and check whether it is also a symmetric function. 
    Note that n is the homogeneous degree of f, and this is likely to be n*(n-1)/2 where
    n is the number of vertices of the graph g. 
    """
    g = YQS(f)
    per_degree = [0 for i in range(n+1)]
    for i in range(n+1):
        for A in Compositions(i):
            c = g.coefficient(A)
            per_degree[i] = per_degree[i] + c*YQS(A)
    for i in range(n+1):
        print("Degree", i,": ", per_degree[i])
        if not per_degree[i]==0:
            symm = per_degree[i].is_symmetric()
            print("\t Is symmetric: ", symm)
            if symm:
                print("\t Schur: ", schur(per_degree[i].to_symmetric_function()))
def print_pos_QS(f,n):
    """
    This code  will nicely separate your quasisymmetric function by degree in the QS basis
    and check positivity. 
    Note that n is the homogeneous degree of f, and this is likely to be n*(n-1)/2 where
    n is the number of vertices of the graph g. 
    """
    g = QS(f)
    print(g)
    for i in range(n+1):
        pos_at_deg_i = True
        for A in Compositions(i):
            c = g.coefficient(A)
            if c<0:
                pos_at_deg_i = False
        if pos_at_deg_i == True:
            print("Pos at deg", i)
        else:
            print("Neg at deg", i)
    return 0
def print_pos_YQS(f,n):
    """
    This code  will nicely separate your quasisymmetric function by degree in the YQS basis
    and check positivity. 
    Note that n is the homogeneous degree of f, and this is likely to be n*(n-1)/2 where
    n is the number of vertices of the graph g. 
    """
    g = YQS(f)
    #print("g =\n", g)
    print(g)
    for i in range(n+1):
        #print("looking at", i)
        pos_at_deg_i = True
        for A in Compositions(i):
            #print("\t looking at", A)
            c = g.coefficient(A)
            #print("\t coeff is", c)
            if c<0:
                pos_at_deg_i = False
        if pos_at_deg_i == True:
            print("Pos at deg", i)
        else:
            print("NEG at deg", i)
    return 0






