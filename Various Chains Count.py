#!/usr/bin/env python
# coding: utf-8



load("Unit_Interval_graphs.py")
QSym = QuasiSymmetricFunctions(QQ)
YQS = QSym.YQS()
Fund = QSym.F()
QS = QSym.QS()
M = QSym.M()

def chains_in_weak_order_to_k_inversions(n,k):
    """
    We build all chains in the week order starting at the bottom element, the identity permutation. 
    We build up the chains as sequences of permutations. 
    We return all chains length k starting from the identity. 
    """
    if k == 0:
        return [[[i+1 for i in range(n)]]]
    Old = chains_in_weak_order_to_k_inversions(n,k-1)
    R = []
    for chain in Old:
        top = chain[-1]
        for i in range(n-1):
            if top[i]<top[i+1]:
                new_top = top.copy()
                new_chain = chain.copy()
                new_top[i], new_top[i+1] = top[i+1], top[i]
                R = R + [new_chain + [new_top]]
    return R
def chains_in_weak_order(n):
    """
    We return all maximal chains of the weak order. 
    The chains are presented as sequences of permutations on 1,2,...,n
    """
    m = n*(n-1)/2
    return chains_in_weak_order_to_k_inversions(n,m)

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

def chain_information(g):
    """
    We will take all chains (given as a list of permutations) in the week order and 
    (1) construct the associated reduced word
    (2) determine the hyperplane sequence (we will just the lower coordinate of the hyperplane)
    (3) using graph g we will note which hyperplanes are considered "deleted"
        (ie. when switching jk in the permutation, there exists a cut set to the right)
    """
    n = g.order()#number of vertices
    all_seq_of_perms = chains_in_weak_order(n)
    R = []
    for chain in all_seq_of_perms:
        word = []
        hyper = []
        delete = []
        for i in range(len(chain)-1):
            before = chain[i]
            after = chain[i+1]
            index = -1 #will look for the s_i that changes before to after
            what_moved_right = -1 #will look for what number got switched right
            is_deleted = -1 #will note if there is a cut vertex to the right of the switch (numbers between j and k)
            what_moved_left = -1
            m = 0
            while index == -1:
                if not before[m]==after[m]:
                    index = m
                    what_moved_right = before[m]
                    what_moved_left = before[m+1]
                m = m + 1
            L = [] #list everything to the right of index that is between what_moved_right and what_moved_left
            for j in range(index,n):
                if before[j]<what_moved_left and before[j]>what_moved_right:
                    L = L + [before[j]]
            if is_cut_vertex_set(g,what_moved_right,what_moved_left,L):
                is_deleted = 1
            else:
                is_deleted = 0
            word = word + [index+1]
            hyper = hyper + [[what_moved_right,what_moved_left]]
            delete = delete + [is_deleted]
        R = R + [[word,hyper,delete]]
    return R
def to_equiv_classes(g):
    """
    Given the chain information output we gather together all reduced words 
    in each equivalence class. 
    Two reduced words are in the same equivalence class if the associated hyperplane sequence
    is equal once you remove the hyperplanes that are considered "deleted".
    We return a list of equivalence classes where the first portion
    is the hyperplane walk with the deleted hyperplanes and the second coordinate
    is a list of all reduced words in that class. 
    """
    all_chains = chain_information(g)
    #print("all_chains:", all_chains)
    R = []
    for chain in all_chains:
        #print("chain:", chain)
        word = chain[0]
        hyper = chain[1]
        delete = chain[2]
        hyper_shortened = [hyper[i][0] for i in range(len(hyper)) if delete[i]==0]
        class_found = False
        for equiv_class in R:
            if equiv_class[0]==hyper_shortened:
                class_found = True
                equiv_class[1] += [chain]
        if class_found == False:
            R = R + [[hyper_shortened,[chain]]]
        #print("R so far:", R)
    return R
def print_shortest_chains(g):
    """
    Given a graph g we will return the shortest maximal chains. 
    We also print out various information.
    """
    E = to_equiv_classes(g)
    n = g.order()
    m = n*(n-1)/2
    p = n*(n-1)/2
    for c in E:
        a = len(c[0])
        if a<m:
            m = a
    print("Shortest chain length:", m)
    count = 0
    print("Class reps:")
    for c in E:
        if len(c[0])==m:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of shortest chains:", count)
    return 0

def print_shortestandlongest_chains(g):
    """
    Given a graph g we will return the shortest maximal chains and longest maximal chains. 
    We also print out various information.
    """
    E = to_equiv_classes(g)
    n = g.order()
    m = n*(n-1)/2#shortest
    M = 0#longest
    p = n*(n-1)/2
    for c in E:
        a = len(c[0])
        if a<m:
            m = a
        if a>M:
            M = a
    print("Shortest chain length:", m)
    count = 0
    print("Class reps of shortest chains:")
    for c in E:
        if len(c[0])==m:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of shortest chains:", count)
    print("longest chain length:", M)
    count = 0
    print("Class reps of longest chains:")
    for c in E:
        if len(c[0])==M:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of longest chains:", count)
    return 0
def print_all_chains(g):
    """
    Given a graph g we will return the all maximal chains. 
    We also print out various information.
    """
    E = to_equiv_classes(g)
    n = g.order()
    m = n*(n-1)/2#shortest
    M = 0#longest
    p = n*(n-1)/2
    print("Total number chains: ", len(E))
    for c in E:
        a = len(c[0])
        if a<m:
            m = a
        if a>M:
            M = a
    count = 0
    print("Class reps of all chains:")
    for c in E:
        if len(c[0])>-1:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of chains:", count)
    print("Shortest chain length:", m)
    count = 0
    print("Class reps of shortest chains:")
    for c in E:
        if len(c[0])==m:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of shortest chains:", count)
    print("longest chain length:", M)
    count = 0
    print("Class reps of longest chains:")
    for c in E:
        if len(c[0])==M:
            print("\t",[c[1][0][1][i] for i in range(p) if c[1][0][2][i]==0])
            count += 1
    print("Number of longest chains:", count)
    return 0





def print_shortest_chain_information_for_unitintervalgraphs(n):
    """
    Looking over all connected unit interval graphs we will print out
    information about the shortest maximal chains. 
    """
    for g in all_unit_interval_graphs_connected(n):
        g.show(figsize=2)
        print_shortest_chains(g)
        print("_________________________")
    return 0
def print_longestshortest_chain_information_for_unitintervalgraphs(n):
    """
    Looking over all connected unit interval graphs we will print out
    information about the shortest and longest maximal chains. 
    """
    for g in all_unit_interval_graphs_connected(n):
        g.show(figsize=2)
        print_shortestandlongest_chains
        print("_________________________")
    return 0
def print_all_chain_information_for_unitintervalgraphs(n):
    """
    Looking over all connected unit interval graphs we will print out
    information about the shortest maximal chains. 
    """
    for g in all_unit_interval_graphs_connected(n):
        g.show(figsize=2)
        print_all_chains(g)
        print("_________________________")
    return 0





