#!/usr/bin/env python
# coding: utf-8

load("Unit_Interval_graphs.py")
def rep(L):
    """
    L is a triple of three objects. We return the second object, which is a list,
    but we drop some items in this list. We only keeps the portion of the list
    where the third object is 0. 
    For example if L[1] is 1234 and L[2] is 0101, then we will return 13. 
    """
    L_rep = []
    for i in range(len(L[0])):
        if L[2][i]==0:
            L_rep = L_rep + [L[1][i]]
    return L_rep
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
def perm_chain_to_hyperplane_walk_noting_equiv(chain, g):
    """
    Input chain, a sequence of permutations associated to a chain in the weak order. 
    Also input a filled graph g. 
    Output the associated sequence of hyperplanes, i.e. the hyper plane walk
    as well as a sequence that notes which hyperplanes are "disappeared". 
    The output e has a 1 at the index for those hyperplanes that are "disappeared".
    """
    m = len(chain)
    h = []
    e = []
    for i in range(m-1):
        a = chain[i]
        b = chain[i+1]
        index = 0
        while(index<len(a) and a[index]==b[index]):
            index +=1
        s = a[index]
        t = a[index+1]
        h += [[s,t]]
        if is_cut_vertex_set(g,s,t,a[index+2:]):
            e += [1]
        else:
            e += [0]
    return [h,e]
            

def hyperplanewalks_noting_equivalence(g):
    """
    Given a filled graph g, will produce all hyperplane walks paired 
    with a binary list that notes which hyperplanes are "disappeared".
    """
    n = g.order()
    C = chains_in_weak_order(n)
    H = []
    for chain in C:
        A = perm_chain_to_hyperplane_walk_noting_equiv(chain, g)
        H += [A]
    return H



def are_related_hyperplanewalks(H1,H2):
    """
    Given two hyperplane walks of the weak order we determine if they are
    adjacent in the weak order. i.e. if they differ by exacly two hyperplanes. 
    """
    c = 0
    for i in range(len(H1)):
        if not H1[i]==H2[i]:
            c += 1
            if c>2:
                return False
    return True
def two_classes_are_related(c1,c2):
    """
    Given two classes of hyperplane walks, we determine if there is a pair 
    of hyper planes, H1 in c1 and H2 in c2 that are related. 
    """
    for i in range(len(c1[1])):
        for j in range(len(c2[1])):
            H1 = c1[1][i]
            H2 = c2[1][j]
            if are_related_hyperplanewalks(H1,H2):
                return True
    return False
def walk_rep(H):
    """
    We input H =[walk, equiv] where walk is the hyperplane walk
    and equiv is a note on which hyperplanes are "disappeared".
    If equiv has a 1 at index i, then we drop the associated hyperplane in walk. 
    If equiv has a 0 at index i, then we keep the associated hyperplane in walk. 
    We return the first number of each hyperplane in walk that we keep. 
    """
    r = []
    walk = H[0]
    #print("\t \t walk:", walk)
    equiv = H[1]
    #print("\t \t equiv:", equiv)
    for i in range(len(walk)):
        if equiv[i]==0:
            r += [walk[i][0]]
    return r

def gather_equiv_classes(g):
    """
    Given a filled graph g, we return all equivalence classes. 
    Each class will be given as [rep, list of hyperplane walks]
    Where rep is the first number of each hyperplane in a hyperplane walk that is 
    not "disappeared". This is constant within an equivalence class. 
    We also return the collection of full hyperplane walks in the weak order in 
    the equivalence class. 
    """
    classes = []
    all_walks = hyperplanewalks_noting_equivalence(g)
    #print("all walks:", all_walks)
    for H in all_walks:
        #print("\t H:", H)
        rep = walk_rep(H)
        #print("\t rep:", rep)
        added = False
        for c in classes:
            #print("\t \t c:", c)
            if rep == c[0]:
                c[1] += [H[0]]
                added = True
        if added == False:
            classes += [[rep, [H[0]]]]
        #print("\t classes so far:", classes)
    return classes

def graph_of_chains_of_maximaltubing_poset(g):
    """
    Given a filled graph g, we will form a graph determined 
    where the vertices are the equivalence classes of hyperplane walks 
    and the edges will be when two classes are related.
    """
    classes = gather_equiv_classes(g)
    #print("classes:", classes)
    m = len(classes)
    G = Graph(m)
    for i in range(m):
        for j in range(i+1,m):
            if two_classes_are_related(classes[i],classes[j]):
                G.add_edge(i,j)
    print("What vertices represent as hyperplane walk equivalence classes:")
    for i in range(m):
        print("Vertex ", i, ":", classes[i][1])
    return G
def graph_of_maximalchains_of_maximaltubing_poset(g):
    """
    Given a filled graph g, we will form a graph determined 
    where the vertices are the equivalence classes of hyperplane walks 
    associated to maximal length chains in the tubing poset 
    (i.e. equivalence classes where no hyperplanes are "disappeared")
    and the edges will be when two classes are related.
    """
    classes = gather_equiv_classes(g)
    n = g.order()
    max_chain_length = n*(n-1)/2
    vertices = []
    #print("classes:", classes)
    m = len(classes)
    G = Graph(0)
    for i in range(m):
        if len(classes[i][0])==max_chain_length:
            G.add_vertex(i)
            vertices += [i]
        for j in range(i+1,m):
            if  (len(classes[i][0])==max_chain_length and len(classes[j][0])==max_chain_length):
                if two_classes_are_related(classes[i],classes[j]):
                    G.add_edge(i,j)
    print("What vertices represent as hyperplane walk equivalence classes:")
    for i in range(m):
        if i in vertices:
            print("Vertex ", i, ":", classes[i][1])
    return G
    



def print_info(n):
    """
    We go through all connected filled graphs on n vertices. 
    We print off various information about the graphs formed
    from equivalence classes of hyperplane walks. 
    """
    for g in all_unit_interval_graphs_connected(n):
        print("Origional unit interval graph:")
        g.show()
        print("Graphs of chains of the tubing poset:")
        g1 = graph_of_chains_of_maximaltubing_poset(g)
        diam1 = g1.diameter()
        print("\t Diameter = ", diam1)
        g1.show()
        print("Graphs of LONGEST chains of the tubing poset:")
        g2 = graph_of_maximalchains_of_maximaltubing_poset(g)
        diam2 = g2.diameter()
        print("\t Diameter = ", diam2)
        g2.show()
        print("___________________________________________")






