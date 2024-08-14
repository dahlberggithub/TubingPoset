#Code for all unit interval graphs

def lift(L):
    """
    We input L, a list of length 2 lists. 
    We will add one to the second coordinate of all the length 2 lists. 
    """
    L1 = []
    for k in range(len(L)):
        L1 = L1 + [[L[k][0],L[k][1]+1]]
    return L1

def inc(L,k):
    """
    L is a list of length 2 list. We will increase all corrdinates by k. 
    """
    L1 = []
    for i in range(len(L)):
        L1 = L1 + [[L[i][0]+k,L[i][1]+k]]
    return L1


def intervals(n):
    """
    Will return list of intervals associated to a unit interval graph g on n vertices.
    For example if the intervals are [1,2],[2,4], then we have edges 12 for the [1,2]
    and edges 23, 24, 34 for the interval [2,4]. This method uses an induction
    similar one for the Catalan numbers. 
    """
    if n == 1:
        #print "begining!"
        return [[[1,1]]]
    I = []
    I_old = intervals(n-1)
    #print "I_old:", I_old
    for k in range(len(I_old)):
        L = I_old[k]
        #print "adding ", [[1,1]]+inc(L,1)
        I = I + [[[1,1]]+inc(L,1)]
    for k in range(1,n-1):
        I_old1 = intervals(k)
        I_old2 = intervals(n-k-1)
        for i in range(len(I_old1)):
            for j in range(len(I_old2)):
                L1 = I_old1[i]
                L2 = I_old2[j]
                #print "from ", L1, "and", L2
                L = lift(L1)+inc(L2,k+1)
                #print "adding", L
                I = I + [L]
    I_old = intervals(n-1)
    for k in range(len(I_old)):
        L = lift(I_old[k])
        I = I + [L]
    return I

def intervalsConnected(n):
    """
    We return all intervals associated to unit interval graphs on n vertices that are connected. 
    """
    I = intervals(n-1)
    J = []
    for i in range(len(I)):
        L = I[i]
        J = J + [lift(L)]
    return J

def unit_interval_graph(L,n):
    """
    Given a list L of intervals associated to a unit interval graph on n vertices, 
    we will return the associated unit interval graph on vertices 0,1,2,...,n-1. 
    """
    g = Graph(0)
    for k in range(len(L)):
        a = L[k][0]+1
        b = L[k][1]+1
        for i in range(a,b):
            for j in range(i+1,b+1):
                g.add_edge(i-1,j-1)
    return g
def all_unit_interval_graphs(n):
    """
    We will return a list of all unit interval graphs on n vertices. 
    The vertices have labels 0,1,...,n-1. 
    """
    G = []
    for L in intervals(n):
        G = G + [unit_interval_graph(L,n)]
    return G
def all_unit_interval_graphs_connected(n):
    """
    We will return a list of all connected unit interval graphs on n vertices. 
    The vertices have labels 0,1,...,n-1. 
    """
    G = []
    for L in intervalsConnected(n):
        G = G + [unit_interval_graph(L,n)]
    return G
print("Unit interval graph code loaded in!")
print("Try all_unit_interval_graphs(n) or all_unit_interval_graphs_connected(n).")