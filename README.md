# TubingPoset

This is a collection of code that creates and does calculations for the Poset of Maximal Tubings for a filled graph G. You can find information about their definition in the following paper by Emily Barnard and Thomas McConville, [Link](https://arxiv.org/pdf/1808.05670). 

The Poset of Maximal Tubings for a filled graph G is a lattice quotient of the Weak Order on permutations of {1,2,...,n}. Esentially we start with the Weak Order on permutations. This is the collection of all permutations with a partial ordering. We say that a permutation sigma is larger than pi if we can find a sequence of adjacent switches of ascents (i.e. when pi(i)<pi(i+1)) that transforms pi into sigma. For example 2134 is less than 2413 becuase we can switch the ascent 34 in 2134 to get 2143 and then we can switch the ascent 14 to get 2413. 

We will identify two permutations if they have the same G-tree. The G-tree of a permutation is defined inductively. Given a permutation pi of {1,2,...,n} we start with the root of our tree, pi(n). We then remove pi(n) from G and potentially get several connected commponents G_1, G_2,.... Suppose when we restrict pi to the vertices in G_1, G_2,... we get w_1,w_2,.... Then w_1 will have an associated G_1-tree, w_1 will have an associated G_2-tree and so on. We connect pi(n) to the roots of the G_1-tree, G_2-tree,.... to get our G-tree. 

We then identify permutations in the Weak Order that have the same G-tree. The resulting object is still a poset. The code in this repository studies this poset.  
