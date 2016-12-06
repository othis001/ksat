# Boolean Formulae, Hypergraphs and Combinatorial Topology

A collection of C++ programs used for the following article: http://arxiv.org/abs/0808.0739 .
The descriptions are given below.

## findEuler

Program used to compute Euler characteristics for the following article: http://arxiv.org/abs/0808.0739 .

As input we take a text file containing a collection of same sized subsets of the positive integers. Each separated with a 0 and a -1 to denote the end of that collection.

For example: 1 2 0 3 4 0 5 6 -1

The program will then form a simplicial complex using this data where the vertices are the subsets and we form an N-1 simplex whenever N subsets have a nonempty intersection. By counting the number of simplices in each dimension and then taking an alternating sum, we determine the Euler characteristic.

## findHomology

Program used to compute rational homology groups for the following article: http://arxiv.org/abs/0808.0739 .

As input we take a text file containing a collection of same sized subsets of the positive integers. Each separated with a 0 and a -1 to denote the end of that collection.

For example: 1 2 0 3 4 0 5 6 -1

The program will then form a simplicial complex using this data where the vertices are the subsets and we form an N-1 simplex whenever N subsets have a nonempty intersection. By using the Linbox library to compute ranks of boundary matrices, we determine the rational homology groups.

## smoosher

Program used to compute homotopy types for the following article: http://arxiv.org/abs/0808.0739 .

As input we take a text file containing a collection of same sized subsets of the positive integers. Each separated with a 0 and a -1 to denote the end of that collection.

For example: 1 2 0 3 4 0 5 6 -1

The program will then form a simplicial complex using this data where the vertices are the subsets and we form an N-1 simplex whenever N subsets have a nonempty intersection. We then perform simple homotopy equivalences or "smooshings" in order to reduce the complex.

There are several command line options available. The default behavior is just to print the number of remaining simplices in each dimension.

-up (starts from the bottom instead of the top)

-simp (displays the simpices remaining)

-matrix (saves the resulting boundary matrices)

-mem (switches to program to memory efficient mode at the cost of performance)


