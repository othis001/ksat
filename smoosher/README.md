# smoosher

Program used to compute homotopy types for the following article: http://arxiv.org/abs/0808.0739 .

As input we take a text file containing a collection of same sized subsets of the positive integers. Each separated with a 0 and a -1 to denote the end of that collection.

For example: 1 2 0 3 4 0 5 6 -1

The program will then form a simplicial complex using this data where the vertices are the subsets and we form an N-1 simplex whenever N subsets have a nonempty intersection. We then perform simple homotopy equivalences or "smooshings" in order to reduce the complex.

There are several command line options available. The default behavior is just to print the number of remaining simplices in each dimension.

-up (starts from the bottom instead of the top)

-simp (displays the simpices remaining)

-matrix (saves the resulting boundary matrices)

-mem (switches to program to memory efficient mode at the cost of performance)
