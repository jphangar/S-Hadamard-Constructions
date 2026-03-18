# PSLQ-Based Approach to Parametrizing Matrices

This directory contains the implementation and usage of Algorithm 1 (the subroutine for PSLQ-based parametrization) described in Section 3.2.1 of the thesis.
Algorithm 1 is implemented in Maple, given in `PSLQ_Param_Subroutine.mpl`.
Note that the implementation differs from the algorithm description slightly as we prefer to keep the parametrized matrix fully symbolic; we assign labels to roots of unity that appear in the parametrization instead of exact values (which is what the algorithm description does).

We provide an implementation of the naive approach (`PSLQ_Try_Param_Naive.mpl`) and the randomized approach (`PSLQ_Try_Param_Rand.mpl`) for searching for parametrizations using the subroutine.
We choose to use the randomized approach, more details about the specific inputs given can be seen in `PSLQ_Main.mpl`.

The 79 order 12 S-Hadamard matrices found by our optimizer are given in `order-12-optimizer-solns.m`.
The file can be read in Maple to get a list called `soln_mats` of length 79 containing numeric order 12 S-Hadamard matrices.

To reproduce results, run `PSLQ_Main.mpl`.
Note that we have hardcoded this script to only try 3 parameters for simplicity.
This script prints the following output and saves the resulting parametrizations as a Maple list called `validMatrices` to the file `pslq_valid_param_mats.m`:

```txt
> read "PSLQ_Main.mpl";
Index      Parameters      Roots      S-Hadamard?
1          3               [[3, w]]   true      
2          3               [[3, w]]   true      
3          3               [[3, w]]   true      
4          3               [[3, w]]   true      
5          3               [[3, w]]   true      
6          3               [[3, w]]   true      
7          3               [[3, w]]   true      
8          3               [[3, w]]   true      
9          3               [[3, w]]   true      
10         3               [[3, w]]   true      
11         3               [[3, w]]   true      
12         3               [[3, w]]   true      
13         3               [[3, w]]   true      
14         3               [[3, w]]   true      
15         3               [[3, w]]   true      
16         3               [[3, w]]   true      
17         3               [[3, w]]   true      
18         3               [[3, w]]   true      
19         3               [[3, w]]   true      
20         3               [[3, w]]   true      
21         3               [[3, w]]   true      
22         3               [[3, w]]   true      
23         3               [[3, w]]   true      
24         3               [[3, w]]   true      
25         3               [[3, w]]   true      
26         3               [[3, w]]   true      
27         3               [[3, w]]   true      
28         3               [[3, w]]   true      
29         3               [[3, w]]   true      
30         3               [[3, w]]   true      
31         3               [[3, w]]   true      
32         3               [[3, w]]   true      
33         3               [[3, w]]   true      
34         3               [[3, w]]   true      
35         3               [[3, w]]   true      
36         3               [[3, w]]   true      
37         3               [[3, w]]   true      
38         3               [[3, w]]   true      
39         3               [[3, w]]   true      
40         3               [[3, w]]   true      
41         3               [[3, w]]   true      
42         3               [[3, w]]   true      
43         3               [[3, w]]   true      
44         3               [[3, w]]   true      
45         3               [[3, w]]   true      
46         3               [[3, w]]   true      
47         3               [[3, w]]   true      
48         3               [[3, w]]   true      
49         3               [[3, w]]   true      
50         3               [[3, w]]   true      
51         3               [[3, w]]   true      
52         3               [[3, w]]   true      
53         3               [[3, w]]   true      
54         3               [[3, w]]   true      
55         3               [[3, w]]   true      
56         3               [[3, w]]   true      
57         3               [[3, w]]   true      
58         3               [[3, w]]   true      
59         3               [[3, w]]   true      
60         3               [[3, w]]   true      
61         3               [[3, w]]   true      
62         3               [[3, w]]   true      
63         3               [[3, w]]   true      
64         3               [[3, w]]   true      
65         3               [[3, w]]   true      
66         3               [[3, w]]   true      
67         3               [[3, w]]   true      
68         3               [[3, w]]   true      
69         3               [[3, w]]   true      
70         3               [[3, w]]   true      
71         3               [[3, w]]   true      
72         3               [[3, w]]   true      
73         3               [[3, w]]   true      
74         3               [[3, w]]   true      
75         3               [[3, w]]   true      
76         3               [[3, w]]   true      
77         3               [[3, w]]   true      
78         3               [[3, w]]   true      
79         3               [[3, w]]   true   
```
