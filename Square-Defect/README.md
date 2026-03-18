# Square Defect Computations

This directory contains Maple functions for computing the square defect of an S-Hadamard matrix, and the values we computed for a select number of non-isolated S-Hadamard BH matrices.
In `defect_functions.mpl`, we provide functions `ComputeDefect` and `ComputeDefectSq` that can be used to compute the defect and square defect.
These functions also return a basis for the kernel of the corresponding linear approximation.
See `compute_defect.mpl` for more detailed information regarding our usage of these functions to compute defect values.

Script `compute_defect.mpl` stores the output objects from the square defect calculations in `output`, see this directory for more details.
The script also stores any S-Hadamard full SLA-parametrizations found in this directory as well.
The expected output of `compute_defect.mpl` is as follows, where for each set of S-Hadamard $BH(n,q)$ matrices we computed the defect values for, we print the index of the matrix in our listing, the square defect value, and if we found a full SLA-parametrization:

```txt
> read "compute_defect.mpl";
Square Defect Calculations for BH(12, 3) Matrices ...
Index      Square Defect   Full SLA  
1          12              false     
2          12              false     
Square Defect Calculations for BH(12, 6) Matrices ...
Index      Square Defect   Full SLA  
1          12              false     
2          3               true      
3          3               true      
4          12              false     
5          3               true      
6          3               true      
Square Defect Calculations for BH(12, 12) Matrices ...
Index      Square Defect   Full SLA  
1          12              false     
2          3               true      
3          3               true      
4          3               true      
5          3               true      
6          3               true      
7          3               true      
8          3               true      
9          12              false     
10         3               true      
11         3               true      
12         3               true      
13         3               true      
14         3               true      
15         3               true      
16         3               true      
17         0               true      
Square Defect Calculations for BH(16, 4) Matrices ...
Index      Square Defect   Full SLA  
1          4               true      
2          6               false     
3          4               true      
4          5               false     
5          4               true      
6          4               true      
7          20              false     
8          9               false     
9          4               true      
10         4               true      
11         5               false     
12         1               true      
13         3               false     
14         3               false     
15         3               false     
Square Defect Calculations for BH(18, 3) Matrices ...
Index      Square Defect   Full SLA  
1          20              false     
2          10              true      
3          10              true      
4          20              false     
5          10              true      
6          10              true      
7          16              false     
8          20              false     
9          16              false     
10         28              false     
11         12              false     
12         12              false     
13         10              true      
14         28              false     
15         10              true      
16         36              false     
17         16              false     
18         20              false     
19         20              false     
20         40              false     
21         4               true      
22         4               true      
23         4               true      
24         4               true      
25         4               true      
26         4               true      
27         28              false     
28         28              false     
29         2               true      
30         16              false     
31         2               true      
32         16              false     
33         40              false     
34         10              false     
35         10              false     
36         10              false     
37         10              false     
38         4               true      
39         20              false     
40         4               true      
41         24              false     
42         8               false     
43         8               false     
44         12              false     
45         8               false     
46         28              false     
47         8               false     
48         12              false     
49         36              false     
50         2               true      
51         16              false     
52         4               true      
53         12              false     
54         10              false     
55         10              false     
56         10              true      
57         10              true      
58         4               true      
59         12              false     
60         10              true      
61         8               false     
62         16              false     
63         8               false     
64         2               true      
65         12              false     
66         10              false     
67         4               true      
68         16              false     
69         10              false     
70         10              true      
71         10              true      
72         4               true      
73         10              true      
74         12              false     
75         8               false     
76         16              false     
77         8               false     
78         4               true      
79         4               true      
80         20              false     
81         24              false     
82         20              false     
83         20              false     
84         4               true      
85         4               true 
```

The script also prints the defect values for $BH(n,q)$ matrices where $q$ is not equal to 3:

```txt
Defect Calculations for BH(12, 6) Matrices ...
Index      Defect         
1          12             
2          10             
3          13             
4          12             
5          10             
6          13             
Defect Calculations for BH(12, 12) Matrices ...
Index      Defect         
1          12             
2          10             
3          10             
4          9              
5          11             
6          10             
7          9              
8          13             
9          12             
10         10             
11         10             
12         9              
13         11             
14         10             
15         9              
16         13             
17         17             
Defect Calculations for BH(16, 4) Matrices ...
Index      Defect         
1          19             
2          31             
3          19             
4          35             
5          17             
6          17             
7          35             
8          21             
9          15             
10         19             
11         21             
12         15             
13         19             
14         27             
15         21  
```
