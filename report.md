# Gorilla or Sea Cucumber

## Results
Our implementation produces the expected results on all pairs of
species.

The closest species to Human is the Gorilla , with the following
optimal alignment

Human--Gorilla: 777
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFKLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH

the most distant species is the Sea-Cucumber with the following optimal alignment:

Human--Sea-Cucumber: 80
M-V--H--LTPEEKSAVTALWGK-V-NVDEVGGEALGRLLVVY-PWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLD-N-LKGTFATLSELHCDKLH-VDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANA-LAHKYH
LAIQAQGDLTLAQKKIVRKTWHQLMRNKTSFVTDVFIRIF-AYDPSAQNKFPQMAGMSA-SQLRSSRQMQAHAIRVSSIMSEYVEELDSDILPELLATLARTH-D-LNKVGADHYNLFAKVLMEALQAELGSDFNEKTRDAWAKAFS-VVQAVLLVKHG

------------------------------


## Implementation details

There are two approaches to the problem, an imperative solution (assambleImp) and
a functional solution (solveFun) wich is using a state Monad the only difference between the two
is that the former is using mutable variables whereas the second one is usnig immutable variables.

The imperative solution use a bottom up approach while creating the table.
It runs through $m*n$ calculations. Each calculation is done in constant time so the total runing time is $O(m*m)$

The recursive implementation that compute the result  from the top to the bottom. 
Our implementation uses $O(nm)$ time since the lookup table has $n*m$ cells and we compute each cell only once. 

We chose a recursive/iterative implementation. For two sequences
of length n and m, respectively, our implementation uses O (( n 3 +
log 2 m ) cos n ) time and O ( 1 ) space
