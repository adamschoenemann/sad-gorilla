---
documentclass: scrartcl
title: Gorilla or Sea Cucumber
author:
- Oscar Felipe Toro
- Adam Schønemann
- Andreas Jakobsen
- Henrik Sloth Schade
subtitle: "Group 14, Algorithm Design"
---

## Results
Our implementation produces the expected results on all pairs of
species.

The closest species to Human is the Gorilla , with the following
optimal alignment

    Human--Gorilla: 777
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
    KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKE
    FTPPVQAAYQKVVAGVANALAHKYH
    
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
    KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFKLLGNVLVCVLAHHFGKE
    FTPPVQAAYQKVVAGVANALAHKYH

the most distant species is the Sea-Cucumber with the following optimal alignment:

    Human--Sea-Cucumber: 80
    M-V--H--LTPEEKSAVTALWGK-V-NVDEVGGEALGRLLVVY-PWTQRFFESFGDLSTPD
    AVMGNPKVKAHGKKVLGAFSDGLAHLD-N-LKGTFATLSELHCDKLH-VDPENFRLLGNVL
    VCVLAHHFGKEFTPPVQAAYQKVVAGVANA-LAHKYH
    
    LAIQAQGDLTLAQKKIVRKTWHQLMRNKTSFVTDVFIRIF-AYDPSAQNKFPQMAGMSA-S
    QLRSSRQMQAHAIRVSSIMSEYVEELDSDILPELLATLARTH-D-LNKVGADHYNLFAKVL
    MEALQAELGSDFNEKTRDAWAKAFS-VVQAVLLVKHG

## Implementation details
There are two approaches to the problem, an imperative solution (`assambleImp`) and
a functional solution (`assambleFun`) which uses a `State` Monad.

- Functional solution
    - No mutable state
    - Recursion instead of iteration
    - Top-down computation
    - Uses a `Map[(Int,Int),Int]` for memoization
    - Uses a State monad
      using `unfoldRight`
- Imperative solution
    - Lots of mutable state
    - Uses nested for loops for iteration
    - Computes the result bottom-up
    - Uses a 2D array for memoization

Both solutions runs through $m ⋅ n$ calcuations, and since each calculation is
done in constant time, the total running time is $O(m ⋅ n)$.