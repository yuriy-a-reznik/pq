# pq
PQ: An algorithm for quantization of discrete probability distributions.

This is an implementation of a lattice quantization algorithm described in:
 *    Y. A. Reznik, "An Algorithm for Quantization of Discrete Probability Distributions," Proc. Data Compression Conference (DCC'11), pp. 333-343, Snowbird, UT, March 2011.

This algorithm receives simplex-constrained quantities, such as discrete probability distributions, 
and quantizes them to a set of points known as "types" in universal source coding literature. 
Plots below show type locations and corresponding Voronoi partitions in 3 dimensions.

<p align="center">
<img width="498" height="198" alt="tl" src="https://github.com/user-attachments/assets/f2f14a90-cbc0-4068-bdf6-5054a082a624" />
</p>

<p>
  Mathematically, the resulting lattice is equivalent to the so-called
  <em>A<sub>n</sub></em> lattice with an additional spacing constraint.
  This constraint enables a simpler construction and a direct combinatorial
  enumeration of reconstruction points.
</p> 

This repository includes both library (pq.h/pq.c and a testing program pq_test_accuracy.c). 
