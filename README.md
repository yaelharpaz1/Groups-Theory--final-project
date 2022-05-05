# Groups-Theory--final-project
Implementation of the paper “L1 rotation averaging using the Weiszfeld algorithm”. 

# Definition of the problem:
The paper considers the problem of rotation averaging under the L1 norm, with comparison to L2 norm. The problem of rotation averaging includes two topics: 
1. Single rotation averaging: Several estimates are obtained of a single rotation, which are then averaged to give the best estimate.
2. Multiple rotation averaging: Relative rotations R_ij are given, and absolute rotations R_i are computed to satisfy the compatibility constraint R_ij R_i=R_j.

# Applications (computer vision):
1.	Structure from motion (SFM): photogrammetric range imaging technique for estimating 3D  structures from 2D image sequences that may be coupled with local motion signals  (By multiple rotation averaging).
2.	Non overlapping camera calibration (By single rotation averaging).
3.	Computing relative rotations between pairs of images (By single rotation averaging).
