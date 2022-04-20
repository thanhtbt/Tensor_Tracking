# TT: Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors

We propose two new provable algorithms for tracking online low-rank approximations of higher-order streaming tensors in the presence of missing data. The first algorithm, dubbed adaptive CP
decomposition (ACP), minimizes an exponentially weighted recursive least-squares cost function to obtain the tensor factors in an efficient way, thanks to the alternative minimization framework and the randomized
sketching technique. Under the Tucker model, the second algorithm called adaptive Tucker decomposition (ATD), which is more flexible than the first one, first tracks the underlying low-dimensional subspaces covering
the tensor factors, and then estimates the core tensor using a stochastic approximation. 

Both algorithms are fast, and require a low computational complexity and memory storage.


## DEMO

+ Run files "demo_ACP_xyz.m" and "demo_ATD_xyz" for synthetic data.
+ Run "demo_real_video_tracking_completion.m" for a practical applications of online tensor completion. Video datasets can be downloaded from Releases 


## State-of-the-art algorithms for comparison
+ TeCPSGD: Mardani Mortezaet et al. "Subspace learning and imputation for streaming big data matrices and tensors." IEEE Trans. Signal Process. 63.10 (2015): 2663-2677.
+ OLSTEC: Kasai Hiroyuki. "Fast online low-rank tensor subspace tracking by CP decomposition using recursive least squares from incomplete observations." Neurocomputing 347 (2019): 177-190.
+ CP-PETRELS:  Minh-Chinh Truong, et al. "Adaptive PARAFAC decomposition for third-order tensor completion." IEEE Int. Conf. Commun. Electr.(IEEE ICCE), 2016.
+ WTucker: Filipović, M., and Jukić, A. "Tucker factorization with missing data with application to low-𝑛-rank tensor completion". Multidim. Syst. Signal Process. 26, 3 (2015), 677–692.
+ iHOOI and ALSaS: Xu, Y. Fast algorithms for higher-order singular value decomposition from incomplete data. J. Comput. Math. 35, 4 (2017), 395–420




## References

This code is free and open source for research purposes. If you use this code, please acknowledge the following papers.

[1] L.T. Thanh, K. Abed-Meraim, N. L. Trung and A. Hafiane. “[*A Fast Randomized Adaptive CP Decomposition For Streaming Tensors*](https://drive.google.comg)”. **Proceedings of the 46th IEEE ICASSP**. [[DOI](https://ieeexplore.ieee.org/document/9413554)].

[2] L.T. Thanh, K. Abed-Meraim, N. L. Trung and A. Hafiane. “[*Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors*](https://drive.google.com/fi)”. **ACM Journal** (submitted). 
