# Tensor Tracking: Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors

We propose two new provable algorithms for tracking online low-rank approximations of higher-order streaming tensors in the presence of missing data. The first algorithm, dubbed adaptive CP decomposition (ACP), minimizes an exponentially weighted recursive least-squares cost function to obtain the tensor factors in an efficient way, thanks to the alternative minimization framework and the randomized
sketching technique. Under the Tucker model, the second algorithm called adaptive Tucker decomposition (ATD), which is more flexible than the first one, first tracks the underlying low-dimensional subspaces covering
the tensor factors, and then estimates the core tensor using a stochastic approximation. 

Both algorithms are fast, and require a low computational complexity and memory storage.


## DEMO

+ Run "demo_ACP_xyz.m" and "demo_ATD_xyz.m" for synthetic experiments.
+ Run "demo_real_video_tracking_completion.m" for the online tensor completion. Video datasets can be downloaded from Releases. 


## State-of-the-art algorithms for comparison
+ TeCPSGD: Mardani Mortezaet et al. “[*Subspace learning and imputation for streaming big data matrices and tensors*](https://ieeexplore.ieee.org/document/7072498)”. **IEEE Trans. Signal Process.** 63, 10 (2015): 2663-2677.
+ OLSTEC: Kasai Hiroyuki. “[*Fast online low-rank tensor subspace tracking by CP decomposition using recursive least squares from incomplete observations*](https://www.sciencedirect.com/science/article/abs/pii/S0925231218313584?via%3Dihub)”. **Neurocomput.** 347 (2019): 177-190.
+ CP-PETRELS:  Minh-Chinh Truong, et al. “[*Adaptive PARAFAC decomposition for third-order tensor completion*](https://ieeexplore.ieee.org/abstract/document/7562652)”. **IEEE Int. Conf. Commun. Electr. (IEEE ICCE),** 2016.
+ WTucker: Filipović, M., and Jukić, A. “[*Tucker factorization with missing data with application to low-𝑛-rank tensor completion*](https://link.springer.com/article/10.1007/s11045-013-0269-9)”. **Multidim. Syst. Signal Process.** 26, 3 (2015), 677–692.
+ iHOOI and ALSaS: Xu, Y. “[*Fast algorithms for higher-order singular value decomposition from incomplete data*](https://global-sci.org/intro/article_detail/jcm/10023.html)”. **J. Comput. Math.** 35, 4 (2017), 395–420.

## Some Results

+ Streaming CP Decomposition

<img src="https://user-images.githubusercontent.com/26319211/167562867-207d050a-7819-462a-8837-73e417ce0bad.png" width="700" height='500'>

+ Streaming Tucker Decomposition

<img src="https://user-images.githubusercontent.com/26319211/167563914-37454381-d12a-4cfb-9fc3-f8773342634b.png" width="700" height='500'>


+ Video Completion

<img src="https://user-images.githubusercontent.com/26319211/167564258-80b54ed2-23ac-43c0-8b9b-a62db1cc7369.png" width="700" height='500'>

 
## References

This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1] **L.T. Thanh**, K. Abed-Meraim, N. L. Trung and A. Hafiane. “[*Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors*](https://drive.google.com/fi)”. **Elsevier Patterns**, 2023, [CellPress](https://www.cell.com/patterns/fulltext/S2666-3899(23)00104-6), [Techrxiv](https://www.techrxiv.org/articles/preprint/Tracking_Online_Low-Rank_Approximations_of_Higher-Order_Incomplete_Streaming_Tensors/19704034)], [[PDF](https://thanhtbt.github.io/files/2023_Patterns_Tensor_Tracking_Draw.pdf)]. 




