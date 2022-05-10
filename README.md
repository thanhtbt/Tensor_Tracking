# Tensor Tracking: Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors

We propose two new provable algorithms for tracking online low-rank approximations of higher-order streaming tensors in the presence of missing data. The first algorithm, dubbed adaptive CP decomposition (ACP), minimizes an exponentially weighted recursive least-squares cost function to obtain the tensor factors in an efficient way, thanks to the alternative minimization framework and the randomized
sketching technique. Under the Tucker model, the second algorithm called adaptive Tucker decomposition (ATD), which is more flexible than the first one, first tracks the underlying low-dimensional subspaces covering
the tensor factors, and then estimates the core tensor using a stochastic approximation. 

Both algorithms are fast, and require a low computational complexity and memory storage.


## DEMO

+ Run "demo_ACP_xyz.m" and "demo_ATD_xyz.m" for synthetic experiments.
+ Run "demo_real_video_tracking_completion.m" for the online tensor completion. Video datasets can be downloaded from Releases. 


## State-of-the-art algorithms for comparison
+ TeCPSGD: Mardani Mortezaet et al. "Subspace learning and imputation for streaming big data matrices and tensors." IEEE Trans. Signal Process. 63.10 (2015): 2663-2677.
+ OLSTEC: Kasai Hiroyuki. "Fast online low-rank tensor subspace tracking by CP decomposition using recursive least squares from incomplete observations." Neurocomputing 347 (2019): 177-190.
+ CP-PETRELS:  Minh-Chinh Truong, et al. "Adaptive PARAFAC decomposition for third-order tensor completion." IEEE Int. Conf. Commun. Electr.(IEEE ICCE), 2016.
+ WTucker: Filipović, M., and Jukić, A. "Tucker factorization with missing data with application to low-𝑛-rank tensor completion". Multidim. Syst. Signal Process. 26, 3 (2015), 677–692.
+ iHOOI and ALSaS: Xu, Y. Fast algorithms for higher-order singular value decomposition from incomplete data. J. Comput. Math. 35, 4 (2017), 395–420

## Some Results

Similated data:

+ Streaming CP Decomposition

<img src="https://user-images.githubusercontent.com/26319211/167562867-207d050a-7819-462a-8837-73e417ce0bad.png" width="700" height='500'>

+ Streaming Tucker Decomposition

<img src="https://user-images.githubusercontent.com/26319211/167563914-37454381-d12a-4cfb-9fc3-f8773342634b.png" width="700" height='500'>


+ Video Completion

<img src="https://user-images.githubusercontent.com/26319211/167564258-80b54ed2-23ac-43c0-8b9b-a62db1cc7369.png" width="700" height='500'>

+ EEG Analysis
+ 
<img src="https://user-images.githubusercontent.com/26319211/167564389-dcd25c01-1210-44f7-8684-48931144589d.JPG" width="700" height='500'>


## References

This code is free and open source for research purposes. If you use this code, please acknowledge the following papers.

[1] L.T. Thanh, K. Abed-Meraim, N. L. Trung and A. Hafiane. “[*A Fast Randomized Adaptive CP Decomposition For Streaming Tensors*](https://drive.google.com/file/d/1DAUTPryASpIoDxUZlRW_jzMSFeOS5EPm/view?usp=sharing)”. **Proceedings of the 46th IEEE ICASSP**, 2021. [[DOI](https://ieeexplore.ieee.org/document/9413554)], [[PDF](https://drive.google.com/file/d/1DAUTPryASpIoDxUZlRW_jzMSFeOS5EPm/view?usp=sharing)].

[2] L.T. Thanh, K. Abed-Meraim, N. L. Trung and A. Hafiane. “[*Tracking Online Low-Rank Approximations of Higher-Order Incomplete Streaming Tensors*](https://drive.google.com/fi)”. **Techrxiv**, 2022, [[DOI](https://www.techrxiv.org/articles/preprint/Tracking_Online_Low-Rank_Approximations_of_Higher-Order_Incomplete_Streaming_Tensors/19704034)], [[PDF](https://drive.google.com/file/d/1BHT80liD97fwKdaugxh5eSs3uEKBRLye/view?usp=sharing)]. 




