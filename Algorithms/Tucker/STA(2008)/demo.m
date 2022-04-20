%% Example Title
% Summary of example objective
clear; clc;


addpath('D:\GoogleDrive\2. MEGA\My Works\My Research\Tensor\Codes\TensorToolbox\Tensor_Toolbox_Kolda2019')

tensor_dim  = [10 15 16 500];
tucker_rank = [5 3 5 3];

fac_noise   = 1e-2;
sampling    = 0.8;
epsilon     = 0*eye(1,tensor_dim(end));
% epsilon(600)= 1;

[Xtrue_cell,Factor_True] = tucker_generate_tensor(tensor_dim,tucker_rank,0,epsilon);
OPTS.Factor_True = Factor_True;
OPTS.Slide_True  = Xtrue_cell;


X_cell = cell(tensor_dim(end),1);
for t  = 1 : tensor_dim(end)
    X_cell{t,1} = Xtrue_cell{t,1} + fac_noise*randn(tensor_dim(1:end-1));
    X_cell{t,1} = X_cell{t,1};
end

[T,S] = STA(X_cell{1,1},tucker_rank(1:end-1));
alpha = 0.1;
N = 4;
Factor_Es = cell(1,N);

t_start = tic;
[Factor_Es,PER1] = DTA_Tracking(X_cell,tucker_rank,OPTS);
toc(t_start);;
t_start = tic;
[Factor_Es,PER2] = STA_Tracking(X_cell,tucker_rank,OPTS);
toc(t_start);;

figure;
for ii  = 1 :  length(tucker_rank)
    subplot(length(tucker_rank),1,ii)
    semilogy(PER1(ii,:));
    hold on;
    semilogy(PER2(ii,:));
    axis([1 tensor_dim(end) 1e-6 1e-0])
end

