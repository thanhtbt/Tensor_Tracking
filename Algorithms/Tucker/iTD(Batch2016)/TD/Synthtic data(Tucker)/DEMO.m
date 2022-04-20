clear;clc;
M = 10 ;
N = 10;
P = 10;

R1 = 3;
R2 = 4;
R3 = 5;

Sampling = 0.9; 
noise_fac = 1e-3;


%% Form tensor
U1=randn(M,R1);
U2=randn(N,R2);
U3=randn(P,R3);
X = tensor(randn(R1,R2,R3));
Y_real = ttm(X,{U1,U2,U3},1:3);
% Y_real = Y_real/norm(Y_real)*sqrt(M*N*P);
noise = tensor(randn(M,N,P));

Y_n = Y_real + noise_fac * noise; %(Y_real + 0*norm(Y_real)/norm(noise)*noise);

Omiga = tensor(rand(M,N,P)>Sampling);
Y = Y_n.*Omiga;

%% IRTD
Max_iter = 10000;
delta = 0.1;
t_start = tic;
Result3 = IRTD(Y,Omiga,Max_iter,[M,N,P],delta,Y_real);
toc(t_start)

er3 = norm((Result3.Y-Y_real)) / norm(Y_real)
%% LRTI - Very fast
Lambda   = 0.3;
Max_iter = 100;
t_start  = tic;
Result1  = LRTI(Y,Omiga,Max_iter,Lambda);
toc(t_start)
er1 = norm((Result1.Y-Y_real)) / norm(Y_real)

%% HaLRTC
Max_iter = 100;
rou = 1e-5;
tic;
t_start  = tic;
Result2 = HaLRTC(Y,Omiga,rou,Max_iter);
toc(t_start)
er2 = norm((Result2.Y-Y_real)) / norm(Y_real)

%% WTucker
X = double(Y);
ind = find(double(Omiga));
rest = 1*[R1,R2,R3];
initial = randn(M,N,P);
topt_maxiter = 300;
poblano_maxiter = 100;
topt_verbose = true;
t_start = tic;
[Xrec_topt_adapt, ~] = ...
    tucker_opt(size(X), X(ind), ind, rest, ...
    'Xinit', initial, 'maxiter', topt_maxiter, ...
    'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
toc(t_start)

er4 =  norm(tensor(Xrec_topt_adapt)-Y_real)/norm(Y_real)
