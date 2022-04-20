clear;clc;
M = 10 ;
N = 10;
P = 10;

T = 500;

R1 = 3;
R2 = 3;
R3 = 3;
R4 = 3;


Sampling  =  0.9; 
noise_fac = 1e-3;


%% Form tensor
U1=randn(M,R1);
U2=randn(N,R2);
U3=randn(P,R3);
U4=randn(T,R4);

X = tensor(randn(R1,R2,R3,R4));
Y_real = ttm(X,{U1,U2,U3,U4},1:4);
noise = tensor(randn(M,N,P,T));

Y_n = Y_real + noise_fac * noise; %(Y_real + 0*norm(Y_real)/norm(noise)*noise);

Omiga = tensor(rand(M,N,P,T)>Sampling);
Y = Y_n.*Omiga;



%% WTucker
X = double(Y);
ind = find(double(Omiga));
rest = 1*[R1,R2,R3,R4];
initial = randn(M,N,P,T);
topt_maxiter = 300;
poblano_maxiter = 100;
topt_verbose = true;
t_start = tic;
[Xrec_topt_adapt,Result] = ...
    tucker_opt(size(X), X(ind), ind, rest, ...
    'Xinit', initial, 'maxiter', topt_maxiter, ...
    'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
toc(t_start)

er4 =  norm(tensor(Xrec_topt_adapt)-Y_real)/norm(Y_real)
Factor = Result.A;
A1 = Factor{1,1};
A2 = Factor{2,1};
A3 = Factor{3,1};
subspace(A1,U1)
subspace(A2,U2)
subspace(A3,U3)



%% IRTD (GOOOD, But SLOW)
% Max_iter = 10000;
% delta    = 0.1;
% t_start  = tic;
% Result3  = IRTD(Y,Omiga,Max_iter,[R1,R2,R3,R4],delta,Y_real);
% toc(t_start)
% 
% er3 = norm((Result3.Y-Y_real)) / norm(Y_real)
% 
% Factor = Result3.A;
% A1 = Factor{1,1};
% A2 = Factor{1,2};
% A3 = Factor{1,3};
% 
% subspace(A1,U1)
% subspace(A2,U2)
% subspace(A3,U3)

