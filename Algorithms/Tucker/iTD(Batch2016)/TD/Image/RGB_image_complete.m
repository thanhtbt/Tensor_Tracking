clear;clc;
%%
MP_set = [0.5,0.8,0.9];
%% Form tensor
Y2 = double(imread('lena.bmp'));
DimX = size(Y2);
Y2 = tensor(Y2);
nm = norm(Y2);
Y_real = Y2/norm(Y2)*sqrt(prod(DimX));
for per = 1 : length(MP_set)
    MP = MP_set(per);
    Omiga = tensor(double(rand(DimX)>MP));
    MM = sum(sum(double(tenmat(~Omiga,1))));
    Y = Y_real.*Omiga;
    %% Completion with LRTI
        Max_iter = 400;
        tic;
        Result1 = LRTI(Y,Omiga,Max_iter,5,Y_real,nm);
        rtime1 = toc;
        MSE1 = norm(~Omiga.*(Y_real-tensor(ktensor(Result1.A)))*nm/sqrt(prod(DimX))/255).^2/(MM);
    %% Completion with HaLRTC
        Max_iter = 100;
        rou = 1e-5;
        t0 = cputime;
        Result2 = HaLRTC(Y,Omiga,rou,Max_iter);
        rtime2 = cputime-t0;
        MSE2 = norm(~Omiga.*(Y_real-Result2.Y)*nm/sqrt(prod(DimX))/255).^2/(MM);
    %% Completion with IRTD
        Max_iter = 400;
        delta = 0.5;
        t0 = cputime;
        Result3 = IRTD(Y,Omiga,Max_iter,DimX,delta,Y_real);
        rtime3 = cputime-t0;
        MSE3 = norm(~Omiga.*(Y_real-Result3.Y)*nm/sqrt(prod(DimX))/255).^2/(MM);
    %% WTucker
        X = double(Y);
        ind = find(double(Omiga));
        rest = [min(DimX(1),15),min(DimX(2),15),min(DimX(3),15)];
        initial = randn(DimX);
        topt_maxiter = 300;
        poblano_maxiter = 100;
        topt_verbose = true;
        t0 = tic;
        [Xrec_topt_adapt, topt_ranks, topt_errors, toptC, toptA] = ...
    tucker_opt(size(X), X(ind), ind, rest, ...
    'Xinit', initial, 'maxiter', topt_maxiter, ...
    'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
         rtime4 = toc;
         MSE4 = norm(~Omiga.*(Y_real-Xrec_topt_adapt)*nm/sqrt(prod(DimX))/255).^2/(MM);
end
