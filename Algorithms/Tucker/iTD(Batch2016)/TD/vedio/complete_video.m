clear;clc;
%%
MP_set = [0.5,0.8,0.9];
%% Form tensor
load VideoData.mat
Y = Data;
TS = size(Y);
Y = tensor(Y);
mean_Y = sum(sum(double(tenmat(Y,1))))/(prod(TS));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
nm = norm(Y);
Y_real = Y/norm(Y)*sqrt(prod(TS));%/128;
for per = 1 : length(MP_set)
    MP = MP_set(per);
    Omiga = tensor(double(rand(TS)>MP));
    MM = sum(sum(double(tenmat(~Omiga,1))));
    Y = Y_real.*Omiga;
    %% Completion with LRTI
        Max_iter = 400;
        tic;
        Result1 = LRTI(Y,Omiga,Max_iter,3,Y_real,nm);
        rtime1 = toc;
        MSE1 = norm(~Omiga.*(Y_real-tensor(ktensor(Result1.A)))*nm/sqrt(prod(TS))/255).^2/(MM);
    %% Completion with HaLRTC
        Max_iter = 100;
        rou = 1e-5;
        t0 = cputime;
        Result2 = HaLRTC_new(Y,Omiga,rou,Max_iter,Y_real);
        rtime2 = cputime-t0;
        MSE2 = norm(~Omiga.*(Y_real-Result2.Y)*nm/sqrt(prod(TS))/255).^2/(MM);
    %% Completion with IRTD
        Max_iter = 100;
        delta = 1;
        t0 = cputime;
        Result3 = IRTD(Y,Omiga,Max_iter,TS,delta,Y_real);
        rtime3 = cputime-t0;
        MSE3 = norm(~Omiga.*(Y_real-Result3.Y)*nm/sqrt(prod(TS))/255).^2/(MM);
    %% WTucker
        X = double(Y);
        ind = find(double(Omiga));
        rest = [min(TS(1),50),min(TS(2),40),min(TS(3),100),min(TS(4),20)];
        initial = randn(TS);
        topt_maxiter = 300;
        poblano_maxiter = 100;
        topt_verbose = true;
        t0 = tic;
        [Xrec_topt_adapt, topt_ranks, topt_errors, toptC, toptA] = ...
    tucker_opt(size(X), X(ind), ind, rest, ...
    'Xinit', initial, 'maxiter', topt_maxiter, ...
    'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
         rtime4 = toc;
         MSE4 = norm(~Omiga.*(Y_real-Xrec_topt_adapt)*nm/sqrt(prod(TS))/255).^2/(MM);
end
