clear;clc;
M = 32;
N = 32;
P = 32; 
R = 6;
Max_count = 100;
SNR = [10,20,30,100];
Rv1 = zeros(8,Max_count,length(SNR));
Rv2 = zeros(6,Max_count,length(SNR));
Rv3 = zeros(12,Max_count,length(SNR));
Rv4 = zeros(6,Max_count,length(SNR));
MP_set = [0.5,0.8];
Lambda_set = 0.2*ones(4,1);%[0.2,1,5,10];
for count = 1 : Max_count
    for SNR_count = 1 : length(SNR)
        [count,SNR(SNR_count)]
        Sigma = 10^(-SNR(SNR_count)/20);
        %% Form tensor
        U1=randn(M,R);
        U2=randn(N,R);
        U3=randn(P,R);
        gamma = ones(R,1);
        Y_real = tensor(zeros(M,N*P),[M N P]);
        for ii = 1 : R
            Y_real = Y_real + ttm(tensor(gamma(ii),[1 1 1]),{U1(:,ii),U2(:,ii),U3(:,ii)},[1:3]);
        end
        Y_real = Y_real/norm(Y_real)*sqrt(M*N*P);
        noise = tensor(randn(M,N,P));
        Y_n = Y_real + Sigma*norm(Y_real)/norm(noise)*noise;
        for per = 1 : 2
            [count,SNR(SNR_count),MP_set(per)]
            MP = MP_set(per);
            Omiga = tensor(rand(M,N,P)>MP);
            Y = Y_n.*Omiga;
            %% Completion with LRTI
            Max_iter = 10000;
            tic;
            Lambda = Lambda_set(SNR_count);
            Result1 = LRTI(Y,Omiga,Max_iter,Lambda);
            Rv1(1+(per-1)*4,count,SNR_count) = toc;
            Rv1(2+(per-1)*4,count,SNR_count) = norm(Y_real-tensor(ktensor(Result1.A)))/norm(Y_real);
            Rv1(3+(per-1)*4,count,SNR_count) = Result1.rank;
            Rv1(4+(per-1)*4,count,SNR_count) = norm((Y_real-tensor(ktensor(Result1.A))).*~Omiga)/norm(Y_real);
            %% Completion with HaLRTC
            Max_iter = 100;
            rou = 1e-5;
            tic;
            Result2 = HaLRTC(Y,Omiga,rou,Max_iter);
            Rv2(1+(per-1)*3,count,SNR_count) = toc;
            Rv2(2+(per-1)*3,count,SNR_count) = norm(Y_real-Result2.Y)/norm(Y_real);
            Rv2(3+(per-1)*3,count,SNR_count) = norm((Y_real-Result2.Y).*~Omiga)/norm(Y_real);
            %% Completion with IRTD
            Max_iter = 10000;
            delta = 0.1;
            tic;
            Result3 = IRTD(Y,Omiga,Max_iter,[M,N,P],delta,Y_real);
            Rv3(1+(per-1)*6,count,SNR_count) = toc;
            Rv3(2+(per-1)*6,count,SNR_count) = norm((Result3.Y-Y_real))/norm(Y_real);
            Rv3(3+(per-1)*6,count,SNR_count) = size(Result3.A{1},2);
            Rv3(4+(per-1)*6,count,SNR_count) = size(Result3.A{2},2);
            Rv3(5+(per-1)*6,count,SNR_count) = size(Result3.A{3},2);
            Rv3(6+(per-1)*6,count,SNR_count) = norm((Y_real-Result3.Y).*~Omiga)/norm(Y_real);
            %% WTucker
            X = double(Y);
            ind = find(double(Omiga));
            rest = 2*[R,R,R];
            initial = randn(M,N,P);
            topt_maxiter = 300;
            poblano_maxiter = 100;
            topt_verbose = true;
            t0 = cputime;
            [Xrec_topt_adapt, topt_ranks, topt_errors, toptC, toptA] = ...
        tucker_opt(size(X), X(ind), ind, rest, ...
        'Xinit', initial, 'maxiter', topt_maxiter, ...
        'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
            Rv4(1+(per-1)*3,count,SNR_count) = cputime - t0;
            Rv4(2+(per-1)*3,count,SNR_count) = norm(tensor(Xrec_topt_adapt)-Y_real)/norm(Y_real);
            Rv4(3+(per-1)*3,count,SNR_count) = norm((tensor(Xrec_topt_adapt)-Y_real).*~Omiga)/norm(Y_real);
        end
    end
end
