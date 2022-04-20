clear;clc;
%%
Max_count = 100;
SNR = [10,20,30,100];
Rv1 = zeros(8,Max_count,length(SNR));
Rv2 = zeros(6,Max_count,length(SNR));
Rv3 = zeros(12,Max_count,length(SNR));
Rv4 = zeros(6,Max_count,length(SNR));
MP_set = [0.5,0.8];
Lambda_set = 0.2*ones(4,1);
for count = 1 : Max_count
    %% Form tensor
   load amino.mat
    Y_real = tensor(reshape(X,DimX));
    for SNR_count = 1 : length(SNR)
        Sigma = 10^(-SNR(SNR_count)/20);
        Y_real = Y_real/norm(Y_real)*sqrt(prod(DimX));
        noise = tensor(randn(DimX));
        Y_t = Y_real + Sigma*noise/norm(noise)*norm(Y_real);
        for per = 1 : 2
            [count,SNR(SNR_count),MP_set(per)]
            MP = MP_set(per);
            Omiga = tensor(double(rand(DimX)>MP));
            Y = tensor(Y_t.*Omiga);
            %% Completion with LRTI
            Max_iter = 5000;
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
            Max_iter = 5000;
    %         Max_iter = 3000;
            tic;
            delta = 0.1;
            init_size = DimX;
            Result3 = IRTD(Y,Omiga,Max_iter,init_size,delta,Y_real);
            Rv3(1+(per-1)*6,count,SNR_count) = toc;
            Rv3(2+(per-1)*6,count,SNR_count) = norm((Result3.Y-Y_real))/norm(Y_real);
            Rv3(3+(per-1)*6,count,SNR_count) = size(Result3.A{1},2);
            Rv3(4+(per-1)*6,count,SNR_count) = size(Result3.A{2},2);
            Rv3(5+(per-1)*6,count,SNR_count) = size(Result3.A{3},2);
            Rv3(6+(per-1)*6,count,SNR_count) = norm((Y_real-Result3.Y).*~Omiga)/norm(Y_real);
        %% WTucker
            X = double(Y);
            ind = find(double(Omiga));
            rest = [min(DimX(1),10),min(DimX(2),10),min(DimX(3),10)];
            initial = randn(DimX);
            topt_maxiter = 300;
            poblano_maxiter = 100;
            topt_verbose = true;
            t0 = tic;
            [Xrec_topt_adapt, topt_ranks, topt_errors, toptC, toptA] = ...
        tucker_opt(size(X), X(ind), ind, rest, ...
        'Xinit', initial, 'maxiter', topt_maxiter, ...
        'poblanomax', poblano_maxiter, 'verbose', topt_verbose);
             Rv4(1+(per-1)*3,count,SNR_count) = toc;
             Rv4(2+(per-1)*3,count,SNR_count) = norm(tensor(Xrec_topt_adapt)-Y_real)/norm(Y_real);
             Rv4(3+(per-1)*3,count,SNR_count) = norm((tensor(Xrec_topt_adapt)-Y_real).*~Omiga)/norm(Y_real);
        end
    end
end
