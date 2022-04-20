function [PER,Xre] = ATD(X_video,Omega,tucker_rank,OPTS)

%% For Evaluation
flag == 1;
Xtrue  = OPTS.Xtrue;
lambda = OPTS.lambda;
%% 
tensor_dim = size(X_video);
T          = tensor_dim(end);

N          = length(tensor_dim); % order
PER        = zeros(1,T);

%% Random Initialization
Tensor_Core = tensor(randn(tucker_rank));
Factor_Es   = cell(1,N);
S           = cell(1,N-1);
for ii = 1 : N-1
    Factor_Es{1,ii} = randn(tensor_dim(ii),tucker_rank(ii));
    S{1,ii}         = 100*eye(tucker_rank(ii));
end
Factor_Es{1,N} = [];

%% 
last_fac = [];
r_N = tucker_rank(N);
K = round(10*r_N*log(prod(tensor_dim(1:end-1)))); 
R = prod(tucker_rank);


Xre = zeros(tensor_dim);
for t = 1 : T
    
    Xt      = X_video(:,:,t);
    Omega_t = Omega(:,:,t);
    idx     = find(Omega_t);
    idx_m   = find(~Omega_t);
    L       = length(idx);
    if L > K 
         idx_sub  = randsample(L,K);
    else
         idx_sub  = 1:L;
    end
    idx_sub  = 1:L;
    idx_t    = idx(idx_sub);
    xt       = Xt(:);
    
    Ht        = ttm(Tensor_Core,Factor_Es(1,[1:N-1]),[1:N-1]);
    Ht_matrix = ten2mat(Ht,N)';
       
    Ht_Omega  = Ht_matrix(idx_t,:);
    xt_Omega  = xt(idx_t,:);
       
    dt        =  Ht_Omega \ xt_Omega; 
    last_fac  = [last_fac;  dt'];
    
    delta_x   = Omega_t(:) .* (xt - Ht_matrix * dt);
    ER        = reshape(delta_x,tensor_dim(1:N-1)); 
    ER        = tensor(ER);
    
    % Factor_Es{1,N} = dt';
    %% Estimation of Factors
    for n =  N-1:-1:1
        ER_unfolding_n  = ten2mat(ER,[n]);
        X_re_n          = ten2mat(tensor(Xt),[n]);
        
        W_n    = pinv(Factor_Es{1,n}) * X_re_n;
        S{1,n} = W_n * W_n' + lambda*S{1,n} + 0.001*eye(tucker_rank(n));
        V_n    = S{1,n} \ W_n;
        Factor_Es{1,n}   = Factor_Es{1,n} + ER_unfolding_n * V_n';    
%         U_tem  = Factor_Es{1,n} + ER_unfolding_n * V_n';
%         Factor_Es{1,n}   = U_tem * (U_tem' * U_tem)^(-0.5);
    end
    
    %% Estimation of Core
    U  = dt;
    for n = N-1 : -1 : 2
        U  = kron(U,(Factor_Es{1,n})');
    end
    ER_Core = pinv(Factor_Es{1,1}) *  ER_unfolding_n * pinv(U);
    Core_Matrix = ten2mat(Tensor_Core,1);
    Core_Matrix = Core_Matrix + ER_Core;
    Core = mat2ten(Core_Matrix,[tucker_rank],1);
    Tensor_Core = tensor(Core);
    
    
    %% Performance Estimation
    
    
    Core_Matrix = ten2mat(Tensor_Core,1);
    X_re     = Factor_Es{1,1} * Core_Matrix * U;
    X_true   = Xtrue(:,:,t); 
    PER(1,t) = norm(X_re(idx_m) - X_true(idx_m)) / norm(X_true(idx_m));
    
    X_re2        = Xt;
    X_re2(idx_m) = X_re(idx_m);
    Xre(:,:,t)   = X_re2;
    
    
end
Factor_Es{1,N} = last_fac;
end