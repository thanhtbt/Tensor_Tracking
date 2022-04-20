function [Factor_Es,PER] = ATD(X_cell,Omega,tucker_rank,OPTS)

%% Missing

if nargin <= 3
    flag = 0; % without performance estimation part
    lambda  = 0.7;
    Factors = []; PER = [];
else
    flag = 1; % with performance estimation part
    if isfield(OPTS,'Factor_True')
        Factor_True  = OPTS.Factor_True;
        flag_factor = 1;
    else
        flag_factor = 0;
    end
    if isfield(OPTS,'Slide_True')
        Slide_True = OPTS.Slide_True;
        flag_slide = 1;
    else
        flag_slide = 0;
    end
    if isfield(OPTS,'lambda')
         lambda = OPTS.lambda;
    else lambda = 0.7;
    end
    
end

T   = length(X_cell);
tensor_dim = size(X_cell{1,1});
N   = length(tensor_dim) + 1;
PER = zeros(N,T);
PER = zeros(N+1,T);

%% Random Initialization
Tensor_Core = tensor(randn(tucker_rank));
Factor_Es = cell(1,N);
S         = cell(1,N-1);
for ii = 1 : N-1
    Factor_Es{1,ii} = randn(tensor_dim(ii),tucker_rank(ii));
    S{1,ii}         = 100*eye(tucker_rank(ii));
end
Factor_Es{1,N} = [];


last_fac = [];
r_N = tucker_rank(N);
K = round(10*r_N*log(r_N)); 
R = prod(tucker_rank);
for t = 1 : T
    t_start = tic;
    Xt      = X_cell{t,1};
    Omega_t = Omega{t,1};
    Xt_data  = Xt.data;
    xt       = Xt_data(:);
    idx     = find(Omega_t);
    L       = length(idx);
    if L > K && L / length(xt) > 0.5
        idx_sub  = randsample(L,K);
    else
         idx_sub  = 1:L;
    end
    % idx_sub  = 1:L;
    idx_t    = idx(idx_sub);
   
    
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
        X_re_n          = ten2mat(Xt,[n]);
        
        W_n    = pinv(Factor_Es{1,n}) * X_re_n;
        S{1,n} = W_n * W_n' + lambda*S{1,n};
        V_n    = S{1,n} \ W_n;
        Factor_Es{1,n}   = Factor_Es{1,n} + ER_unfolding_n * V_n';
   
      
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
    
    t_end = toc(t_start);
    
    %% Performance Evaluation Part
    
    if flag == 1
        Core_Matrix = ten2mat(Tensor_Core,1);
        X_true = Slide_True{t,1};
        X_re   = Factor_Es{1,1} * Core_Matrix * U;
        X_true = ten2mat(X_true,1);
        PER(N,t) = norm(X_re - X_true) / norm(X_true(:)); 
        for n = 1 : N-1
            PER(n,t) = sub_est_per(Factor_Es{1,n},Factor_True{t,n},'SEP');
        end
        PER(N+1,t)   = t_end;
    else
        
    end
    
end
Factor_Es{1,N} = last_fac;
end