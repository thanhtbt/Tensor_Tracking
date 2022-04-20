function [Factor_Es,PER] = ACP(X_cell,Omega,cp_rank,OPTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author     : LE Trung Thanh 
% Contact    : letrungthanhtbt@gmail.com or trung-thanh.le@univ-orleans.fr 
% Reference:  [1] Le Trung Thanh, Karim Abed-Meraim, Nguyen Linh-Trung and Adel Hafiane
%             "A fast randomized adaptive CP decomposition	for streaming tensors"
%             IEEE-ICASSP, 2021. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <= 3
    flag = 0;     % without performance estimation part
    lambda  = 0.5;
    Factors = []; PER = [];
else
    flag = 1;      % with performance estimation part
    Slide_True   = OPTS.Slide_True;
    Factor_True  = OPTS.Factor_True; 
    if isfield(OPTS,'lambda')
         lambda = OPTS.lambda;
    else lambda = 0.5;
    end
    
end

R          = cp_rank;
T          = length(X_cell);
tensor_dim = size(X_cell{1,1});
N          = length(tensor_dim) + 1;
PER        = zeros(N,T);

%% Random Initalization
N = length(tensor_dim);
Factor_es = cell(1,N+1);
for n = 1 : N
    Factor_Es{1,n}    = randn(tensor_dim(n),R);
    Factor_Delta{1,n} = zeros(tensor_dim(n),R);
    S{1,n}            = 100*eye(R);
end

Factor_Es{1,N+1} = [];
last_fac         = [];
K    = round(10*R*log(prod(tensor_dim)));
Ht   = khatrirao(Factor_Es(1,1:N),'r');
beta = lambda;

for t = 1 : T
    
    Xt  = X_cell{t,1};
    xt  = Xt(:); 
 
    Omega_t  = Omega{t,1};
    idx      = find(Omega_t(:));
    xt_Omega = xt(idx,:);
    Ht_Omega = Ht(idx,:); 
    L        = length(idx);
    if (L > K) 
        idx_re = randperm(L,K); 
        dt     = Ht_Omega(idx_re,:) \ xt_Omega(idx_re,:);
    else
        dt     = Ht_Omega \ xt_Omega;
    end
    last_fac  = [last_fac; dt'];
    
    X_re = Ht * dt;
    X_re = reshape(X_re,tensor_dim(1:N));
    ER   = Omega_t .* (Xt - X_re);
    ER   = tensor(ER);
    
    Factor_Es{1,N+1} = dt';
    for n = 1 : N  
        ER_unfolding_n  = ten2mat(ER,[n]);      
        X_re_n          = ten2mat(tensor(Xt),[n]);
        W_n             = pinv(Factor_Es{1,n}) * X_re_n;
        S{1,n}          = W_n * W_n' + beta*S{1,n};
        V_n             = S{1,n} \ W_n;  
        Factor_Es{1,n}  = Factor_Es{1,n} + ER_unfolding_n * V_n';
    end
    Ht  = khatrirao(Factor_Es(1,1:N),'r');

    %% Performance Evaluation Part 
    if flag == 1
        for n = 1 : N
           A_es = Factor_Es{1,n};
           A_tr = Factor_True{t,1}{1,n};
           [~,~,~,erU] = solve_perm_scale(A_es,A_tr);
           PER(n,t)    = erU;
        end   
        X_re = Ht * dt;
        X_tr = OPTS.Slide_True{t,1};
        Er = X_re - X_tr(:); 
        PER(N+1,t) = norm(Er) / norm(X_tr(:));
    else
    end
end
Factor_Es{1,N+1} = last_fac; 
end

