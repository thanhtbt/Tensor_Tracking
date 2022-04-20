function [PER,Xre] = ACP(X_video,Omega,R,OPTS)
% "Randomized Adaptive CP Decomposition for Three-way Tensors"
% Author: LE Trung Thanh
% Contract: trung-thanh.le@univ-orleans.fr
% Last UpDate: 6/1/2020

flag == 1;
Xtrue  = OPTS.Xtrue;
lambda = OPTS.lambda;
%% 
tensor_dim = size(X_video);
T          = tensor_dim(end);

N          = length(tensor_dim) - 1; 
PER        = zeros(1,T);

%% Random Initalization
Factor_es = cell(1,N+1);
for n = 1 : N
    Factor_Es{1,n}    = randn(tensor_dim(n),R);
    Factor_Delta{1,n} = zeros(tensor_dim(n),R);
    S{1,n} = 100*eye(R);
end
Factor_Es{1,N+1} = [];
last_fac = [];
% K  = round(10*R*log( prod(tensor_dim(1:end-1))));
K  = round(10*R*log(  R ));
Ht  = khatrirao(Factor_Es(1,1:N),'r');


for t = 1 : T
    
    Xt  = X_video(:,:,t); 
    xt  = Xt(:); 
    
    Omega_t = Omega(:,:,t);
    idx = find(Omega_t);
    idx_m = find(~Omega_t);
    
    
    xt_Omega = xt(idx,:);
    Ht_Omega = Ht(idx,:); 
    L        = length(idx);
    if L/(length(xt)) >= 0.5  ||   L > K 
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
        W_n    = pinv(Factor_Es{1,n}) * X_re_n;
        S{1,n} = W_n * W_n' + 0.01*eye(R); 
        V_n    = S{1,n} \ W_n;
        Factor_Es{1,n}   = Factor_Es{1,n} + ER_unfolding_n * V_n';
        
        
    end
    Ht  = khatrirao(Factor_Es(1,1:N),'r');

    %% Performance Estimation
    X_re     = Ht * dt;
    X_re     = reshape(X_re,[tensor_dim(1:2)]);
    X_true   = Xtrue(:,:,t);
    PER(1,t) = norm(X_re(:) - X_true(:)) / norm(X_true(:));
    
    X_re2        = Xt;
    X_re2(idx_m) = X_re(idx_m);
    Xre(:,:,t)   = X_re2;
end
Factor_Es{1,N+1} = last_fac;

end

