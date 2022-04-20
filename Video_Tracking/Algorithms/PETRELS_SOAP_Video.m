function [PER,Xre] = PETRELS_SOAP(X,Omega,R,OPTS)

[I,J,T] = size(X);


%% Performance Evaluation

 
Xtrue      = OPTS.Xtrue;
PER        = zeros(1,T);


%% Batch Initialization
t_train = max(I,J);
X_train = X(:,:,1:t_train);
[P_train,~,~] = cp_wopt_mod(tensor(X_train),Omega(:,:,1:t_train),R);
At =  P_train.U{1}; % randn(I,R); 
Bt =  P_train.U{2}; % randn(J,R);
Ct =  P_train.U{3}; % randn(t_train,R);

%% Main Algorithm

Rinv = repmat(100*eye(R),1,I*J);
OPTS_PETRELS.Rinv = Rinv;
OPTS_PETRELS.lambda = 0.999;
Ht = khatrirao(Bt,At);

for t   =  1 : T
    Xt      = X(:,:,t);
    xt      = Xt(:);
    idx     = find(xt);
    idx_m     = find(~xt); 
    %% Estimate ct and Ht = khatrirao(Bt,At) using PETRELS 
    [Ht,ct,OPTS_PETRELS] = Petrels_Subspace(xt,Ht,OPTS_PETRELS);
    
    %% Estimate At, Bt from Ht
    A1 = zeros(I,R);
    B1 = zeros(J,R);
    for r = 1:R
        Hr = reshape(Ht(:,r),I,J);       % ar is the left principal sing vector and conj(cr) the left one
        b1 = Hr'*At(:,r);
        B1(:,r) = conj(b1);
        a1 = Hr*b1;
        A1(:,r) = a1/norm(a1,'fro');     % Normalize it because sing value was included in c1
    end
    At  = A1;
    Bt  = B1;
    %% Re-estimate ct
    Ht_es    = khatrirao(Bt,At);
    ct       = Ht_es(idx,:) \ xt(idx,:);
    Ct = [Ct; ct'];
    
    %% Performance Evaluation
   

    X_true = Xtrue(:,:,t);                 
    xt_re   = Ht_es * ct;
    x_re    =  xt;
    x_re(idx_m) = xt_re(idx_m);
    X_re    = reshape(x_re,[I J]);
    
    Xre(:,:,t) = X_re;
    ER        = X_re(idx_m) - X_true(idx_m);
    PER(1,t)  =  norm(ER)/norm(X_true(idx_m));
    
end

Factor.A = At;
Factor.B = Bt;
Factor.C = Ct;
end


function [U,coeff,OPTS] = Petrels_Subspace(xt,U_old,OPTS)
if isfield(OPTS,'lambda'), % Forgetting factor
    lambda = OPTS.lambda;
else lambda = 0.98;
end
Rinv = OPTS.Rinv;

r = size(U_old,2);

%% Main Algorithm

U = U_old;
idx = find(xt);
x_omega = xt(idx,:);
U_omega = U_old(idx,:); 

coeff = U_omega \ x_omega;
x_re  = U_omega * coeff ;
residual = x_omega - x_re;

for ii=1:length(idx) 
    Tinv = Rinv(:,(idx(ii)-1)*r+1 : idx(ii)*r); 
    ss   = Tinv *coeff;
    Tinv = Tinv-ss*ss'*1/(lambda + coeff'*Tinv*coeff);
    U(idx(ii),:) =  U_omega(ii,:) + lambda^(-1) * residual(ii) * coeff'* Tinv;
    Rinv(:,(idx(ii)-1)*r+1 : idx(ii)*r) = Tinv;
end

Rinv = lambda^(-1) * Rinv;

OPTS.Rinv = Rinv;
end






