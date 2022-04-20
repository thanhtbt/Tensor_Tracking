function [PER,Factor] = PETRELS_SOAP(X,Omega,R,OPTS)

[I,J,T] = size(X);

%% Performance Evaluation
 
Factors    = OPTS.Factor_True;
Slide_true = OPTS.Slide_True;

PER.A  = zeros(1,T);
PER.B  = zeros(1,T);
PER.X  = zeros(1,T);
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
OPTS_PETRELS.lambda = 0.9;
Ht = khatrirao(Bt,At);

for t   =  t_train + 1 : T
    Xt      = X(:,:,t);
    xt      = Xt(:);
    idx     = find(xt);
      
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
    Ct       = [Ct; ct'];
    
    %% Performance Evaluation
   
    Factors_t  = Factors{t,1};
    A_true     = Factors_t{1,1};
    B_true     = Factors_t{1,2};
    Xt_true    = Slide_true{t,1};                
                    
    xt_true = Xt_true(:);
    xt_re   = Ht_es * ct;
    er      = xt_re - xt_true;
    
    [~,~,~,erA] = solve_perm_scale(At,A_true);
    [~,~,~,erB] = solve_perm_scale(Bt,B_true);
    
    PER.A(1,t)  = PER.A(1,t)  + erA;
    PER.B(1,t)  = PER.B(1,t)  + erB;
    PER.X(1,t)  = PER.X(1,t)  + norm(er)/norm(xt_true);
end

Factor.A = At;
Factor.B = Bt;
Factor.C = Ct;
end

