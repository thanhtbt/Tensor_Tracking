function [U_es,P,D,err] = perm_scale_solve(U_hat,U_tr)

[I_tr,R_tr]   = size(U_tr);
[I_hat,R_hat] = size(U_hat);
if I_tr ~= I_hat || R_tr ~= R_hat
    error('permsolve:InvalidSize',['Input matrices must have the same size']);
end
U_in = U_tr;
R    = R_tr;
% Normalization of all columns
if det(U_hat' * U_hat) > 0
    for r = 1 : R
        U_hat(:,r) = U_hat(:,r) / norm(U_hat(:,r));
        U_tr(:,r)  = U_tr(:,r) / norm(U_tr(:,r));
    end 
    prod_init = abs(pinv(U_hat)*U_tr);
    prod      = prod_init;
    c1 = 1 : R;  p1 = c1;
    c2 = 1 : R;  p2 = c2;
   
    for r = 1 : R
        [row,col]  = find(prod==max(max(prod)));
        row=row(1); 
        col=col(1);
        p1(r)   = c1(col);
        p2(r)   = c2(row);
        c1(col) = [];
        c2(row) = [];
        prod(row,:) = [];
        prod(:,col) = [];
    end
    [temp p2] = sort(p2,'ascend');
    perm_col  = p1(p2);
    
    % Permutation 
    P = eye(R);
    P = P(:,perm_col);
    % permuted estimate
    U_hat(:,perm_col) = U_hat;
  
    % Scaling 
    scal=zeros(1,R); % vectors to store scaling factors
    for r=1:R
        scal(r) = (U_in(:,r)'*U_in(:,r))/(U_in(:,r)'*U_in(:,r));
    end
    U_es =  U_hat * diag(scal);
    D    =  diag(1./scal);
    err  = norm(U_es - U_in,'fro')/norm(U_in,'fro');
    
    
else
    err  = norm(U_hat - U_in,'fro')/norm(U_in,'fro');
    U_es = [];
    P    = [];
    D    = [];
end




