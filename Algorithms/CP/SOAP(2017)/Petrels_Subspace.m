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




