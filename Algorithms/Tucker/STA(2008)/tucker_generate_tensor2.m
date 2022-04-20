function [X,Factors] = tucker_generate_tensor(tensor_dim,tensor_rank,fac_noise,epsilon)
% Author     : LE Trung-Thanh
% Affiliation: University of Orleans, France
% Contact    : trung-thanh.le@univ-orleans.fr // letrungthanhtbt@gmail.com
% Date       : 4/2/2019

N = length(tensor_dim);
T = tensor_dim(N);
Tensor_Core  = tensor(randn(tensor_rank));
Factors      = cell(T,N); 
X            = cell(T,1);
for n = 1 : N-1
    Factors{1,n} = randn(tensor_dim(n),tensor_rank(n));
end
Factors{1,N} = randn(1,tensor_rank(N));
X1      = ttm(Tensor_Core,Factors(1,:),[1:N]);
X{1,1}  = tensor(X1.data);
for t = 2 : T
    for n = 1 : N-1
        Factors{t,n} = (1-epsilon(t)) * Factors{t-1,n} + epsilon(t)*randn(tensor_dim(n),tensor_rank(n));
    end
    Factors{t,N} = randn(1,tensor_rank(N));
    X_slice      = ttm(Tensor_Core,Factors(t,:),[1:N]); 
    X_slice      = X_slice.data + fac_noise*randn([tensor_dim(1:N-1)]);
    X{t,1}       = tensor(X_slice);
end

end
