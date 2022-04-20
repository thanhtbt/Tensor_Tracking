function X_stream = online_tensor_generator(Size_Fixed,Rank,Num_Slides,epsilon)
%% At each time instant t
% X_stream{1,:}: Slides of Tensor
% X_stream{2,:}: Factors of Tensor
% X_t = U_1 o U_2 o ... o U_{n-1} o U_n(t,:)
% factors are changing slowly or static
% U(t) = U(t-1) + epsilon*randn(Size);


%% Main Algorithms
L = length(Size_Fixed);
Factors = cell(1,L+1);
for n = 1 : L
    Factors{1,n} = randn(Size_Fixed(n),Rank);
end

N = L+1; % Order of tensor
Factors{1,N} = [];


for t = 1 : Num_Slides
    Factors{1,N}  = randn(1,Rank);
    Xt_mat        = khatrirao(Factors,'r');
    Xt            = sum(Xt_mat')';
    X_stream{t,1} = reshape(Xt,Size_Fixed);
    X_stream{t,2} = Factors;
      
    p = mod(t + Rank - 2, Rank - 1) + 1;
    Q = eye(Rank);
     
    alpha = epsilon(t);
    T = [cos(alpha) -sin(alpha);
        sin(alpha)   cos(alpha)];
    Q(p:p+1,p:p+1) = T;
    
    if alpha  == 1
        for n = 1 : L
            Factors{1,n} =  randn(Size_Fixed(n),Rank);;
        end
    else
        for n = 1 : L
            Factors{1,n} =   Factors{1,n} * Q; 
        end
    end
end







