function X_stream = cp_nft_generator(Size_Fixed,Rank,Num_Slides,epsilon)
%% At each time instant t
% X_stream{1,:}: Slides of Nonnegative Tensor
% X_stream{2,:}: Factors of Nonnegative Tensor
% X_t = U_1 o U_2 o ... o U_{n-1} o U_n(t,:)
% factors are changing slowly or static
% U(t) = U(t-1) + epsilon*randn(Size);


%% Main Algorithms
L = length(Size_Fixed);
Factors = cell(1,L+1);
for n = 1 : L
     Factors{1,n} = max(0,randn(Size_Fixed(n),Rank)); % Possitive Sparse
end

N = L+1; % Order of tensor
Factors{1,N} = [];


for t = 1 : Num_Slides
    Factors{1,N}  = max(0,randn(1,Rank)); 
    Xt_mat        = khatrirao(Factors,'r');
    Xt            = sum(Xt_mat')';
    X_stream{t,1} = reshape(Xt,Size_Fixed);
    X_stream{t,2} = Factors;
    
    sigma = epsilon(1,t);
    if sigma == 1 %% an abrupt change
        for n = 1 : L
             Factors{1,n} = max(0,randn(Size_Fixed(n),Rank));
        end
    else
        for n = 1 : L 
            Noise        =  randn(Size_Fixed(n),Rank);
            Factors{1,n} =  Factors{1,n} + sigma*Noise;
            Factors{1,n} =  max(0, Factors{1,n});
        end
    end

end







