function Result = LRTI(Y,Omiga,Max_iter,Lambda)   
    I = size(Y);
    Nway = length(I);
    %initial
    R = max(I);%init_rank
    eps = 1e-10;
    for ii = 1 : Nway
        Yn = double(tenmat(Y,ii));
        [U,S,~] = svd(Yn,'econ');
        US = U*sqrt(S);
        L = size(US,2);
        if L>=R
            A{ii} = US(:,1:R);
        else
            A{ii}(:,1:L) = US;
            A{ii}(:,L+1:R) = randn(I(ii),R-L) + 1i*randn(I(ii),R-L);
        end
    end
    E_X = tensor(zeros(I));
    for iter = 1 : Max_iter
%         [iter,R]
        EX_old = E_X;
        for ii = 1 : Nway
            gamma{ii} = (sum(abs(A{ii}).^2)).';
            ind = sort(setdiff(Nway:-1:1,ii),'descend');
            B = ones(1,R);
            for jj = ind
                B = KR_P(B,A{jj});
            end
            Y_n = double(tenmat(Y,ii));
            Index = logical(double(tenmat(Omiga,ii)));
            for jj = 1 : I(ii)
                A{ii}(jj,:) = Lambda*Y_n(jj,Index(jj,:))*conj(B(Index(jj,:),:))/(Lambda*conj(B(Index(jj,:),:)'*B(Index(jj,:),:)) + eye(size(B,2)));
            end
        end
        % Pruning
        gamma_t = zeros(R,1) + eps;
        for ii = 1 : Nway
            gamma_t = gamma_t + gamma{ii};
        end
        ind = gamma_t > max(gamma_t)/50;
        for ii = 1 : Nway
            A{ii} = A{ii}(:,ind);
        end
        gamma_t = gamma_t(ind);
        R = length(gamma_t);
        E_X = cktensor(A,ones(R,1));
        if norm(EX_old - E_X)<1e-3
            break;
        end
    end
    Result.Y = cktensor(A,ones(R,1));
    Result.A = A;
    Result.gamma = gamma;
    Result.rank = length(gamma{1});
    Result.Lambda = Lambda;
end
 function Y = KR_P(A,B)
    [I,T] = size(A);
    J = size(B,1);
    Y = zeros(I*J,T);
    for jj = 1 : T
        Y(:,jj) = kron(A(:,jj),B(:,jj));
    end
 end
function X = cktensor(A,lambda)
    [~,Nway] = size(A);
    R = length(lambda);
    X = tendiag(lambda,R*ones(1,Nway));
    X = ttm(X,{A{1:Nway}});
end