function Result = HaLRTC(Y,Omiga,rou,Max_iter)
    X = Y;
    TS = size(Y);
    M = TS(1);
    N = TS(2);
    P = TS(3);
    alpha = 1;
    t = 1.15;
    Y1{1} = tensor(randn(M,N,P));
    Y1{2} = tensor(randn(M,N,P));
    Y1{3} = tensor(randn(M,N,P));
    M1{1} = tensor(zeros(M,N,P));
    M1{2} = tensor(zeros(M,N,P));
    M1{3} = tensor(zeros(M,N,P));
    index = 1 : 3;
    for iter = 1 : Max_iter
        X_o = tensor(zeros(M,N,P));
        % Update M
        for ii = 1 : 3
            T = double(tenmat(X + 1/rou*Y1{ii},ii));
            [U,S,V] = svd(T,'econ');
            S = max(S-alpha/rou,0);
            ind = sort(setdiff(index,ii),'ascend');
            M1{ii} = tensor(tenmat(U*S*V',ii,ind,[M,N,P]));
        end
        % Update X
        for ii = 1 : 3
            X_o = tensor(X_o) + M1{ii} - Y1{ii}/rou;
        end
        X = X.*Omiga + X_o.*~Omiga/3;
        % Update Y
        for ii = 1 : 3
            Y1{ii} = Y1{ii} - rou*(M1{ii} - X);
        end
        % Update rou
        rou = rou*t;
    end
    Result.Y = X;

end
