function Result = HaLRTC(Y,Omiga,rou,Max_iter)
    X = Y;
    TS = size(Y);
    Nway = length(TS);
    alpha = 1;
    t = 1.15;
    for ii = 1 : Nway
        Y1{ii} = tensor(randn(TS));
        M1{ii} = tensor(zeros(TS));
    end
    index = 1 : Nway;
    for iter = 1 : Max_iter
        X_o = tensor(zeros(TS));
        % Update M
        for ii = 1 : Nway
            T = double(tenmat(X + 1/rou*Y1{ii},ii));
            [U,S,V] = svd(T,'econ');
            S = max(S-alpha/rou,0);
            ind = sort(setdiff(index,ii),'ascend');
            M1{ii} = tensor(tenmat(U*S*V',ii,ind,TS));
        end
        % Update X
        for ii = 1 : Nway
            X_o = tensor(X_o) + M1{ii} - Y1{ii}/rou;
        end
        X = X.*Omiga + X_o.*~Omiga/Nway;
        % Update Y
        for ii = 1 : Nway
            Y1{ii} = Y1{ii} - rou*(M1{ii} - X);
        end
        % Update rou
        rou = rou*t;
    end
    Result.Y = X;
end