function Result = IRTD(Y,Omiga,Max_iter,init_size,delta,Y_real)
%--------------------------------------------------------------------------
% Author : Linxiao Yang
% Date : 10/9/2015
% Intro : Yang, Linxiao, et al. 
%         "An iterative reweighted method for tucker decomposition of incomplete tensors."
%         IEEE Transactions on Signal Processing 64.18 (2016): 4817-4829. 
%--------------------------------------------------------------------------
%% Initialization
    TS = size(Y);
    Nway = length(TS);
    C_S = init_size;
    % parameter setting
    a = 0;
    b = 1e-10;
    sig = 0.1;
    FISTA_count = 2;
    C = 100;
    M = sum(sum(double(tenmat(~Omiga,1))));
    % initializing core tensor and factor matrices
    for ii = 1 : Nway
        Y_m = double(tenmat(Y,ii));
        [U,S,V] = svd(Y_m,'econ');
        A{ii} = U(:,1:C_S(ii));
    end
    X = ttm(Y,{A{1:Nway}},'t');
    X2 = X.^2;
    for ii = 1 : Nway
        X_m = double(tenmat(X2,ii));
        alpha{ii} = (1 + a)./(sum(X_m,2)+b);
    end
    Y_E = tensor(zeros(TS));
%% Main Loop
    for iter = 1 : Max_iter
        Y_E_old = Y_E;
        % Update X
        t = 1;
        % Form D
        TTemp = tensor(ones(C_S));
        D = tensor(zeros(C_S));
        for ii = 1 : Nway
            D = D + ttm(TTemp,diag(alpha{ii}),ii);
        end
        % Caculate Lipschitz constant
        L_x = 1;
        for ii = 1 : Nway
            A_t = svd(A{ii},'econ');
            L_x = L_x*A_t(1).^2;
        end
        L_x = 2*delta*L_x;
        F_old = delta*sum(sum(double(tenmat((Omiga.*(Y-ttm(X,{A{1:Nway}}))).^2,1))))+sum(sum(double(tenmat((X.^2).*D,1))));
        Xt = X;
        for ii = 1 : FISTA_count
            to = t;
            X_old = X;
            % Caculate gradient
            ERR = ttm(Xt,{A{1:Nway}},1:Nway)-Y;
            f_d = 2*delta*ttm(ERR.*Omiga,{A{1:Nway}},'t');
            % Caculte stepsize
            beta_x = (2-sig)/L_x;
            % Update Z
            Z = (Xt-beta_x*f_d)./(2*beta_x*D + 1);
            FZ = delta*sum(sum(double(tenmat((Omiga.*(Y-ttm(Z,{A{1:Nway}}))).^2,1))))+sum(sum(double(tenmat((Z.^2).*D,1))));
            if FZ<F_old;%FX
                X = Z;
                F_old = FZ;
            else
                X = X_old;
            end
            t = (1+sqrt(1+4*to^2))/2;
            Xt = X + to/t*(Z-X)+(to-1)/t*(X-X_old)+to/t*(1-sig)*(Xt-Z);
        end
    % Update factor matrices
        for ii = randperm(Nway)
            ind = sort(setdiff(1:Nway,ii),'ascend');
            Q_temp = double(tenmat(ttm(X,{A{ind}},ind),ii));
            Omiga_m = logical(double(tenmat(Omiga,ii)));
            Y_m = double(tenmat(Y,ii));
            for jj = randperm(TS(ii))
                Q = Q_temp(:,Omiga_m(jj,:));
                Sigma_A = eye(C_S(ii))/(delta*(Q*Q')+eye(C_S(ii)));
                A{ii}(jj,:) = delta*Y_m(jj,Omiga_m(jj,:))*Q'*Sigma_A;
            end
        end
    % Update alpha
        X2 = X.^2;
        for ii = 1 : Nway
            X_m = double(tenmat(X2,ii));
            alpha{ii} = (1 + a)./(sum(X_m,2)+b);
        end
    % Pruning
        for ii = 1 : Nway
            th = C*min(alpha{ii});
            Samp_M = eye(C_S(ii));
            if (size(A{ii},2)>2) && (max(alpha{ii}>th))
                A{ii} = A{ii}(:,alpha{ii}<th);
                Samp_M = Samp_M(alpha{ii}<th,:);
                X = ttm(X,Samp_M,ii);
                C_S(ii) = sum(alpha{ii}<th);
                alpha{ii} = alpha{ii}(alpha{ii}<th);
            end
        end
    % Evaluate the stop criterion
        Y_E = ttm(X,{A{1:Nway}},1:Nway);
        if norm(Y_E-Y_E_old)/sqrt(prod(TS))<1e-5
            break;
        end
    end
    Result.X = X;
    Result.A = A;
    Result.alpha = alpha;
    Result.delta = delta;
    Result.Y = Y_E;
end