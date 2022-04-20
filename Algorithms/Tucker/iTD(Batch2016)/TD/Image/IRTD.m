function Result = IRTD(Y,Omiga,Max_iter,init_size,delta,Y_real)
%--------------------------------------------------------------------------
% Author : Linxiao Yang
% Date : 10/9/2015
% Intro : This is the final version of iterative reweighted Tucker
%         decomposition for incomplete tensor. 
%--------------------------------------------------------------------------
%% Initialization
    TS = size(Y);
    Nway = length(TS);
    C_S = init_size;
    % parameter setting
    a = 0;
    b = 1e-10;
    sig = 0.1;
    FISTA_count = 100;
    C = 1e+9;
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
        iter
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
            F_old = delta*sp_F2(Y-ttm(X,{A{1:Nway}}),Omiga)+sum(sum(double(tenmat((X.^2).*D,1))));
            Xt = X;
            for ii = 1 : FISTA_count
                to = t;
                X_old = X;
                % Caculate gradient
                ERR = sptensor((ttm(Xt,{A{1:Nway}},1:Nway)-Y).*Omiga);
                f_d = 2*delta*ttm(ERR,{A{1:Nway}},'t');
                % Caculte stepsize
                beta_x = (2-sig)/L_x;
                % Update Z
                Z = (Xt-beta_x*f_d)./(2*beta_x*D + 1);
                FZ = delta*sp_F2(Y-ttm(X,{A{1:Nway}}),Omiga)+sum(sum(double(tenmat((Z.^2).*D,1))));
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
                FN = sqrt(diag(A{ii}'*A{ii}));
                X = ttm(X,diag(FN),ii);
                A{ii} = A{ii}*diag(1./FN);
            end
        % Update alpha
            X2 = X.^2;
            for ii = 1 : Nway
                X_m = double(tenmat(X2,ii));
                alpha{ii} = (1 + a)./(sum(X_m,2)+b);
            end
        % Pruning
            for ii = 1 : Nway
                th = C;%*min(alpha{ii});
%                 th = 1e+0;
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
            norm((Y_E-Y_real).*~Omiga)^2/M
            C_S
    end
    Result.X = X;
    Result.A = A;
    Result.alpha = alpha;
    Result.delta = delta;
    Result.Y = Y_E;
end
function R = sp_F2(Y,Omega)
    Y1 = double(tenmat(Y,1));
    Omega1 = double(tenmat(Omega,1));
    y = Y1(:);
    ind = logical(Omega1(:));
    R = sum(y(ind).^2);
end