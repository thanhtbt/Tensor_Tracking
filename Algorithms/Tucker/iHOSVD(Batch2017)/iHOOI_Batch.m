function [Factor_Es,PER] = iHOOI_Batch(X_cell,tucker_rank,OPTS)

%% Missing

flag = 1; % with performance estimation part
if isfield(OPTS,'Factor_True')
    Factor_True  = OPTS.Factor_True;
    flag_factor = 1;
else
    flag_factor = 0;
end
if isfield(OPTS,'Slide_True')
    Slide_True = OPTS.Slide_True;
    flag_slide = 1;
else
    flag_slide = 0;
end




T   = length(X_cell);
tensor_dim = size(X_cell{1,1});
N   = length(tensor_dim) + 1;

PER = zeros(N+1,T);
Factor_Es = [];

for t = 1 : T
    t_start  = tic;
    Xt       = X_cell{t,1};
    
    known = find(Xt.data);
    data  = Xt(known);
    
    Nway  = tensor_dim;
    EstCoreNway = [1 1 1];
    opts = [];
    opts.maxit = 100; opts.tol = 1e-4;
    opts.rank_adj = ones(1,N-1);
    opts.rank_max = 1*tucker_rank(1:N-1);
    
    [A,C,Out] = iHOOI(data,known,Nway,EstCoreNway,opts);
    
    t_end = toc(t_start);
    %% Performance Estimation
    
    if flag == 1
        X_true = Slide_True{t,1};
        X_re   = full(ttensor(C,A));
        Err    = X_re - X_true;
        PER(N,t) = norm(Err) / norm(X_true(:));
        for n = 1 : N-1
           PER(n,t) = sub_est_per(A{1,n},Factor_True{t,n},'SEP');
        end
        PER(N+1,t) = t_end;
    else
        
    end
    
end
% Factor_Es{1,N} = last_fac;
end