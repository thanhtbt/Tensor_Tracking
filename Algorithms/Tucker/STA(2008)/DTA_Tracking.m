function [Factor_Es,PER] = Tracking_STA(X_cell,tucker_rank,OPTS)

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

%% 
L   = length(X_cell);
tensor_dim = size(X_cell{1,1});
N   = length(tensor_dim) + 1;
PER = zeros(N,L);
Factor_Es = cell(1,N);

[T,S] = DTA(X_cell{1,1},tucker_rank(1:end-1));
alpha = 1;
for t = 2 : L
    
    Xt      = X_cell{t,1};
    [T,S]   = DTA(Xt,tucker_rank(1:end-1),S,alpha);
   
    %% Performance Estimation
    Tensor_Core  = T.core;
    for n = 1 : N-1
        Factor_Es{1,n} = T.U{n};
    end
  
    X_re   = ttm(Tensor_Core,Factor_Es(1,1:N-1));
    X_true = Slide_True{t,1};
    ER = X_re - X_true;
    err = norm(ER)/norm(X_true);
    PER(N,t) = err;
    for n = 1 : N - 1
        PER(n,t) = sub_est_per(Factor_Es{1,n},Factor_True{t,n},'SEP');
    end   
    
    
end
Factor_Es{1,N} = [];
end