function [Factor_Es,PER] = WTucker(X_cell,Omega,tucker_rank,OPTS);

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


T          = length(X_cell);
tensor_dim = size(X_cell{1,1});
N          = length(tensor_dim) + 1;

PER = zeros(N,T);
Factor_Es = [];


for t = 1 : T
    t
    Xt    = X_cell{t,1};
    Omega_t = Omega{t,1};   
    
    Max_iter = 100;
    delta    = 0.1;
    Result3  = IRTD(Xt,Omega_t,Max_iter,tucker_rank(1:end-1),delta);
    
    X_re = Result3.Y;
    Factor_Es = Result3.A;

    
   
    
    %% Performance Estimation
    
    if flag == 1
        X_true = Slide_True{t,1};
        
        Err    = X_re - X_true;
        PER(N,t) = norm(Err) / norm(X_true(:));
        for n = 1 : N-1
            PER(n,t) = sub_est_per(Factor_Es{1,n},Factor_True{t,n},'SEP');
        end
        
    else
        
    end
    
end





end 