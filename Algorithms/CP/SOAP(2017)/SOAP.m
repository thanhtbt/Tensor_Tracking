function [PER,Factor] = SOAP(X,R,OPTS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by NGUYEN Viet-Dung, Karim Abed-Meraim and Nguyen Linh-Trung
% modified from programs by D. Nion and N.D. Sidiropoulos
%
% If you find our program useful, please cite to: 
% Nguyen Viet-Dung, Karim Abed-Meraim and Nguyen Linh-Trung
% "Second-order optimization based adaptive PARAFAC decomposition of
% three-way tensors"
% Digital Signal processing, 63:100-111, April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin < 3
    flag = 0; % without performance estimation part
    Factors = [];
    PER     = [];
else
    flag = 1; % with performance estimation part
    if isfield(OPTS,'Factor_True')
        Factors  = OPTS.Factor_True;
        flag_factor = 1;
    else
        flag_factor = 0;
    end
    if isfield(OPTS,'Slide_True')
        Slide_true = OPTS.Slide_True;
        flag_slide = 1;
    else
        flag_slide = 0;
    end
end

[I,J,T] = size(X);
%% Performance Evaluation
PER.A  = zeros(1,T);
PER.B  = zeros(1,T);
PER.X  = zeros(1,T);

%% Training
t_train = max(I,J);
X_train = X(:,:,1:t_train);

[A1, B1 ,C1]=cp3_alsls(X_train,R);
At = A1;
Bt = B1;
Ct = C1;
H4 = kat_rao(Bt,At);
H4inv = pinv(H4);
Rinv4 = pinv(Ct'*Ct);  

for t   =  1 : T
    Xt      = X(:,:,t);
    x       = Xt(:);
    %% Estimate ct 
    [ct,At,Bt,H4,H4inv,Rinv4] = complexFastSOAParafac(x,At,Bt,H4,H4inv,Rinv4,I,J,R,t);
    Ct = [Ct; ct];        
    
    
    %% Performance Evaluation
    
     if flag == 1  
        % Data
        if flag_slide == 1
            Xt_true    = Slide_true{t,1};
        else
            Xt_true = Xt;
        end
            Ht      = khatrirao(Bt,At);
            xt_re   = Ht * ct';
            xt_true = (Xt_true(:));
            er      = (xt_re - xt_true);
            PER.X(1,t)  = PER.X(1,t)  + norm(er)/norm(xt_true);
       
        % Factors 
        if flag_factor == 1
            Factors_t   = Factors{t,1};
            A_true      = Factors_t{1,1};
            B_true      = Factors_t{1,2};
            [~,~,~,erA] = solve_perm_scale(At,A_true);
            [~,~,~,erB] = solve_perm_scale(Bt,B_true);
            PER.A(1,t)  = PER.A(1,t)  + erA;
            PER.B(1,t)  = PER.B(1,t)  + erB;
        else
        end
     end
    
     

end
Factor.A = At;
Factor.B = Bt;
Factor.C = Ct; 
     

end
