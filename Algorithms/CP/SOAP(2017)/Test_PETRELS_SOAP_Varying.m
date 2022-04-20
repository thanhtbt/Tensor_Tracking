clc;
clear all;

run_path;

n_exp     = 1;
std_brt   = 1e-3;
sampling  = 0.9;
outlier_den = 0.1;
outlier_fac = 0;
tensor_dims = [10 12 1001];
I = tensor_dims(1);
J = tensor_dims(2);
K = tensor_dims(3);
R = 3;

PER.A = zeros(1,K);
PER.SE = zeros(1,K);
PER.H  = zeros(1,K);
PER.x  = zeros(1,K);
PER.B  = zeros(1,K);



for jj = 1 : n_exp
    fprintf('RUN %d/%d \n ',jj,n_exp)
    
    A = randn(I,R);
    B = randn(J,R);
    C = randn(K,R);
    H = khatrirao(B,A);
    
    X_noisefree = zeros(I,J,K);
    for ii = 1:K
        X_noisefree(:,:,ii) = A * diag(C(ii,:)) * B';
    end
    
    X_true = X_noisefree;
    X = X_true;
    Omega  = rand(I,J,K);
    Omega  = 1 .* (Omega < sampling);
    
    disp('  Online algorithm ')
    
    
    
    
    etaA = 1e-4;
    etaB = 1e-4;
    A_true = A;
    B_true = B;
    
    
    t_train = 1;
   
    At =   randn(I,R); 
    Bt =   randn(J,R); 
    Ct =   randn(t_train,R);
    
    
    %% Main Algorithm
    
    t_start = tic;
    
    Rinv = repmat(100*eye(R),1,I*J);
    OPTS_PETRELS.Rinv = Rinv;
    OPTS_PETRELS.lambda = 0.98;
    Ht = khatrirao(Bt,At);
    for t   =  t_train + 1 : K
        
        if t == 600
            etaA = 1e-2;
            etaB = 1e-2;
        else
            etaA = 1e-4;
            etaB = 1e-4;
        end
        
        A_true = (1-etaA) * A_true + etaA*randn(I,R);
        B_true = (1-etaB) * B_true + etaB*randn(J,R);
        H_true = khatrirao(B_true,A_true);
        Xtrue  = A_true * diag(C(t,:)) * B_true';
        
        Xt  = Xtrue; % + std_brt*randn(I,J);
        % Outliers
        S  = zeros(I,J);
        mS = max(abs(Xt(:)));
        p  = randperm(I*J);
        L  = round(outlier_den*I*J);
        S(p(1:L))  = outlier_fac * mS * rand(L,1);
        Xt = Xt + S;
        
        % Missing
        Xt  = Omega(:,:,t).*Xt;
        
        
%         Omega2  = rand(I,J);
%         Omega2  = 1 .* (Omega2 < 0.8);
%         Xt = Xt.*Omega2;
        
        xt  = vec(Xt);
        Omega_t = Omega(:,:,t);
        
        %% Estimate ct and Ht = khatrirao(Bt,At) using PETRELS
        
        [Ht,ct,OPTS_PETRELS] = Petrels_Subspace(xt,Ht,OPTS_PETRELS);
        
        %% Estimate At, Bt from Ht
        A1 = zeros(I,R);
        B1 = zeros(J,R);
        for r = 1:R
            Hr = reshape(Ht(:,r),I,J);       % ar is the left principal sing vector and conj(cr) the left one
            b1 = Hr'*At(:,r);
            B1(:,r) = conj(b1);
            a1 = Hr*b1;
            A1(:,r) = a1/norm(a1,'fro');     % Normalize it because sing value was included in c1
        end
        At  = A1;
        Bt  = B1;
        %% Re-estimate ct
%         Ht_es    = khatrirao(Bt,At);
%         Ht_omega = Ht_es(idx,:);
%         xt_omega = xt(idx,:);
        %     ct = Ht_omega \ xt_omega;
        Ct = [Ct; ct'];
        
        
        %% Performance Evaluation
        Ht = khatrirao(Bt,At);
        xt_re = Ht * ct;
        
        xt_true = vec(Xtrue);
        er = (xt_re - xt_true);
        [~,~,~,erA]=solve_perm_scale(At,A_true);
        [~,~,~,erB]=solve_perm_scale(Bt,B_true);
        
        erC= sqrt( 1 - ( (C(t,:)*ct)/( norm(ct)*norm(C(t,:)) ) )^2);
        
        PER.A(1,t)  = PER.A(1,t)  + erA;
        PER.SE(1,t) = PER.SE(1,t) + sub_est_per(At,A_true,'Angle');
        PER.H(1,t)  = PER.H(1,t)  + sub_est_per(Ht,H_true,'Angle');
        PER.B(1,t)  = PER.B(1,t)  + erB;
        PER.x(1,t)  = PER.x(1,t)  + norm(er);
        
    end
    toc(t_start);
    
    fprintf('Online Tracking:   SE(A,At)     = %d | RE(A,At)     = %d \n',sub_est_per(At,A,'Angle'),erA);
    
    
    
end

RE_A = (PER.A)/n_exp;
SE_A = (PER.SE)/n_exp;
SE_H = (PER.H)/n_exp;
RE_x = (PER.x)/n_exp;
RE_B = (PER.B)/n_exp;

%% PLOT
makerSize = 14;
numbMarkers = 50;
LineWidth = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];

blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];

fig = figure;
k = 2;
d1 = semilogy(1:k:K,SE_A(1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;
d2 = semilogy(1:k:K,SE_H(1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;

lgd = legend([d1 d2],'Loading Factor $\mathbf{A}_t$','$\mathbf{H}_t = \mathbf{B}_t \odot \mathbf{A}_t $');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin(\theta(\mathbf{U}_{est}, \mathbf{U}_{true} ))$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 22);
axis([0 K 1e-4 1e1]);
grid on;


fig2 = figure;
k = 10;
d1 = semilogy(1:k:K,RE_A(1:k:end),'linestyle','-','color',red_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;

d2 = semilogy(1:k:K,RE_B(1:k:end),'linestyle','-','color',blue_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;

d3 = semilogy(1:k:K,RE_x(1:k:end),'linestyle','-','color',gree_n,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;

lgd = legend([d1 d2 d3],'$\mathbf{A}_t$','$\mathbf{B}_t$','$\mathbf{X}_t$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$ \displaystyle \frac{\| \mathbf{U}_{est} - \mathbf{U}_{true}  \|}{\| \mathbf{U}_{true}\|} $','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig2, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 22);
axis([0 K 1e-4 1e1]);
grid on;





