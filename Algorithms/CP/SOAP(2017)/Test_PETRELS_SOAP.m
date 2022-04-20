clc;
clear all; close all;

n_exp     = 10;
std_brt   = 1e-3;
sampling  = 0.9;

tensor_dims = [10 12 1000];
I = tensor_dims(1);
J = tensor_dims(2);
K = tensor_dims(3);
Rank = 3; R = Rank;

PER_SOAP.A  = zeros(1,K);
PER_SOAP.X  = zeros(1,K);
PER_SOAP.B  = zeros(1,K);
PER_SOAP.time = 0;




%% Algorithms
for jj = 1 : n_exp
    fprintf('RUN %d/%d \n',jj,n_exp)
    
    A_true = randn(I,Rank);
    B_true = randn(J,Rank);
    C_true = randn(K,Rank);
    
    X_true = zeros(I,J,K);
    for ii = 1:K
        X_true(:,:,ii) = A_true * diag(C_true(ii,:)) * B_true';
    end
    X = X_true + std_brt*randn(I,J,K);
    
    Omega  = rand(I,J,K);
    Omega  = 1 .* (Omega < sampling);
    X = X.*Omega;
    
    %% PETRELS-SOAP Algorithm
    OPTS.A = A_true;
    OPTS.B = B_true;
    OPTS.X = X_true;
    t_start = tic;
    [PER,~] = PETRELS_SOAP(X,Omega,Rank,OPTS);
    t_end = toc(t_start);
    PER_SOAP.time = PER_SOAP.time + t_end; 
    fprintf('+ PETRELS-SOAP: %f(s) \n',t_end)

    PER_SOAP.A = PER.A + PER_SOAP.A;
    PER_SOAP.B = PER.B + PER_SOAP.B;
    PER_SOAP.X = PER.X + PER_SOAP.X;
    
    
    
end
PER_SOAP.A = PER_SOAP.A / n_exp;
PER_SOAP.B = PER_SOAP.B / n_exp;
PER_SOAP.X = PER_SOAP.X / n_exp;
PER_SOAP.time = PER_SOAP.time / n_exp;



%% PLOT
makerSize = 14;
numbMarkers = 50;
LineWidth = 1.5;
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
k = 1;
d1 = semilogy(1:k:K,PER_SOAP.A(1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;
d2 = semilogy(1:k:K,PER_SOAP.B(1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
hold on;
d3 = semilogy(1:k:K,PER_SOAP.X(1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...

lgd = legend([d1 d2 d3],'$\mathbf{A}$','$\mathbf{B}$','$\mathbf{X}$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathbf{U}_t, \mathbf{U})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 22);
axis([0 K 1e-5 1e1]);
grid on;


