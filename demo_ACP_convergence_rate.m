%% This demo illustrates the convergence rate of ACP

clc; clear all; close all;
run_path;

n_exp         = 3;                      % Number of independent runs
fac_noise     = 1e-3;                   % Noise factor
tensor_dim    = [20 20 20 1000];        % Dimension of streaming tensors
R             = 2;                      % CP Rank
Sampling      = [0.9 0.7 0.5];          % Sampling rate >< Missing density
time_varying  = 0;                      % Time-varying factor
epsilon    = time_varying*ones(1,tensor_dim(end));

PER_Our   = zeros(length(tensor_dim),tensor_dim(end));
VAR_U_Our = zeros(length(tensor_dim),tensor_dim(end));
VAR_f_Our = zeros(1,tensor_dim(end));

RESULT = cell(length(Sampling),2);
%% Main Program
for m =  1 : length(Sampling),
    sampling   = Sampling(m);
    fprintf('\n Sampling %f \n', sampling)
    for ii = 1: n_exp
        
        fprintf('\n + Run [%d/%d]: ', ii, n_exp)
        %% Generate True Tensor
        N = length(tensor_dim);
        Factor = cell(1,N);
        for n = 1 : N
            Factor{1,n} = randn(tensor_dim(n),R);
        end
        
        % True Data
        Size_Fixed = tensor_dim(1:end-1);
        Num_Slides = tensor_dim(end);
        X_stream   = cp_online_generator_rotation(Size_Fixed,R,Num_Slides,epsilon);
        
        OPTS.Factor_True = X_stream(:,2);
        OPTS.Slide_True  = X_stream(:,1);
        OPTS.lambda      = 0.5;
        
        % Add Missing + Gaussian Noise
        X_cell = cell(tensor_dim(end),1);
        Omega  = cell(tensor_dim(end),1);
        for t  = 1 : tensor_dim(end)
            Omega_t = rand(tensor_dim(1:end-1));
            Omega_t = 1 .* (Omega_t <= sampling);
            Omega{t,1}  = Omega_t;
            X_cell{t,1} = OPTS.Slide_True{t,1}  + fac_noise * randn(tensor_dim(1:end-1));
            X_cell{t,1} = Omega_t .* X_cell{t,1};
        end
        
        t_start = tic;
        [~,PER,VAR_U,VAR_f] = ACP_C(X_cell,Omega,R,OPTS);
        toc(t_start);
        PER_Our   = PER_Our + PER;
        VAR_U_Our = VAR_U_Our +  VAR_U;
        VAR_f_Our = VAR_f_Our +  VAR_f;

    end
    PER_Our   = PER_Our/n_exp;
    VAR_U_Our = VAR_U_Our/n_exp;
    VAR_f_Our = VAR_f_Our/n_exp;

    RESULT{m,1} = PER_Our;
    RESULT{m,2} = VAR_U_Our;
    RESULT{m,3} = VAR_f_Our;

end





%% PLOT RESULTS
makerSize    = 14;
numbMarkers  = 50;
LineWidth    = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = 'g';%[0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];

blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];

%%
fig = figure;
k = 10;
K = tensor_dim(end);
KL = 200;

hold on;

d2 = semilogy(1:k:K,RESULT{2,2}(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,RESULT{2,2}(1,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,RESULT{2,2}(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,RESULT{3,2}(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,RESULT{3,2}(1,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,RESULT{3,2}(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d1 = semilogy(1:k:K,RESULT{1,2}(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,RESULT{1,2}(1,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,RESULT{1,2}(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



lgd = legend([d12 d22 d32 ],'$\rho = 10\%$','$\rho = 30\%$','$\rho = 50\%$');
lgd.FontSize = 16;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\| \mathbf{U}_{t+1} - \mathbf{U}_{t} \|_F$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;


%%
fig2 = figure;
k = 10;
K = tensor_dim(end);
KL = 200;

hold on;

d2 = semilogy(1:k:K,RESULT{2,3}(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,RESULT{2,3}(1,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,RESULT{2,3}(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,RESULT{3,3}(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,RESULT{3,3}(1,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,RESULT{3,3}(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d1 = semilogy(1:k:K,RESULT{1,3}(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,RESULT{1,3}(1,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,RESULT{1,3}(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



lgd = legend([d12 d22 d32],'$\rho = 10\%$','$\rho = 30\%$','$\rho = 50\%$');
lgd.FontSize = 16;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('Objective Value $f(\mathbf{U}_t)$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;


