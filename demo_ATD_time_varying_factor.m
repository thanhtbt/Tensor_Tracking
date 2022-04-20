%% This demo illustrates the tracking ability of ATD in time-varying environments

clear;clc; close all;
run_path;


n_exp       = 3;                      % Number of independent runs
tensor_dim  = [20 20 20 500];         % Dimension of streaming tensors
tucker_rank = [3  3  3  3];           % Tucker rank
fac_noise   = 1e-3;                   % Noise level
sampling    = 0.9;                    % 1 - Sampling = Missing density
MAG         = [1e-4 1e-3 1e-2 1e-1];  % Time-varying factor
N           = length(tensor_dim);
T           = tensor_dim(end);

%% Main Program
PER_U = cell(length(Lambda),length(MAG));
PER_X = cell(length(Lambda),length(MAG));

for jj = 1 : length(MAG)
    PER_U{1,jj} = zeros(N-1,T);
    PER_X{1,jj} = zeros(1,T);
end

for jj = 1 : length(MAG)
    mag = MAG(jj);
    epsilon        = mag*ones(1,tensor_dim(end));
    epsilon(300)   = 1; % Create an abrupt change at t = 300';
    fprintf('\n \n // Time-varying factor = %f \n',mag)
    PER_U_n = zeros(N-1,T);
    PER_X_n = zeros(1,T);
    
    for n = 1 : n_exp
        fprintf('\n  +  Run [%d/%d]: ',n,n_exp);
        %% Generate Tensor
        [Xtrue_cell,Factor_True] = tucker_generate_tensor(tensor_dim,tucker_rank,0,epsilon);
        OPTS.Factor_True = Factor_True;
        OPTS.Slide_True  = Xtrue_cell;
        
        % Missing + Gaussian Noise
        X_cell = cell(tensor_dim(end),1);
        Omega  = cell(tensor_dim(end),1);
        for t  = 1 : tensor_dim(end)
            Omega_t     = rand(tensor_dim(1:end-1));
            Omega_t     = 1 .* (Omega_t <= sampling);
            Omega{t,1}  = Omega_t;
            X_cell{t,1} = Xtrue_cell{t,1} + fac_noise*randn(tensor_dim(1:end-1));
            X_cell{t,1} = Omega_t .* X_cell{t,1};
        end
        
        %% Algorithm
        t_start = tic;
        [~,PER] = ATD(X_cell,Omega,tucker_rank,OPTS);
        toc(t_start);
        PER_U_n = PER_U_n + PER(1:N-1,1:end);
        PER_X_n = PER_X_n + PER(N,1:end);
    end
    PER_U{1,jj} = PER_U_n/n_exp;
    PER_X{1,jj} = PER_X_n/n_exp;
end

%% PLOT RESULTS

PER1 = PER_U{1,1};
PER2 = PER_U{1,2};
PER3 = PER_U{1,3};
PER4 = PER_U{1,4};

makerSize = 14;
numbMarkers = 50;
LineWidth = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o = [0, 0.5, 0];
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


subplot(212);

k = 1;
K = tensor_dim(end);
hold on;

d1 = semilogy(1:k:K,PER1(end,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER1(end,1:100:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER1(end,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d2 = semilogy(1:k:K,PER2(end,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER2(end,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER2(end,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER3(end,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,PER3(end,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d32 = semilogy(1:1,PER3(end,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d4 = semilogy(1:k:K,PER4(end,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d41 = plot(1:100:K,PER4(end,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER4(end,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);


ylabel('SEP $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(gca, 'YScale', 'log');
set(gca,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(gca,'Xtick',0:100:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(gca,'FontSize', 24);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);


subplot(211);

hold on;

d1 = semilogy(1:k:K,PER_X{1,1}(end,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_X{1,1}(end,1:100:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_X{1,1}(end,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d2 = semilogy(1:k:K,PER_X{1,2}(end,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_X{1,2}(end,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER_X{1,2}(end,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_X{1,3}(end,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,PER_X{1,3}(end,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d32 = semilogy(1:1,PER_X{1,3}(end,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d4 = semilogy(1:k:K,PER_X{1,4}(end,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d41 = plot(1:100:K,PER_X{1,4}(end,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_X{1,4}(end,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

lgd = legend([d12 d22 d32 d42],'$\epsilon = 10^{-4}$','$\epsilon = 10^{-3}$','$\epsilon = 10^{-2}$', ...
    '$\epsilon = 10^{-1}$');
lgd.FontSize = 22;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
ylabel('RE $(\mathcal{X}_{tr}, \mathcal{X}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:100:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 24);
%axis([0 K 5e-4 1e0]);
grid on;
box on;







