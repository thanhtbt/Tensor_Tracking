%% This demo illustrates the tracking ability of ACP in time-varying environments

clc; clear all; close all;
run_path;

n_exp      = 3;                      % Number of independent runs
fac_noise  = 1e-4;                   % Noise factor
sampling   = 0.9;                    % (1-sampling) = missing density 
tensor_dim = [20 20 20 1000];        % Tensor Dimension
R          = 5;                      % CP rank
Mag        = [1e-1 1e-2 1e-3 1e-4];  % time-varying factor

PER_ACP   = cell(length(Mag),1);
for m = 1 : length(Mag)
    PER_ACP{m,1} = zeros(length(tensor_dim),tensor_dim(end));
end

%% Main Program
N = length(tensor_dim);
for m = 1 : length(Mag);
    mag = Mag(m);
    fprintf('\n Time-varying factor %f: \n',mag);
    epsilon      = mag*ones(1,tensor_dim(end));
    epsilon(600) = 1;    
    PER_Our   = zeros(length(tensor_dim),tensor_dim(end));
    for ii = 1 : n_exp
        fprintf(' + Run [%d/%d]: \n', ii, n_exp)
        
        %% Generate True Tensor
        Factor = cell(1,N);
        for n = 1 : N
            Factor{1,n} = randn(tensor_dim(n),R);
        end
        
        % True Data
        Size_Fixed = tensor_dim(1:end-1);
        Num_Slides = tensor_dim(end);
        X_stream   = cp_online_generator(Size_Fixed,R,Num_Slides,epsilon);
        
        OPTS.Factor_True = X_stream(:,2);
        OPTS.Slide_True  = X_stream(:,1);
      
        
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
          
        %% Our Method
        t_start = tic;
        [~,PER] = ACP(X_cell,Omega,R,OPTS);
        t_end   = toc(t_start);
        PER_Our = PER_Our + PER;
        
    end
    PER_Our      = PER_Our/n_exp;
    PER_ACP{m,1} =  PER_Our;

    
end

% For four values of time-varying factor
PER_Our1 =  PER_ACP{1,1};
PER_Our2 =  PER_ACP{2,1};
PER_Our3 =  PER_ACP{3,1};
PER_Our4 =  PER_ACP{4,1};


%% PLOT RESULTS
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

%%
fig = figure;
hold on;
k = 1;
K = tensor_dim(end);

d1 = semilogy(1:k:K,PER_Our1(3,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_Our1(3,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_Our1(3,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



d2 = semilogy(1:k:K,PER_Our2(3,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_Our2(3,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER_Our2(3,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_Our3(3,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,PER_Our3(3,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,PER_Our3(3,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);



d4 = semilogy(1:k:K,PER_Our4(3,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d41 = plot(1:100:K,PER_Our4(3,1:100:end),...
    'marker','p','markersize',makerSize,...
    'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_Our4(3,1:1),...
    'marker','p','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);



lgd = legend([d12 d22 d32  d42],'$\epsilon = 10^{-1}$','$\epsilon = 10^{-2}$','$\epsilon = 10^{-3}$','$\epsilon = 10^{-4}$');
lgd.FontSize = 22;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathcal{X}_{tr}, \mathcal{X}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(gca, 'YScale', 'log')
set(h,'Xtick',0:200:K,'FontSize',24,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
axis([0 K 5e-5 1e1]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);

%% 

fig2 = figure;
hold on;
k = 1;
K = tensor_dim(end);

d1 = semilogy(1:k:K,PER_Our1(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_Our1(1,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_Our1(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d2 = semilogy(1:k:K,PER_Our2(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_Our2(1,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER_Our2(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_Our3(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d31 = plot(1:100:K,PER_Our3(1,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,PER_Our3(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d4 = semilogy(1:k:K,PER_Our4(1,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
% 'marker','d','markersize',makerSize,'markerfacecolor','w',...
d41 = plot(1:100:K,PER_Our4(1,1:100:end),...
    'marker','p','markersize',makerSize,...
    'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_Our4(1,1:1),...
    'marker','p','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);



lgd = legend([d12 d22 d32 d42],'$\epsilon = 10^{-1}$','$\epsilon = 10^{-2}$','$\epsilon = 10^{-3}$','$\epsilon = 10^{-4}$');
lgd.FontSize = 22;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig2, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(gca, 'YScale', 'log')
set(h,'Xtick',0:200:K,'FontSize',24,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');

axis([0 K 5e-5 1e1]);
grid on;
box on;
set(fig2, 'units', 'inches', 'position', [0.5 0.5 8 7]);



