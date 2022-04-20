%% This demo illustrates the effect of missing data on performance of ACP

clc; clear all; close all;
run_path;

n_exp         = 2;                        % Number of independent runs
fac_noise     = 1e-3;                     % Noise factor
SAMPLING      = [0.9 0.7 0.5];            % (1-sampling) = missing density 
time_varying  = [1e-4 1e-3 1e-2];         % Time-varying factor
tensor_dim    = [20 20 20 1000];          % Tensor Dimension
R             = 2;                        % CP rank

%% Main Program
PER_ACP = cell(length(SAMPLING),length(time_varying));
index = 1;
for m = 1 : length(time_varying)
    mag        = time_varying(m);
    epsilon    = mag*ones(1,tensor_dim(end));
    epsilon(600)= 1; % Aim to create an abrupt change
    
    for s = 1:length(SAMPLING),
        sampling   = SAMPLING(s);
        PER_Our = zeros(length(tensor_dim),tensor_dim(end));
        fprintf('\n \n // Case %d : [time-varying = %f | sampling = %f] \n',index, mag,sampling)
        index = index + 1;
        for ii = 1: n_exp
            
            fprintf('\n + Run [%d/%d]: ', ii,n_exp)
            %% Generate True Tensor
            N = length(tensor_dim);
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
            [~,PER] = ACP(X_cell,Omega,R,OPTS);
            toc(t_start);
            PER_Our = PER_Our + PER;
            
        end
        PER_ACP{s,m} = PER_Our/n_exp;
        
    end
end

X11 = PER_ACP{1,1};
X12 = PER_ACP{1,2};
X21 = PER_ACP{2,1};
X22 = PER_ACP{2,2};
X31 = PER_ACP{3,1};
X32 = PER_ACP{3,2};

%% PLOT RESULTS
makerSize   = 14;
numbMarkers = 50;
LineWidth   = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];
blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(end,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350 0.580 0.2840];

%%
fig = figure;
k = 1;
K = tensor_dim(end);

%%
subplot(211);
hold on;

d1 = semilogy(1:k:K,X11(end,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:100:K,X11(end,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,X11(end,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d2 = semilogy(1:k:K,X21(end,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,X21(end,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,X21(end,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,X31(end,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d31 = plot(1:100:K,X31(end,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d32 = semilogy(1:1,X31(end,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d4 = semilogy(1:k:K,X12(end,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d41 = plot(1:100:K,X12(end,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,X12(end,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d5 = semilogy(1:k:K,X22(end,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d51 = plot(1:100:K,X22(end,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d52 = semilogy(1:1,X22(end,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d6 = semilogy(1:k:K,X32(end,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d61 = plot(1:100:K,X32(end,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d62 = semilogy(1:1,X32(end,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);



lgd = legend([d12 d22 d32],'$\rho = 10\%$','$\rho = 30\%$','$\rho = 50\%$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

ylabel('RE $(\mathcal{X}_{tr}, \mathcal{X}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;

%% 
subplot(212);
hold on;
d1 = semilogy(1:k:K,X11(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,X11(1,1:100:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,X11(1,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d2 = semilogy(1:k:K,X21(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,X21(1,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,X21(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,X31(1,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d31 = plot(1:100:K,X31(1,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d32 = semilogy(1:1,X31(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d4 = semilogy(1:k:K,X12(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d41 = plot(1:100:K,X12(1,1:100:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,X12(1,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


d5 = semilogy(1:k:K,X22(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d51 = plot(1:100:K,X22(1,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d52 = semilogy(1:1,X22(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d6 = semilogy(1:k:K,X32(1,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d61 = plot(1:100:K,X32(1,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d62 = semilogy(1:1,X32(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);


