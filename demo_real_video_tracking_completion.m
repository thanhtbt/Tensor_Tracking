%% This demo shows the effectiveness of ACP and ATD for online video completion

clear;clc
run_path;

%% Data

load X_Lobby.mat
%% Missing
sampling = 0.9;
tensor_dim = size(X_video);

Omega  = rand(tensor_dim);
Omega  = 1 .* (Omega <= sampling);
Noise  = randn(tensor_dim);
fac    = 0; 
X_miss = Omega .* (X_video + fac*Noise);

%% Inputs
tucker_rank = [10 10 10]; % [15 15 15]; %
R_CP        = 10;
OPTS.lambda = 0.9;
OPTS.Xtrue  = X_video;


disp('Processing ....')
%% Algorithms
%% Ours
    t_start = tic;
    [PER_ATD,Xre_ATD] = ATD_Video(X_miss,Omega,tucker_rank,OPTS);
    t_end = toc(t_start);
    fprintf('+ ATD (2021): %f(s)\n',t_end)

%% Ours
    t_start = tic;
    [PER_ACP,Xre_ACP] = ACP_Video(X_miss,Omega,R_CP,OPTS);
    t_end = toc(t_start);
    fprintf('+ ACP (2021): %f(s)\n',t_end)
    
    
%% TESGD (2015)
    maxepochs               = 1;
    verbose                 = 2;
    tolcost                 = 1e-8;
    permute_on              = false;
    clear options;
    options.maxepochs       = maxepochs;
    options.tolcost         = tolcost;
    options.permute_on      = permute_on;
    options.stepsize        = 0.1;
    options.lambda          = 0.001;  % Forgetting paramter
    options.mu              = 0.01;   % Regualization paramter
    options.tw_flag         = 0;     % 0:Exponential Window, 1:Truncated Window (TW)
    options.tw_len          = 10;    % Window length for Truncated Window (TW) algorithm
    options.store_subinfo   = true;
    options.store_matrix    = false;
    options.verbose         = verbose;
   
    t_start = tic;
    [PER_TeCPSGD,Xre_TeCPSGD] = TeCPSGD_Video(X_miss,Omega,R_CP,options,OPTS);
    t_end = toc(t_start);   
    fprintf('+ TeSGD (2015): %f(s)\n',t_end)
    
%% OLSTEC (2019)    
    clear options;
    options.maxepochs       = maxepochs;
    options.tolcost         = tolcost;
    options.permute_on      = permute_on;    
    options.lambda          = 0.7;  % Forgetting paramter
    options.mu              = 0.1;  % Regualization paramter
    options.tw_flag         = 0;    % 0:Exponential Window, 1:Truncated Window (TW)
    options.tw_len          = 10;   % Window length for Truncated Window (TW) algorithm
    options.store_subinfo   = true;     
    options.store_matrix    = false; 
    options.verbose         = verbose; 

    t_start = tic;
    [PER_OLSTEC,Xre_OLSTEC] = OLSTEC_Video(X_miss,Omega,R_CP,options,OPTS);
    t_end = toc(t_start);
    fprintf('+ OLSTEC (2019): %f(s)\n',t_end)

    
%% CP-PETRELS   
    t_start = tic;
    [PER_PETRELS,Xre_PETRELS] = PETRELS_SOAP_Video(X_miss,Omega,R_CP,OPTS);
    t_end = toc(t_start);
    fprintf('+ CP-PETRELS: %f(s)\n',t_end)

%% Test
t  = 500;
X_t          = X_video(:,:,t);
Xmiss_t      = X_miss(:,:,t);
Xret_ATD     = Xre_ATD(:,:,t);
Xret_ACP     = Xre_ACP(:,:,t);
Xret_TeCPSGD = Xre_TeCPSGD(:,:,t);
Xret_OLSTEC  = Xre_OLSTEC(:,:,t);
Xret_PETRELS  = Xre_PETRELS(:,:,t);


% 
figure; 
subplot(231); imagesc(Xmiss_t);axis off; 
subplot(232); imagesc(Xret_ACP);axis off 
subplot(233); imagesc(Xret_ATD);axis off
subplot(234); imagesc(Xret_TeCPSGD);axis off
subplot(235); imagesc(Xret_OLSTEC);axis off
subplot(236); imagesc(Xret_PETRELS);axis off

colormap('gray');   

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

%%

K = tensor_dim(end);

fig = figure;
hold on;
k = 5;
step = 200;

d2 = semilogy(1:k:K,PER_TeCPSGD(1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d21 = plot(1:step:K,PER_TeCPSGD(1:step:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d22 = semilogy(1:1,PER_TeCPSGD(1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_PETRELS(1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d31 = plot(1:step:K,PER_PETRELS(1:step:end),...
    'marker','^','markersize',makerSize,...
    'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d32 = semilogy(1:1,PER_PETRELS(1:1),...
    'marker','^','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);


d4 = semilogy(1:k:K,PER_OLSTEC(1:k:end),...
    'linestyle','-','color',blue_n,'LineWidth',LineWidth);
d41 = plot(1:step:K,PER_OLSTEC(1:step:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',blue_n,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_OLSTEC(1:1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',blue_n,'LineWidth',LineWidth);


d0 = semilogy(1:k:K,PER_ACP(1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d01 = plot(1:step:K,PER_ACP(1:step:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d02 = semilogy(1:1,PER_ACP(1:1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);

d1 = semilogy(1:k:K,PER_ATD(1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:step:K,PER_ATD(1:step:end),...
    'marker','p','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_ATD(1:1),...
    'marker','p','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



lgd = legend([d22 d42 d32 d02 d12],'TeCPSGD','OLSTEC','CP-PETRELS','ACP (Proposed)','ATD (Proposed)');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Frame Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathbf{X}_t, \mathbf{X}_{true})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
h=gca;
set(gca, 'YScale', 'log')
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:300:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 22);
grid on;
box on;
axis([1 K 5e-2 2])
%
