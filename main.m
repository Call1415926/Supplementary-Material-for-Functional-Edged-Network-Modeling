%% Functional-Edged Network Modeling: end-to-end simulation driver
% This is the only script users need to run. Backend functions are organized
% under src/, and Tensor Toolbox is bundled under external/.

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root,'src')));
addpath(fullfile(project_root,'external','tensor_toolbox'));

results_dir = fullfile(project_root,'results');
if ~exist(results_dir,'dir')
    mkdir(results_dir);
end


warning_state = warning;
myoutput = struct([]);

%% User-adjustable simulation and optimization settings
% The current reduced-size setting is convenient for checking installation.
% For the large-scale setting underlying Figure 5, use m=50, s1=15, K=25.
rep = 20; % number of Monte Carlo replications
m = 50; % number of nodes
n = m;
s1 = 15; % node-mode Tucker rank
s2 = s1;
K = 25; % functional-mode Tucker rank
sigma = sqrt(0.2); % standard deviation of observation noise
w = 0.4; % missing rate
N = 20; % upper scale parameter used to generate the core tensor
t_start = -0.98;
t_end = 1;
seq = 0.02;
r = [s1, s1, K];
max_m = 5; 
tol = 1e-6; %tolerance
maxiter = 1000;
seqtime = t_start:seq:t_end;
N_min = round(length(seqtime)*(1-w));
N_max = round(length(seqtime)*(1-w));
mydim = [n, m, length(seqtime)]; 
Sigma = dif_mat(mydim(3)); 
alpha = 0.1 * ones(1, r(3)); % smoothness penalty
Alpha = diag(alpha); 

for tt = 1:rep
    % Generate the latent functional adjacency tensor and its true factors.
    [X_true,U_true,G_true] = generate_data_true(...
        n,m,s1,s2,K,N,t_start,t_end,seq);
    myoutput(tt).G_true = G_true;
    myoutput(tt).U_true = U_true;
    
    % Add Gaussian noise and independently mask functional observations.
    [X_observe,Omiga,seqtime,X_observe_full,X_observe_full_nornd] = ...
        generate_data_observe(X_true,N_min,N_max,t_start,t_end,seq,sigma);
    myoutput(tt).Omiga = Omiga;
    myoutput(tt).X_observe_full = X_observe_full;
    myoutput(tt).X_observe_full_nornd = X_observe_full_nornd;
    
    X_observe_central = X_observe;
    myoutput(tt).X_observe_central = X_observe_central;
    
    [G_old, U_old, X_old] = Initialization_tucker(X_observe_central, r);

    % Initialize and optimize the shared-Phi functional Tucker model.
    temp = grad_Euc(G_old, U_old, X_old, X_observe_central, Omiga, Sigma, Alpha); 
    epsilon_old = Euc_to_Rie(temp, G_old, U_old);
    eta_old = -epsilon_old;
    alphA_old = find_alpha(X_old,eta_old,Omiga,X_observe_central);
    m_old = find_m(X_old,alphA_old,U_old,eta_old,epsilon_old,r,...
        X_observe_central,Omiga,Sigma,alpha,max_m);
    [G_new, U_new,X_new] = Retraction(X_old,2^(-m_old)*alphA_old*eta_old,r);
    score_new = target_f(X_new, X_observe_central,Omiga, U_new{3},Sigma, alpha);
    [G_next,U_next,X_next,score_next,score_list] = conjugate_dradient(G_old, U_old, X_old,eta_old,...
        G_new, U_new, X_new,score_new,X_observe_central, Omiga, Sigma, Alpha,alpha,max_m,r, tol, maxiter);

    % Factor alignment is used only for reporting U and G. The tensor used
    % for prediction and evaluation remains the optimizer output X_next.
    [~, G_output, U_output] = adjust(X_next, G_next, U_next, X_true, G_true, U_true,t_start, t_end, seq);
    X_output = X_next;
    [MSE,R2,Rec_esti, Rec_true,Cyc_esti, Cyc_true,InTwo_esti, ... 
        InTwo_true, OutTwo_esti,OutTwo_true,Densi_esti,Densi_true] = evlua(X_output, X_observe, X_observe_full_nornd,Omiga);
    
    score_list = [score_list, target_f(X_output, X_observe_central,Omiga, U_next{3},Sigma, alpha)];
    
    myoutput(tt).G_output = G_output;
    myoutput(tt).U_output = U_output;
    myoutput(tt).X_output = X_output;
    myoutput(tt).MSE = MSE;
    myoutput(tt).R2 = R2;
    myoutput(tt).Rec_esti = Rec_esti;
    myoutput(tt).Rec_true = Rec_true;
    myoutput(tt).Cyc_esti = Cyc_esti;
    myoutput(tt).Cyc_true = Cyc_true;
    myoutput(tt).InTwo_esti = InTwo_esti;
    myoutput(tt).InTwo_true = InTwo_true;
    myoutput(tt).OutTwo_esti = OutTwo_esti;
    myoutput(tt).OutTwo_true = OutTwo_true;
    myoutput(tt).Densi_esti = Densi_esti;
    myoutput(tt).Densi_true = Densi_true;
    myoutput(tt).score_list = score_list;

    % Compute residual diagnostics using the complete fitted residual tensor.
    res = double(X_output - X_observe_full);
    res = res(:);
    warning('off','all');
    [h,p] = lillietest(res);
    myoutput(tt).normal_p = p;
    e = res;
    lag = 20;
    [~,pValue,~,~] = lbqtest(e,'Lags',lag,'Alpha',0.05);
    myoutput(tt).LjungBox = pValue;
    sprintf('estimate %d-th sample', tt)
end

figure_handles = gobjects(7,1);

figure_handles(1) = figure;
boxplot([vertcat(myoutput.Rec_esti),vertcat(myoutput.Rec_true)])
set(gca, 'XTickLabel', {'Rec_esti','Rec_true'});
title('Comparison of Estimated and True Reciprocity');

figure_handles(2) = figure;
boxplot([vertcat(myoutput.Cyc_esti),vertcat(myoutput.Cyc_true)])
set(gca, 'XTickLabel', {'Cyc_esti','Cyc_true'});
title('Comparison of Estimated and True Cyclic triads');

figure_handles(3) = figure;
boxplot([vertcat(myoutput.InTwo_esti),vertcat(myoutput.InTwo_true)])
set(gca, 'XTickLabel', {'InTwo_esti','InTwo_true'});
title('Comparison of Estimated and True In-two-stars');

figure_handles(4) = figure;
boxplot([vertcat(myoutput.OutTwo_esti),vertcat(myoutput.OutTwo_true)])
set(gca, 'XTickLabel', {'OutTwo_esti','OutTwo_true'});
title('Comparison of Estimated and True Out-two-stars');

figure_handles(5) = figure;
boxplot([vertcat(myoutput.Densi_esti), vertcat(myoutput.Densi_true)])
set(gca,'XTickLabel',{'Densi_esti','Densi_true'});
title('Comparison of Estimated and True Edge Density');

figure_handles(6) = figure;
boxplot([vertcat(myoutput.normal_p)])
title('Boxplot of p_{normal}');

figure_handles(7) = figure;
boxplot([vertcat(myoutput.LjungBox)])
title('Boxplot of p_{ind}');

fprintf('MSE = %.4f(%.4f)\n', mean(vertcat(myoutput.MSE)),std(vertcat(myoutput.MSE)))
fprintf('R2 = %.4f(%.4f)\n', mean(vertcat(myoutput.R2)),std(vertcat(myoutput.R2)))

%% Save numerical outputs and all seven Figure 5 diagnostic panels
settings = struct('rep',rep,'m',m,'n',n,'s1',s1,...
    's2',s2,'K',K,'sigma',sigma,'missing_rate',w,'core_scale',N,...
    't_start',t_start,'t_end',t_end,'seq',seq,'rank',r,...
    'max_backtracking',max_m,'tolerance',tol,'max_iterations',maxiter,...
    'alpha',alpha);
summary = struct('MSE_mean',mean(vertcat(myoutput.MSE)),...
    'MSE_sd',std(vertcat(myoutput.MSE)),...
    'R2_mean',mean(vertcat(myoutput.R2)),...
    'R2_sd',std(vertcat(myoutput.R2)));

save(fullfile(results_dir,'figure5_results.mat'),...
    'myoutput','settings','summary','-v7.3');

figure_names = {'reciprocity','cyclic_triads','in_two_stars',...
    'out_two_stars','edge_density','residual_normality',...
    'residual_independence'};
for ii = 1:numel(figure_handles)
    savefig(figure_handles(ii),fullfile(results_dir,[figure_names{ii},'.fig']));
    print(figure_handles(ii),fullfile(results_dir,[figure_names{ii},'.png']),...
        '-dpng','-r300');
end

summary_file = fullfile(results_dir,'run_summary.txt');
fid = fopen(summary_file,'w');
if fid ~= -1
    fprintf(fid,'Functional-Edged Network Modeling simulation summary\n');
    fprintf(fid,'Replications: %d\n',rep);
    fprintf(fid,'Tensor dimensions: %d x %d x %d\n',n,m,length(seqtime));
    fprintf(fid,'Tucker rank: [%d, %d, %d]\n',r(1),r(2),r(3));
    fprintf(fid,'Missing rate: %.4f\n',w);
    fprintf(fid,'Noise variance: %.4f\n',sigma^2);
    fprintf(fid,'MSE: %.6f (%.6f)\n',summary.MSE_mean,summary.MSE_sd);
    fprintf(fid,'R2: %.6f (%.6f)\n',summary.R2_mean,summary.R2_sd);
    fclose(fid);
end

warning(warning_state);

% myoutput(tt).X_output is the final estimate for replication tt.
