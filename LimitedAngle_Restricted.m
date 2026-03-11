close all;
clear;

use_corrected=true;
run("LoadProblemSetup_restricted.m")

sinogram_slice_flat = reformatted_sinogram_slice(:);
As = sum(A,2);
disp(sum((As==0) & (sinogram_slice_flat>0)));

% ---------------------------- ART-TV   Testing----------------------------
pm_ARTTV.Iterations     = 200; 
pm_ARTTV.GradDescSteps  = 30;
pm_ARTTV.initial        = zeros(ProblemSetup.N*ProblemSetup.N,1);
pm_ARTTV.lambda         = 0.05; 
pm_ARTTV.alpha          = .1; 
pm_ARTTV.epsilon        = 1e-6;

tic; arttv_recons = ART_TV(ProblemSetup, pm_ARTTV); toc; 

% ---------------------------- SART-ATV Testing----------------------------
pm_SARTATV.initial      = zeros(ProblemSetup.N*ProblemSetup.N, 1);
pm_SARTATV.alphas       = atv_angles; 
pm_SARTATV.iterations   = 200;
pm_SARTATV.epsilon      = 1e-6; % for gradient calculation
pm_SARTATV.lambda       = .075; % for SART stepsize
pm_SARTATV.tvStep       = 0.01; % stepsize for gradient
pm_SARTATV.N_tv         = 30;
pm_SARTATV.alphas       = atv_angles; 

tic; sartatv_recons = SART_ATV(ProblemSetup, pm_SARTATV); toc; 

% ---------------------------- L1dL2    Testing----------------------------
r1 = 1; r2 = 1;b=1; l = 0.05;
pm_L1dL2.rho1 = r1; 
pm_L1dL2.rho2 = r2; %perhaps .1 both
pm_L1dL2.beta = b; %maybe 1 or 0.1 
pm_L1dL2.lambda = l; 
pm_L1dL2.maxit = 200; 
pm_L1dL2.Imaxit = 5; 
pm_L1dL2.Itol = 1e-4; 
pm_L1dL2.tol = 1e-5;
pm_L1dL2.upper_bound = 3;
pm_L1dL2.lower_bound = 0;
pm_L1dL2.initial = zeros(N*N,1);

tic; l1dl2_recons = L1L2_ADMM(ProblemSetup, pm_L1dL2); toc;

% -------------------------------------------------------------------------
% evaluate sart-tv
% figure(6); imshow(reshape(arttv_recons(:,end),N,N),[],'InitialMagnification','fit'); hold on;
% arttv_metrics = eval_sensitom(arttv_recons, ProblemSetup.N, pm_ARTTV.Iterations, EvalParameters);
% 
% figure(7); f1 = tiledlayout(3,2);
% nexttile(1); plot([arttv_metrics.air1_recon]); title('air1'); ylim([0,1]);
% nexttile(2); plot([arttv_metrics.air2_recon]); title('air2'); ylim([0,1]);
% nexttile(3); plot([arttv_metrics.delr_error]); title('delr'); ylim([0,1]);
% nexttile(4); plot([arttv_metrics.ldpe_error]); title('ldpe'); ylim([0,1]);
% nexttile(5); plot([arttv_metrics.pmp_error]); title('pmp'); ylim([0,1]);
% nexttile(6); plot([arttv_metrics.tefl_error]); title('tefl'); ylim([0,1]);
% title(f1,'ARTTV-Eval')
% 
% % evaluate sart-atv
% sartatv_metrics = eval_sensitom(sartatv_recons, ProblemSetup.N, pm_SARTATV.iterations, EvalParameters);
% figure(8); f1 = tiledlayout(3,2);
% nexttile(1); plot([sartatv_metrics.air1_recon]); title('air1'); ylim([0,1]);
% nexttile(2); plot([sartatv_metrics.air2_recon]); title('air2'); ylim([0,1]);
% nexttile(3); plot([sartatv_metrics.delr_error]); title('delr'); ylim([0,1]);
% nexttile(4); plot([sartatv_metrics.ldpe_error]); title('ldpe'); ylim([0,1]);
% nexttile(5); plot([sartatv_metrics.pmp_error]); title('pmp'); ylim([0,1]);
% nexttile(6); plot([sartatv_metrics.tefl_error]); title('tefl'); ylim([0,1]);
% title(f1,'SARTATV-Eval')

% evaluate l1dl2
l1dl2_metrics = eval_sensitom(l1dl2_recons, ProblemSetup.N, size(l1dl2_recons,2),EvalParameters);
figure(9); f1 = tiledlayout(3,2);
nexttile(1); plot([l1dl2_metrics.air1_recon]); title('air1'); ylim([0,1]);
nexttile(2); plot([l1dl2_metrics.air2_recon]); title('air2'); ylim([0,1]);
nexttile(3); plot([l1dl2_metrics.delr_error]); title('delr'); ylim([0,1]);
nexttile(4); plot([l1dl2_metrics.ldpe_error]); title('ldpe'); ylim([0,1]);
nexttile(5); plot([l1dl2_metrics.pmp_error]); title('pmp'); ylim([0,1]);
nexttile(6); plot([l1dl2_metrics.tefl_error]); title('tefl'); ylim([0,1]);
title(f1,'L1dL2-Eval')
% evaluate art-tv
% arttv_metrics = eval_sensitom(arttv_recons, ProblemSetup.N, ,EvalParameters);
figure(10); f1 = tiledlayout(3,2); 
nexttile(1); hold on; title('Air (Top)'); ylim([0,1]);
plot([arttv_metrics.air1_recon],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.air1_recon],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.air1_recon],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit); 
legend({'ARTTV','SARTATV','L1dL2'});
hold off;
nexttile(2); hold on; title('Air (Bottom)'); ylim([0,1]);
plot([arttv_metrics.air2_recon],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.air2_recon],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.air2_recon],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit); 
legend({'ARTTV','SARTATV','L1dL2'});
hold off;
nexttile(3); hold on; title('Delrin'); ylim([0,1]);
plot([arttv_metrics.delr_error],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.delr_error],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.delr_error],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit); 
legend({'ARTTV','SARTATV','L1dL2'});
hold off;
nexttile(4); hold on; title('LDPE'); ylim([0,1]);
plot([arttv_metrics.ldpe_error],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.ldpe_error],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.ldpe_error],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit); 
legend({'ARTTV','SARTATV','L1dL2'});
hold off;
nexttile(5); hold on; title('PMP'); ylim([0,1]);
plot([arttv_metrics.pmp_error],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.pmp_error],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.pmp_error],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit); 
legend({'ARTTV','SARTATV','L1dL2'});
hold off;
nexttile(6); hold on; title('Teflon'); ylim([0,1]);
plot([arttv_metrics.tefl_error],'Color','k','LineWidth',2,'MarkerIndices',1:25:pm_ARTTV.Iterations); 
plot([sartatv_metrics.tefl_error],'Color','g','LineWidth',2,'Marker','x','MarkerIndices',1:25:pm_SARTATV.iterations); 
plot([l1dl2_metrics.tefl_error],'Color','b','LineWidth',2,'Marker','o','MarkerIndices',1:25:pm_L1dL2.maxit);
legend({'ARTTV','SARTATV','L1dL2'});
hold off;

% -------------------------------------------------------------------------
% visualize all results
figure(11); t = tiledlayout(1,3); 
nexttile(1); imshow(reshape(arttv_recons(:,end),N,N),[],'InitialMagnification','fit'); title('ART-TV Reconstruction');
nexttile(2); imshow(reshape(sartatv_recons(:,end),N,N),[],'InitialMagnification','fit'); title('SART-ATV Reconstruction');
nexttile(3); imshow(reshape(l1dl2_recons(:,end),N,N),[],'InitialMagnification','fit'); title('L1dL2 Reconstruction');