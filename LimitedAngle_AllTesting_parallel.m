%% Testing on existing pCT data
% problem setup
% close all;
% clc;clear;

addpath("data/", "HelperFunctions/", "ReconScripts/");
run("LoadProblemSetup.m")

[sinogram,sinogram_slice] = read_pct_sinogram('data/pct_sinogram.txt', ...
    ProblemSetup.vbins, ProblemSetup.tbins, ProblemSetup.angle_count, ProblemSetup.target_slice);
% % hull = findHull(sinogram, angles);

% hull = FindHull(sinogram, projections, 100);


if size(sinogram,1)~=ProblemSetup.angle_count * ProblemSetup.vbins || ...
        size(sinogram,2)~=ProblemSetup.tbins
    disp('mismatch of sinogram')
    return;
end

desired_angles = 1:ProblemSetup.angle_count;
remove_angles = 60:120;
angle_subset = desired_angles; desired_angles(remove_angles) = [];
tic;

A = FormA(ProblemSetup);

phan320 = phantom(320);
phansin = A*phan320(:);

toc;
hull = FindHull(ProblemSetup);
hull = hull(:)';

% A = A*spdiags(hull, 0, length(hull), length(hull));

% verification
sinogram_slice_flat = sinogram_slice(:);
N = ProblemSetup.N;
As = sum(A,2);
disp(sum((As==0) & (sinogram_slice_flat>0)));

ProblemSetup.projections = sinogram_slice_flat;
    ProblemSetup.projections = phansin;
ProblemSetup.A = A;
ProblemSetup.hull = hull;

pm_ARTTV.Iterations = 10;
pm_ARTTV.GradDescSteps = 30;
pm_ARTTV.initial = zeros(ProblemSetup.N*ProblemSetup.N,1);

pm_ARTTV.lambda = 5; %?? too big? 
pm_ARTTV.alpha = 0.6;
pm_ARTTV.epsilon = 0.1;

% tic; recons = ART_TV(ProblemSetup, pm_ARTTV); toc;
    


pm_ARTATV.initial = zeros(ProblemSetup.N*ProblemSetup.N, 1);
pm_ARTATV.N_tv = 20;
pm_ARTATV.alphas = linspace(0,2*pi, 9); pm_ARTATV.alphas(end) = []; 
        %use N_alpha=8, so 8 angles within valid scan range
pm_ARTATV.iterations = 2000;
pm_ARTATV.epsilon = .1; 
pm_ARTATV.lambda = 1.1;
pm_ARTATV.tvStep = 0.05;
ProblemSetup.angles = linspace(0,2*pi, ProblemSetup.angle_count+1); ProblemSetup.angles(end) = [];
recons = ART_ATV(ProblemSetup, pm_ARTATV);

pm_ARTTV.lambda = 0.75; %?? too big? 
pm_ARTTV.alpha = 0.6;
pm_ARTTV.epsilon = 0.1;

tic; recons = ART_TV(ProblemSetup, pm_ARTTV); toc;
 
% pm_ARTATV.initial = zeros(1, ProblemSetup.N*ProblemSetup.N);
% pm_ARTATV.N_tv = 20;
% pm_ARTATV.alphas = linspace(0,2*pi, 9); pm_ARTATV.alphas(end) = []; 
%         %use N_alpha=8, so 8 angles within valid scan range
% pm_ARTATV.iterations = 10;
% pm_ARTATV.epsilon = .1; 
% pm_ARTATV.lambda = .8;
% pm_ARTATV.a = 0.05;

