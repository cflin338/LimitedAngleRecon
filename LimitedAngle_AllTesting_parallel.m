%% 
% problem setup
close all;
clc;clear;

run("LoadProblemSetup.m")

% below, some variables to loop through different reconstruction scenarios
init_list               = "zeros";
% if ProblemSetup.AngularShift==0
%     init_list = [init_list,"guo"];
% end

max_angles = [90, 120, 150];
angle_list = cell(1,length(max_angles));
% Addiing the structure to the variables collecting all results
for max_angle = max_angles
    field_name = ['a', num2str(max_angle)];
    all_metrics.(field_name) = struct();
    for init_type = init_list
        all_metrics.(field_name).(init_type) = struct('arttv',struct(), ...
                                             'artatv',struct(),...
                                             'l1dl2',struct());
        all_reconstructions.(field_name).(init_type) = struct('arttv',struct('sub',struct(),'full',struct()), ...
                                                     'artatv',struct('sub',struct(),'full',struct()),...
                                                     'l1dl2',struct('sub',struct(),'full',struct()));
    end
end

% creating empty struct result holder for parallelization
for mang_idx = 1:length(max_angles)
    max_angle = max_angles(mang_idx);
    angle_field_name = ['a', num2str(max_angle)];
    for init_type = init_list

        all_reconstructions.(angle_field_name).(init_type).('l1dl2') = [];
        all_samples.(angle_field_name).(init_type).('l1dl2') = [];
        all_metrics.(angle_field_name).(init_type).('l1dl2').('full') = [];
        all_metrics.(angle_field_name).(init_type).('l1dl2').('sub') = [];

        all_reconstructions.(angle_field_name).(init_type).('arttv') = [];
        all_samples.(angle_field_name).(init_type).('arttv') = [];
        all_metrics.(angle_field_name).(init_type).('arttv').('full') = [];
        all_metrics.(angle_field_name).(init_type).('arttv').('sub') = [];

        all_reconstructions.(angle_field_name).(init_type).('artatv') = [];
        all_samples.(angle_field_name).(init_type).('artatv') = [];
        all_metrics.(angle_field_name).(init_type).('artatv').('full') = [];
        all_metrics.(angle_field_name).(init_type).('artatv').('sub') = []; 
    end
end

%%
% iterate through each range of scans; perform reconstructions
for mang_idx = 1:length(max_angles)
    max_angle = max_angles(mang_idx);
    angle_field_name = ['a', num2str(max_angle)];
    angle_list{mang_idx} = angle_field_name;
    
    Max_angle       = max_angle; 
    Min_angle       = 0;
    num_projections = max_angle;
    Angle_Range     = Max_angle - Min_angle;
    Initial_angle   = Angle_Range/num_projections - ProblemSetup.AngularShift;
    Angle_increment = ProblemSetup.Angle_increment; 
    Max_angle       = Max_angle - ProblemSetup.AngularShift;
    angles          = Initial_angle:Angle_increment:Max_angle;

    ProblemSetup.angles      = angles;
    [A, projections, img]    = PRtomo(ProblemSetup);

    ProblemSetup.A           = A;
    ProblemSetup.projections = projections;
    ProblemSetup.img         = img;
    
    %% 
    % select inital guess
    initial_guesses.zeros    = zeros(ProblemSetup.N*ProblemSetup.N,1);
    if ProblemSetup.AngularShift==0
        initial_guesses.guo  = Guo_InitialGuess_v2(ProblemSetup,10,10, false);
    end

    %removing empty rows
    rowsToRemove = sum(A, 2) == 0;
    A(rowsToRemove, :) = [];
    projections(rowsToRemove) = [];
    
    ProblemSetup.A           = A;
    ProblemSetup.projections = projections;
    ProblemSetup.img         = img;
    % randomizing order of performing projection (ART step in Chen and Sidky)
    ProblemSetup.proj_order  = randperm(length(ProblemSetup.projections));
    %---------------------------------------
    %%
    % miscelleneous settings
    visualize_11dl2 = false;
    visualize_sidky = false;
    visualize_chen  = false;
    chen_prog       = false;
    sidky_prog      = false;
    l1dl2_prog      = false;
    
    %% 
    % specific setups-L1dL2, ADMM
    pm_L1dL2.maxit       = max_iters; % from paper
    pm_L1dL2.tol         = 1e-5; %from code
    pm_L1dL2.Itol        = 1e-4; %from wang code
    pm_L1dL2.Imaxit      = 5; % paper uses 1, 3, 5, 10
    pm_L1dL2.upper_bound = 1; % due to phantom min/max
    pm_L1dL2.lower_bound = 0; % due to phantom min/max
    pm_L1dL2.sample_rate = sample_rate;
    switch Max_angle
        % below are values taken from Wang:
        % https://github.com/wangcmath/limited_angle_CT_L1dL2/blob/main/demo_noiseless_SL.m
        % currently restricted to angles: 30, 45, 60, 150, 180. May be
        %       unoptimal for others. "Otherwise" values arbitrarily selected
        case 30
            pm_L1dL2.lambda =0.05;
            pm_L1dL2.rho1 = 1;pm_L1dL2.rho2 = 1;pm_L1dL2.beta = 1; 
        case 45
            pm_L1dL2.lambda =0.05;pm_L1dL2.rho1 = 1;pm_L1dL2.rho2 = 1;pm_L1dL2.beta = 1; % 1.733e-02
        case 60
            pm_L1dL2.lambda =0.05;pm_L1dL2.rho1 = 1;pm_L1dL2.rho2 = 1;pm_L1dL2.beta = .5; 
        case 90
            pm_L1dL2.lambda =0.05; pm_L1dL2.beta = .1; pm_L1dL2.rho1 = .1;pm_L1dL2.rho2 = pm_L1dL2.rho1;
        case 150
            pm_L1dL2.lambda =0.05; pm_L1dL2.beta = 1; pm_L1dL2.rho1 = 1;pm_L1dL2.rho2 = pm_L1dL2.rho1;
        case 180
            pm_L1dL2.lambda =.5;pm_L1dL2.rho1 = 10;pm_L1dL2.rho2 = 1; pm_L1dL2.beta = 10;
        otherwise
            pm_L1dL2.lambda =0.05;pm_L1dL2.rho1 = 1;pm_L1dL2.rho2 = 1; pm_L1dL2.beta = 1;
    end
    
    pm_L1dL2.disp_prog  = l1dl2_prog;
    pm_L1dL2.visualize   = visualize_11dl2;
    
    %% 
    % specific setups-Sidky, ARTTV
    pm_ART_TV.alpha         = 0.2;      %paper value
    pm_ART_TV.Iterations    = max_iters;      %paper value for lim-ang  example
    pm_ART_TV.NData         = numel(projections);
    pm_ART_TV.GradDescSteps = 20;       %paper value
    pm_ART_TV.exit_error    = 1e-5; %10e-8;    %not in paper, taken from online code
    pm_ART_TV.epsilon       = 10e-8;    %paper value
    pm_ART_TV.sample_rate   = sample_rate;
    pm_ART_TV.disp_prog     = sidky_prog;
    pm_ART_TV.visualize     = visualize_sidky;
    %%
    % specific setups-Chen, ARTATV
    Na = 8;
    pm_ART_ATV.alphas = Initial_angle:Angle_Range/(Na-1):Max_angle;
    pm_ART_ATV.alphas = pm_ART_ATV.alphas + 90;
    
    pm_ART_ATV.iterations    = max_iters;           % paper value, uses between 200 and 10000 in tables
    pm_ART_ATV.N_tv          = 20;            % paper value            
    pm_ART_ATV.N_projections = num_projections;
    pm_ART_ATV.alpha         = 0.15;          % paper value, .2 in other paper, for grad descent
    if max_angle==150
        pm_ART_ATV.lambda        = 0.8;
    elseif max_angle==120
        pm_ART_ATV.lambda        = 1;
    elseif max_angle == 90
        pm_ART_ATV.lambda        = 1.6; %what if i do higher than 1.4? changing to 1.6
    else
        pm_ART_ATV.lambda        = 0.75;
    end
    pm_ART_ATV.epsilon       = 1e-5; %1e-8;          % value from other code
    pm_ART_ATV.tv_epsilon    = 1e-8;          % same value from Sidky
    pm_ART_ATV.thetas        = angles;
    pm_ART_ATV.sample_rate   = sample_rate;
    
    pm_ART_ATV.disp_prog     = chen_prog;
    pm_ART_ATV.visualize     = visualize_chen;
    
    %%
    for init_idx = 1:length(init_list)
        init_type = init_list(init_idx);
        fprintf('--%i, %s--\n', max_angle, init_type);
        pm_L1dL2.initial                = initial_guesses.(init_type);
        pm_ART_TV.initial               = initial_guesses.(init_type);
        pm_ART_ATV.initial              = initial_guesses.(init_type);

        par_recon = cell(1,3);
        par_met = cell(1,3);
        par_samp = cell(1,3);
        
        parfor i = 1:3
            if mod(i,3)==1
                % corresponds with l1dl2
                disp('L1dL2')
                tic;
                [l1dl2_recon, l1dl2_metrics, l1dl2_samples]    = L1L2_ADMM(ProblemSetup,pm_L1dL2);
                par_recon{i} = l1dl2_recon;
                par_met{i}   = l1dl2_metrics;
                par_samp{i}  = l1dl2_samples;
                elapsed_time = toc;
                fprintf('L1dL2 %.2f\n', elapsed_time);
            elseif mod(i,3)==2
                % corresponds with ARTTV
                disp('ART-TV')
                tic;
                [arttv_recon, arttv_metrics, arttv_samples]  = ART_TV(ProblemSetup, pm_ART_TV);
                par_recon{i} = arttv_recon;
                par_met{i}   = arttv_metrics;
                par_samp{i}  = arttv_samples;
                elapsed_time = toc;
                fprintf('ART-TV %.2f\n', elapsed_time);
            else
                % corresponds with ARTATV
                disp('ART-ATV')
                tic;
                [artatv_recon,artatv_metrics, artatv_samples] = ART_ATV(ProblemSetup, pm_ART_ATV);
                par_recon{i} = artatv_recon;
                par_met{i}   = artatv_metrics;
                par_samp{i}  = artatv_samples;
                elapsed_time = toc;
                fprintf('ART-ATV %.2f\n', elapsed_time);
            end
        end
        
        all_reconstructions.(angle_field_name).(init_type).('l1dl2') = par_recon{1};
        all_samples.(angle_field_name).(init_type).('l1dl2') = par_samp{1};
        all_metrics.(angle_field_name).(init_type).('l1dl2').('full') = par_met{1}.full;
        all_metrics.(angle_field_name).(init_type).('l1dl2').('sub') = par_met{1}.sub;

        all_reconstructions.(angle_field_name).(init_type).('arttv') = par_recon{2};
        all_samples.(angle_field_name).(init_type).('arttv') = par_samp{2};
        all_metrics.(angle_field_name).(init_type).('arttv').('full') = par_met{2}.full;
        all_metrics.(angle_field_name).(init_type).('arttv').('sub') = par_met{2}.sub;

        all_reconstructions.(angle_field_name).(init_type).('artatv') = par_recon{3};
        all_samples.(angle_field_name).(init_type).('artatv') = par_samp{3};
        all_metrics.(angle_field_name).(init_type).('artatv').('full') = par_met{3}.full;
        all_metrics.(angle_field_name).(init_type).('artatv').('sub') = par_met{3}.sub;
    end
end
%%
date_str                = datestr(datetime('today'), 'mmddyyyy');  
if length(max_angles)==3
    ang_str = "all_angles";
else
    ang_str = int2str(max_angles(end));
end
filename                = sprintf('./Results/full_run_%ishift_%s_%s.mat', ProblemSetup.AngularShift, date_str, ang_str);
save( filename, 'all_reconstructions', 'all_samples', 'all_metrics', 'ProblemSetup');
methods                 = ["l1dl2", "arttv", "artatv"];
sample_rates            = struct(); 
sample_rates.l1dl2      = pm_L1dL2.sample_rate;
sample_rates.arttv     = pm_ART_TV.sample_rate;
sample_rates.artatv    = pm_ART_ATV.sample_rate;

% ResultEvaluation(all_samples, angle_list, methods, init_list, ...
%                                 8, ProblemSetup, sample_rates);
                            
                            
% run("ViewFinalReconMetrics.m");
