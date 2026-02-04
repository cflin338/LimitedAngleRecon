% loading in generic variables 
seed = 1; rng(seed); 

N                          = 128;
max_iters                  = 2000;% 10000;
sample_rate                = 1;

ProblemSetup.N             = N;
ProblemSetup.bins          = 2*N;
ProblemSetup.Angle_increment = 1;
ProblemSetup.width         = round(sqrt(2)*ProblemSetup.N)-1;
ProblemSetup.empty_sub     = struct('t1', [], 't2',[], 't3',[], 't4',[], ...
                                    't5',[], 't6',[],'t7',[], 't8',[]);
ProblemSetup.MetricMode    = "loganshep";
ProblemSetup.sample_rate   = sample_rate;
ProblemSetup.SystemMatrix  = "ChordLength";
% ProblemSetup.SystemMatrix  = "Intercept";
ProblemSetup.AngularShift  = 0;

RegionSearchSettings.N          = ProblemSetup.N;
RegionSearchSettings.large_rad  = 8;
RegionSearchSettings.large_step = 4;
RegionSearchSettings.med_rad    = 8;
RegionSearchSettings.med_step   = 4;
RegionSearchSettings.small_rad  = 4;
RegionSearchSettings.small_step = 1;

[targ_r, targ_c] = FindLogShepRegions(RegionSearchSettings);
ProblemSetup.target_rows = targ_r;
ProblemSetup.target_cols = targ_c;

% max_angle -> 
%              initialization type -> 
%                                     method ->
%                                               sub or overall
all_metrics         = struct(); 
all_reconstructions = struct();
all_samples         = struct();

