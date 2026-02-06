% loading in generic variables 
seed = 1; rng(seed); 

ProblemSetup.angle_count=360;

ProblemSetup.ssd_tsize = 35;
ProblemSetup.tbinsize=0.1;
ProblemSetup.tbins = floor(ProblemSetup.ssd_tsize / ProblemSetup.tbinsize+0.5);
ProblemSetup.ssd_vsize = 6;
ProblemSetup.vbinsize=0.25;
ProblemSetup.vbins = floor(ProblemSetup.ssd_vsize / ProblemSetup.vbinsize+0.5);

ProblemSetup.recon_cyl_radius = 10;
ProblemSetup.voxel_width = 0.0625;
ProblemSetup.voxel_height = 0.0625;
ProblemSetup.voxel_thickness = 0.25; 
ProblemSetup.N = round(2*ProblemSetup.recon_cyl_radius / ProblemSetup.voxel_width); 
% reconstruction width

ProblemSetup.source_radius = 260;
ProblemSetup.target_slice = 12;