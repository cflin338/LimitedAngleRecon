% loading in generic variables 
addpath("data/", "HelperFunctions/", "ReconScripts/");

seed = 1; rng(seed); 

ProblemSetup.angle_count=360;

ProblemSetup.ssd_tsize = 35; % cm
ProblemSetup.tbinsize=0.1;   % cm
ProblemSetup.tbins = floor(ProblemSetup.ssd_tsize / ProblemSetup.tbinsize+0.5);
ProblemSetup.ssd_vsize = 6;  % cm
ProblemSetup.vbinsize=0.25;
ProblemSetup.vbins = floor(ProblemSetup.ssd_vsize / ProblemSetup.vbinsize+0.5);

ProblemSetup.recon_cyl_radius = 8;   % cm
ProblemSetup.voxel_width = 0.0625;   % cm
ProblemSetup.voxel_height = 0.0625;  % cm
ProblemSetup.voxel_thickness = 0.25; % cm
% voxels
ProblemSetup.N = round(2*ProblemSetup.recon_cyl_radius / ProblemSetup.voxel_width); 
% reconstruction width

ProblemSetup.source_radius = 260; % cm
ProblemSetup.source_2_object = 15*ProblemSetup.source_radius;
ProblemSetup.target_slice = 12;   

N = ProblemSetup.N;

[sinogram,sinogram_slice] = read_pct_sinogram('data/pct_sinogram.txt', ...
    ProblemSetup.vbins, ProblemSetup.tbins, ProblemSetup.angle_count, ProblemSetup.target_slice);
% top 2 and bottom 2 slices all read as zero; only middle 20 contain data

% verifying data, problem setup
[reconFBP,reconFBP_slice] = ReadRecon('data/FBP.txt', ProblemSetup.target_slice, ProblemSetup.N);
[reconx30,reconx30_slice] = ReadRecon('data/x_30.txt', ProblemSetup.target_slice, ProblemSetup.N);

if size(sinogram,1)~=ProblemSetup.angle_count * ProblemSetup.vbins || ...
        size(sinogram,2)~=ProblemSetup.tbins
    disp('mismatch of sinogram')
    return;
end

A = FormA(ProblemSetup);  
simulated_sino = A*reconFBP_slice(:);

reformatted_sinogram_slice = flip(sinogram_slice,2)'; 
reformatted_sinogram_slice = reformatted_sinogram_slice(:);

