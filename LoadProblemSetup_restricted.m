addpath("data/", "HelperFunctions/", "ReconScripts/", "EvalFunctions/");
gan_bin_size = 4; %1 or 4;
ProblemSetup.angle_count=round(360/gan_bin_size);

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

ProblemSetup.N = round(2*ProblemSetup.recon_cyl_radius / ProblemSetup.voxel_width); 
N = ProblemSetup.N;

ProblemSetup.source_radius = 260; % cm
% for gan_bin = 4, there are 95 empty bins on both sides; need to increase
% distance from source until this is correctly lined up
ProblemSetup.source_2_object = 15*ProblemSetup.source_radius; 
ProblemSetup.target_slice = 12;   

if gan_bin_size==1
    [sinogram,sinogram_slice] = read_pct_sinogram('data/sinogram_pileup.txt', ...
        ProblemSetup.vbins, ProblemSetup.tbins, ProblemSetup.angle_count, ProblemSetup.target_slice);

    [reconFBP,reconFBP_slice] = ReadRecon('data/FBP_pileup.txt', ProblemSetup.target_slice, ProblemSetup.N);
    [recon10,reconx10_slice] = ReadRecon('data/x_10_pileup.txt', ProblemSetup.target_slice, ProblemSetup.N);
    ProblemSetup.PseuTarget = reconx10_slice(:);
elseif gan_bin_size==4
    if use_corrected
        [sinogram,sinogram_slice] = read_pct_sinogram('data/sinogram_corrected90.txt', ...
            ProblemSetup.vbins, ProblemSetup.tbins, ProblemSetup.angle_count, ProblemSetup.target_slice);
    
        [reconFBP,reconFBP_slice] = ReadRecon('data/FBP_corrected90.txt', ProblemSetup.target_slice, ProblemSetup.N);
        [reconx15,reconx15_slice] = ReadRecon('data/x_15_corrected90.txt', ProblemSetup.target_slice, ProblemSetup.N);
        ProblemSetup.PseuTarget = reconx15_slice(:);
    else
        [sinogram,sinogram_slice] = read_pct_sinogram('data/sinogram_pileup4deg.txt', ...
            ProblemSetup.vbins, ProblemSetup.tbins, ProblemSetup.angle_count, ProblemSetup.target_slice);

        [reconFBP,reconFBP_slice] = ReadRecon('data/FBP_pileup4deg.txt', ProblemSetup.target_slice, ProblemSetup.N);
        [reconx15,reconx15_slice] = ReadRecon('data/x_15_pileup4deg.txt', ProblemSetup.target_slice, ProblemSetup.N);
        ProblemSetup.PseuTarget = reconx15_slice(:);
    end
    
end
ProblemSetup.scan = 'incomplete';
A = FormA(ProblemSetup);  
simulated_sino = A*reconFBP_slice(:);

reformatted_sinogram_slice = flip(sinogram_slice,2)'; 

% verification
figure(4);
tiledlayout(1,2);
nexttile(1); imagesc(reshape(reformatted_sinogram_slice,ProblemSetup.tbins,ProblemSetup.angle_count )); title('data');
nexttile(2); imagesc(reshape(simulated_sino,ProblemSetup.tbins, ProblemSetup.angle_count)); title('simulated sinogram');

if gan_bin_size==1
    reformatted_sinogram_slice(:,:) = [];
    simulated_sino = reshape(simulated_sino,ProblemSetup.tbins, ProblemSetup.angle_count);
    simulated_sino(:,:) = [];
elseif gan_bin_size==4
    % ignore 32-38, 77-83
    reformatted_sinogram_slice(:,76:85) = [];
    reformatted_sinogram_slice(:,31:40) = [];
    simulated_sino = reshape(simulated_sino,ProblemSetup.tbins, ProblemSetup.angle_count);
    simulated_sino(:,76:85) = [];
    simulated_sino(:,31:40) = [];
    
    thetas = linspace(0,2*pi, ProblemSetup.angle_count+1); thetas(end) = [];
    thetas(76:85)=[]; % removing angles 304-340
    thetas(31:40) = []; % removing angles 124-160
    % atv_angles = mod([-10, 5, 20, 40, 40+ 40/3, 40+ 80/3,80]./360*2*pi, 2*pi);
    % atv_angles = mod([-10,0,10,20, 40,50,60,70,80]./360*2*pi, 2*pi);
    % atv_angles = mod([80,90,100,110, 130,140,150,160,170]./360*2*pi, 2*pi);
    atv_angles = mod([5, 50, 95, 170, 215, 260, 310, 355]./360*2*pi, 2*pi);
    % atv_angles = mod([0, 90,280,270]./360*2*pi, 2*pi);
    % remove A rows:
    A(75*ProblemSetup.tbins+1 : 85*ProblemSetup.tbins,:) = [];
    A(30*ProblemSetup.tbins+1 : 40*ProblemSetup.tbins,:) = [];
    ProblemSetup.angle_count = 70;
end

% reformatted_sinogram_slice = reformatted_sinogram_slice(:);

hull = FindHull(ProblemSetup, 115);
hull = hull(:)';

ProblemSetup.A = A;
ProblemSetup.projections = reformatted_sinogram_slice(:);
ProblemSetup.angles = thetas;
ProblemSetup.hull = hull; %ProblemSetup.A = A.*hull;

% metrics used for evaluation
% note, these are mirror from what is on the sample image, because our
% sinogram is flipped
Insert_Locations256 = struct();
Insert_Locations256.("air1") = [108, 43];
Insert_Locations256.("air2") = [147, 214];
Insert_Locations256.("tefl") = [67, 64];
Insert_Locations256.("pmp") = [155, 43];
Insert_Locations256.("ldpe") = [214, 108];
Insert_Locations256.("delr") = [43, 149];

insert_material = ["pmp", "ldpe", "poly", "acry", "delr", "tefl", "air1", "air2"]; %check teflon and pmp, within 1%
insert_value    = [.883,  .997,   1.024,   1.16,   1.359, 1.79,   0.00, 0.00];
% ldpe should be .979?
    % truth values from series2 paper
Insert_Truth = dictionary(insert_material, insert_value);

EvalParameters.Insert_Locations = Insert_Locations256;
EvalParameters.InsertRadius     = 5;
EvalParameters.Insert_Truth     = Insert_Truth;

