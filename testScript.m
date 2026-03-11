% file 1
clc; clear;
slice = 11;
num_slices = 24;
num_angles = 90;
slice_bins = 350;
min_bin = (slice-1)*num_angles*slice_bins;
max_bin = slice*num_angles*slice_bins-1;

% file = '/Users/cflin/Downloads/histories_orig.bin';
% 
% fid = fopen(file,'rb','ieee-le');
% n = fread(fid,1,'int32');
% t_0 = fread(fid, n, 'single');
% t_1 = fread(fid, n, 'single');
% t_2 = fread(fid, n, 'single');
% t_3 = fread(fid, n, 'single');
% v_0 = fread(fid, n, 'single');
% v_1 = fread(fid, n, 'single');
% v_2 = fread(fid, n, 'single');
% v_3 = fread(fid, n, 'single');
% u_0 = fread(fid, n, 'single');
% u_1 = fread(fid, n, 'single');
% u_2 = fread(fid, n, 'single');
% u_3 = fread(fid, n, 'single');
% WEPL_orig         = fread(fid, n, 'single'); WEPL_orig = WEPL_orig./10;
% gantry_angle = fread(fid, n, 'single');
% 
% fclose(fid);
% neg_WEPL = WEPL_orig<0;

% file 2
file = '/Users/cflin/Downloads/history_info.bin';
fid = fopen(file,'rb','ieee-le');
n = fread(fid,1,'int32');
bin_num_vector = fread(fid, n, 'int32');
WEPL_fbp = fread(fid, n, 'single');
x_entry_vector = fread(fid, n, 'single');
y_entry_vector = fread(fid, n, 'single');
z_entry_vector = fread(fid, n, 'single');
x_exit_vector = fread(fid, n, 'single');
y_exit_vector = fread(fid, n, 'single');
z_exit_vector = fread(fid, n, 'single');
xy_entry_angle_vector = fread(fid, n, 'single');
xz_entry_angle_vector = fread(fid, n, 'single');
xy_exit_angle_vector = fread(fid, n, 'single');
xz_exit_angle_vector = fread(fid, n, 'single');
fclose(fid);

slice_restrict = (bin_num_vector>=min_bin) & (bin_num_vector<=max_bin);

bin_num_vector = bin_num_vector(slice_restrict);
WEPL_fbp = WEPL_fbp(slice_restrict);
x_entry_vector = x_entry_vector(slice_restrict);
y_entry_vector = y_entry_vector(slice_restrict);
z_entry_vector = z_entry_vector(slice_restrict);
x_exit_vector = x_exit_vector(slice_restrict);
y_exit_vector = y_exit_vector(slice_restrict);
z_exit_vector = z_exit_vector(slice_restrict);
xy_entry_angle_vector = xy_entry_angle_vector(slice_restrict);
xz_entry_angle_vector = xz_entry_angle_vector(slice_restrict);
xy_exit_angle_vector = xy_exit_angle_vector(slice_restrict);
xz_exit_angle_vector = xz_exit_angle_vector(slice_restrict);

% restricting histories to only those that are binned to our slice
full_data = [bin_num_vector, WEPL_fbp, x_entry_vector, y_entry_vector, ...
                z_entry_vector, x_exit_vector, y_exit_vector, z_exit_vector, ...
                xy_entry_angle_vector, xz_entry_angle_vector, ...
                xy_exit_angle_vector, xz_exit_angle_vector];

file = '/Users/cflin/Downloads/sinogram.txt';
[sino, sino_slice] = read_pct_sinogram(file, 0, 0, 90, 11); %should be 2160x90=756000 bins?
file_name = '/Users/cflin/Downloads/bin_counts_post.txt';
bin_counts = readmatrix(file_name);

negative_bins = []; 
neg_sino_values = [];
neg_bin_counts = [];
for kk = 10*90+1:90*11%size(sino,1)
    for jj = 1:size(sino,2)
        if sino(kk,jj)<0
            b = (kk-1)*350 + jj - 1;
            negative_bins = [negative_bins, b];
            neg_sino_values = [neg_sino_values, sino(kk,jj)];
            neg_bin_counts = [neg_bin_counts, bin_counts(kk,jj)];
        end
    end
end
% min of 326030, max of 343887 -> 24828 histories remaining
neg_bin_hists = full_data(ismember(full_data(:,1),negative_bins),:);

% what happens if we restrict to valid histories?
sub1 = neg_bin_hists(neg_bin_hists(:,2) >= -20,:); % 17897 remaining

% check the angles: this displays all the histories
figure; hold on;
for ii = 1:100%size(metadata,1)
    plot([sub1(ii,3), sub1(ii,6)], [sub1(ii,4), sub1(ii,7)]);
end
hold off;

% look at total counts for these:



pause;
% loading 
file_name = '/Users/cflin/Downloads/OddHistories_Metadata.txt';
metadata = readmatrix(file_name);
figure; hold on;
for ii = 1:100%size(metadata,1)
    plot([metadata(ii,4), metadata(ii,7)], [metadata(ii,5), metadata(ii,8)]);
end
hold off;
% WEPL_vector (1), bin_num_vector (2), gantry_angle_orig (3), x_entry_vector (4), 
%     y_entry_vector (5), z_entry_vector (6), x_exit_vector (7), y_exit_vector (8), z_exit_vector (9),
%     xy_entry_angle_vector (10), xz_entry_angle_vector (11), xy_exit_angle_vector (12), 
%     xz_exit_angle_vector (13) , rel_ut_angle_vector (14), rel_uv_angle_vector (15)
ang = 48;
tmp = metadata(metadata(:,3) == ang,:);
figure; hold on;

for ii = 1:size(tmp,1)
    plot([tmp(ii,4), tmp(ii,7)], [tmp(ii,5), tmp(ii,8)]);
end
hold off;

neg_bins = []; 


neg_histories = metadata(ismember(metadata(:,2), negative_bins),:);
figure; hold on;
for ii = 1:size(neg_histories)
    if neg_histories(ii,1)>=-20
    plot([neg_histories(ii,4), neg_histories(ii,7)], [neg_histories(ii,5), neg_histories(ii,8)]);
    disp([neg_histories(ii,1), neg_histories(ii,2), neg_histories(ii,3)]); %display wepl, bin number, gantry
    end
end
hold off;

file_name = '/Users/cflin/Downloads/bin_counts_post.txt';
bin_counts = readmatrix(file_name);
bin_ct_slice = bin_counts((slice-1)*num_angles+1:slice*num_angles,:);
negative_bins = []; %keep slice 1-15
neg_sino_values = [];
neg_sinbin_counts = [];
for kk = 10*90+1:90*11%size(sino,1)
    for jj = 1:size(sino,2)
        if sino(kk,jj)<0
            b = (kk-1)*350 + jj - 1;
            negative_bins = [negative_bins, b];
            neg_sino_values = [neg_sino_values, sino(kk,jj)];
            neg_sinbin_counts = [neg_sinbin_counts, bin_counts(kk,jj)];
        end
    end
end