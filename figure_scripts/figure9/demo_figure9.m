% demo_figure9.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/10/2025, Last modified: 07/10/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
package_path = 'D:\cartesian_maxgirf_epi_3d';

%% Add source directories to search path
addpath(genpath(package_path));

%% Define the full path of an output directory
%--------------------------------------------------------------------------
% Fully sampled 2D SE-EPI (coronal)
%--------------------------------------------------------------------------
output_path1 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding0_phc0_cfc0_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% Accelerated 2D SE-EPI (axial)
%--------------------------------------------------------------------------
output_path2 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% Diffusion-weighted 2D SE-EPI (axial)
%--------------------------------------------------------------------------
output_path3 = 'D:\cartesian_maxgirf_epi_2d\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% Long-ETL 3D GRE-EPI (axial)
%--------------------------------------------------------------------------
output_path4 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% High-resolution 3D GRE-EPI (axial)
%--------------------------------------------------------------------------
output_path5 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_highres_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% High-resolution 3D GRE-EPI in-vivo data (axial)
%--------------------------------------------------------------------------
output_path6 = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc1_sfc0_gnc1_topup0';

%% Read a .cfl file
Nkx1 = 256;
Nky1 = 128;
nr_slices1 = 20;

parabolic_shift1 = zeros(Nkx1, Nky1 ,nr_slices1, 'single');

x1 = zeros(Nkx1, Nky1, nr_slices1, 'single');
y1 = zeros(Nkx1, Nky1, nr_slices1, 'single');
z1 = zeros(Nkx1, Nky1, nr_slices1, 'single');

for slice_number = 1:nr_slices1
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    if nr_slices1 > 1 % multi-slice
        if mod(nr_slices1,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices1 / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices1 / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end

    %----------------------------------------------------------------------
    % parabolic_shift (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('parabolic_shift_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    parabolic_shift1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('x_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('y_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('z_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Remove readout oversampling
Nx1 = Nkx1 / 2;
Ny1 = Nky1;
Nz1 = nr_slices1;

idx1_range = (-floor(Nx1/2):ceil(Nx1/2)-1).' + floor(Nkx1/2) + 1;

parabolic_shift1 = parabolic_shift1(idx1_range, :, :);

x1 = x1(idx1_range, :, :);
y1 = y1(idx1_range, :, :);
z1 = z1(idx1_range, :, :);

%% Read a .cfl file
Nkx2 = 512;
Nky2 = 256;
nr_slices2 = 56;

parabolic_shift2 = zeros(Nkx2, Nky2 ,nr_slices2, 'single');

x2 = zeros(Nkx2, Nky2, nr_slices2, 'single');
y2 = zeros(Nkx2, Nky2, nr_slices2, 'single');
z2 = zeros(Nkx2, Nky2, nr_slices2, 'single');

for slice_number = 1:nr_slices2
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    if nr_slices2 > 1 % multi-slice
        if mod(nr_slices2,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices2 / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices2 / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end

    %----------------------------------------------------------------------
    % parabolic_shift (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path2, sprintf('parabolic_shift_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    parabolic_shift2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path2, sprintf('x_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path2, sprintf('y_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path2, sprintf('z_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Remove readout oversampling
Nx2 = Nkx2 / 2;
Ny2 = Nky2;
Nz2 = nr_slices2;

idx1_range = (-floor(Nx2/2):ceil(Nx2/2)-1).' + floor(Nkx2/2) + 1;

parabolic_shift2 = parabolic_shift2(idx1_range, :, :);

x2 = x2(idx1_range, :, :);
y2 = y2(idx1_range, :, :);
z2 = z2(idx1_range, :, :);

%% Read a .cfl file
Nkx3 = 352;
Nky3 = 176;
nr_slices3 = 20;

parabolic_shift3 = zeros(Nkx3, Nky3 ,nr_slices3, 'single');

x3 = zeros(Nkx3, Nky3, nr_slices3, 'single');
y3 = zeros(Nkx3, Nky3, nr_slices3, 'single');
z3 = zeros(Nkx3, Nky3, nr_slices3, 'single');

for slice_number = 1:nr_slices3
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    if nr_slices3 > 1 % multi-slice
        if mod(nr_slices3,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices3 / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices3 / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end

    %----------------------------------------------------------------------
    % parabolic_shift (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path3, sprintf('parabolic_shift_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    parabolic_shift3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path3, sprintf('x_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path3, sprintf('y_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path3, sprintf('z_slc%d_flat', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Remove readout oversampling
Nx3 = Nkx3 / 2;
Ny3 = Nky3;
Nz3 = nr_slices3;

idx1_range = (-floor(Nx3/2):ceil(Nx3/2)-1).' + floor(Nkx3/2) + 1;

parabolic_shift3 = parabolic_shift3(idx1_range, :, :);

x3 = x3(idx1_range, :, :);
y3 = y3(idx1_range, :, :);
z3 = z3(idx1_range, :, :);

%% Read a .cfl file
%--------------------------------------------------------------------------
% parabolic_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'parabolic_shift');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'x_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'y_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'z_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Remove readout oversampling
[Nkx4,Nky4,Nkz4] = size(parabolic_shift4);

Nx4 = Nkx4 / 2;
Ny4 = Nky4;
Nz4 = Nkz4;

idx1_range = (-floor(Nx4/2):ceil(Nx4/2)-1).' + floor(Nkx4/2) + 1;

parabolic_shift4 = parabolic_shift4(idx1_range, :, :);

x4 = x4(idx1_range, :, :);
y4 = y4(idx1_range, :, :);
z4 = z4(idx1_range, :, :);

%% Read a .cfl file
%--------------------------------------------------------------------------
% parabolic_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path5, 'parabolic_shift');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path5, 'x_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path5, 'y_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path5, 'z_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Remove readout oversampling
[Nkx5,Nky5,Nkz5] = size(parabolic_shift5);

Nx5 = Nkx5 / 2;
Ny5 = Nky5;
Nz5 = Nkz5;

idx1_range = (-floor(Nx5/2):ceil(Nx5/2)-1).' + floor(Nkx5/2) + 1;

parabolic_shift5 = parabolic_shift5(idx1_range, :, :);

x5 = x5(idx1_range, :, :);
y5 = y5(idx1_range, :, :);
z5 = z5(idx1_range, :, :);

%% Read a .cfl file
%--------------------------------------------------------------------------
% parabolic_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'parabolic_shift');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift6 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'x_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x6 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'y_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y6 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'z_shift_flat');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z6 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Remove readout oversampling
[Nkx6,Nky6,Nkz6] = size(parabolic_shift6);

Nx6 = Nkx6 / 2;
Ny6 = Nky6;
Nz6 = Nkz6;

idx1_range = (-floor(Nx6/2):ceil(Nx6/2)-1).' + floor(Nkx6/2) + 1;

parabolic_shift6 = parabolic_shift6(idx1_range, :, :);

x6 = x6(idx1_range, :, :);
y6 = y6(idx1_range, :, :);
z6 = z6(idx1_range, :, :);

%%
Gu = zeros(6, 1, 'double'); % [mT/m]
Gu(1) = 12.13;
Gu(2) = 23.93;
Gu(3) = 14.19;
Gu(4) = 14.57;
Gu(5) = 18.60;
Gu(6) = 15.94;

t_ramp = zeros(6, 1, 'double'); % [msec]
t_ramp(1) = 80 * 1e-3;
t_ramp(2) = 140 * 1e-3;
t_ramp(3) = 90 * 1e-3;
t_ramp(4) = 90 * 1e-3;
t_ramp(5) = 110* 1e-3;
t_ramp(6) = 90 * 1e-3;

t_plateau = zeros(6, 1, 'double'); % [msec]
t_plateau(1) = 900 * 1e-3;
t_plateau(2) = 870 * 1e-3;
t_plateau(3) = 1070 * 1e-3;
t_plateau(4) = 1540 * 1e-3;
t_plateau(5) = 1540 * 1e-3;
t_plateau(6) = 1840 * 1e-3;

t_esp = zeros(6, 1, 'double'); % [msec]
t_esp(1) = ( 80 +  900 +  80) * 1e-3;
t_esp(2) = (140 +  870 + 140) * 1e-3;
t_esp(3) = ( 90 + 1070 +  90) * 1e-3;
t_esp(4) = ( 90 + 1540 +  90) * 1e-3;
t_esp(5) = (110 + 1540 + 110) * 1e-3;
t_esp(6) = ( 90 + 1840 +  90) * 1e-3;

bandwidth = zeros(6, 1, 'double'); % [Hz/Px]
bandwidth(1) = 1002;
bandwidth(2) = 977;
bandwidth(3) = 861;
bandwidth(4) = 610;
bandwidth(5) = 610;
bandwidth(6) = 528;

acceleration_factor = zeros(6, 1, 'double');
acceleration_factor(1) = 1;
acceleration_factor(2) = 2;
acceleration_factor(3) = 2;
acceleration_factor(4) = 1;
acceleration_factor(5) = 1;
acceleration_factor(6) = 1;

interleaves = zeros(6, 1, 'double');
interleaves(1) = 1;
interleaves(2) = 1;
interleaves(3) = 1;
interleaves(4) = 3;
interleaves(5) = 5;
interleaves(6) = 5;

TE = zeros(6, 1, 'double');
TE(1) = 149;
TE(2) = 88;
TE(3) = 80;
TE(4) = 90;
TE(5) = 58;
TE(6) = 80;

%% Display results
close all;

baby_blue = [193 220 243] / 255;
blue      = [0   173 236] / 255;
orange    = [239 173 127] / 255;
green     = [205 235 188] / 255;
yellow    = [253 234 155] / 255;

orange_siemens = [236 99 0] / 255;
green_siemens = [3 153 153] / 255;

red_color = [201 37 31] / 255;
blue_color = [86 120 191] / 255;

color_order{1} = '#1f77b4';
color_order{2} = '#ff7f0e';
color_order{3} = '#2ca02c';
color_order{4} = '#d62728';
color_order{5} = '#9467bd';
color_order{6} = '#8c564b';
color_order{7} = '#e377c2';
color_order{8} = '#7f7f7f';
color_order{9} = '#bcbd22';
color_order{10} = '#17becf';

color_order_rgb = hex2rgb(color_order);

cmap = turbo(256);

figure('Color', 'w', 'Position', [3 5 950 987]);

%--------------------------------------------------------------------------
% Fully sampled 2D SE-EPI (coronal)
%--------------------------------------------------------------------------
s1 = Nx1;
s2 = 1;
s3 = Nz1;

c1 = floor(Nx1/2) + 1;
c2 = floor(Ny1/2) + 1;
c3 = floor(Nz1/2) + 1;

ax1 = subplot(3,3,1);

hold on;
surf(reshape(x1(:,:,s3) * 1e3, [Nx1 Ny1]), reshape(y1(:,:,s3) * 1e3, [Nx1 Ny1]), reshape(z1(:,:,s3) * 1e3, [Nx1 Ny1]), reshape(parabolic_shift1(:,:,s3) * 1e3, [Nx1 Ny1]), 'EdgeColor', 'none'); % axial
surf(reshape(x1(s1,:,:) * 1e3, [Ny1 Nz1]), reshape(y1(s1,:,:) * 1e3, [Ny1 Nz1]), reshape(z1(s1,:,:) * 1e3, [Ny1 Nz1]), reshape(parabolic_shift1(s1,:,:) * 1e3, [Ny1 Nz1]), 'EdgeColor', 'none'); % sagittal
surf(reshape(x1(:,s2,:) * 1e3, [Nx1 Nz1]), reshape(y1(:,s2,:) * 1e3, [Nx1 Nz1]), reshape(z1(:,s2,:) * 1e3, [Nx1 Nz1]), reshape(parabolic_shift1(:,s2,:) * 1e3, [Nx1 Nz1]), 'EdgeColor', 'none'); % coronal
axis image;
set(ax1, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'YDir', 'reverse', 'ZDir', 'reverse');
hXLabel1 = xlabel(ax1, 'x [mm]', 'Interpreter', 'latex', 'FontSize', 12);
hYLabel1 = ylabel(ax1, 'y [mm]', 'Interpreter', 'latex', 'FontSize', 12, 'Position', [178.1953 -67.1093 147.5499]);
hZLabel1 = zlabel(ax1, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
view(ax1, 45, 29);
colormap(ax1, cmap);
hc1 = colorbar('Location', 'eastoutside', 'TickLabelInterpreter', 'latex');
set(hc1, 'Position', [0.3093 0.7568 - 0.1287 - 0.0122 0.0095 0.1185]);
hTitle1_mm = title(hc1, '[mm]', 'Interpreter', 'latex', 'Position', [-4.8656 91.4996 0]);
clim(ax1, [0 10]);

hText1_label = text(ax1, 0, 0, '(A)', 'FontWeight', 'bold', 'FontSize', 12, 'Position', [-118.6616 44.5212 -202.4979]);

ax4 = subplot(3,3,4);
hold on;
plot(squeeze(x1(c1,:,1)) * 1e3, squeeze(parabolic_shift1(c1,:,1)) * 1e3, 'LineWidth', 1, 'Color', color_order{1});
plot(squeeze(x1(c1,:,end)) * 1e3, squeeze(parabolic_shift1(c1,:,end)) * 1e3, 'LineWidth', 1, 'Color', color_order{2});
grid on; grid minor;
xlim(ax4, [-150 150]);
ylim(ax4, [0 10]);
legend(ax4, '$$\delta_{\mathrm{pe}}(x,-15.8)$$', '$$\delta_{\mathrm{pe}}(x,-72.8)$$', 'Interpreter', 'latex', 'Location', 'best');
set(ax4, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax4, 'x [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax4, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
subtitle(ax4, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(1), t_ramp(1), t_plateau(1)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(1), acceleration_factor(1), interleaves(1))}, 'Interpreter', 'latex', 'FontSize', 12);

hText4_label = text(ax4, -145, 10, '(D)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hTitle1_1 = text(ax4, 0, 31.96 - 0.6, {'Figure 2 (coronal)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');
hTitle1_2 = text(ax4, 0, 31.96 - 2.0, {'Parabolic shift $$\delta_{\mathrm{pe}}(x,y)$$'}, 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'normal', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Accelerated 2D SE-EPI (axial)
%--------------------------------------------------------------------------
s1 = Nx2;
s2 = Ny2;
s3 = 1;

c1 = floor(Nx2/2) + 1;
c2 = floor(Ny2/2) + 1;
c3 = floor(Nz2/2) + 1;

ax2 = subplot(3,3,2);
hold on;
surf(reshape(x2(:,:,s3) * 1e3, [Nx2 Ny2]), reshape(y2(:,:,s3) * 1e3, [Nx2 Ny2]), reshape(z2(:,:,s3) * 1e3, [Nx2 Ny2]), reshape(parabolic_shift2(:,:,s3) * 1e3, [Nx2 Ny2]), 'EdgeColor', 'none'); % axial
surf(reshape(x2(s1,:,:) * 1e3, [Ny2 Nz2]), reshape(y2(s1,:,:) * 1e3, [Ny2 Nz2]), reshape(z2(s1,:,:) * 1e3, [Ny2 Nz2]), reshape(parabolic_shift2(s1,:,:) * 1e3, [Ny2 Nz2]), 'EdgeColor', 'none'); % sagittal
surf(reshape(x2(:,s2,:) * 1e3, [Nx2 Nz2]), reshape(y2(:,s2,:) * 1e3, [Nx2 Nz2]), reshape(z2(:,s2,:) * 1e3, [Nx2 Nz2]), reshape(parabolic_shift2(:,s2,:) * 1e3, [Nx2 Nz2]), 'EdgeColor', 'none'); % coronal
axis image;
set(ax2, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'YDir', 'reverse', 'ZDir', 'reverse');
hXLabel2 = xlabel(ax2, 'x [mm]', 'Interpreter', 'latex', 'FontSize', 12);
hYLabel2 = ylabel(ax2, 'y [mm]', 'Interpreter', 'latex', 'FontSize', 12);
hZLabel2 = zlabel(ax2, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
view(ax2, 45, 29);
colormap(ax2, cmap);
hc2 = colorbar('Location', 'eastoutside', 'TickLabelInterpreter', 'latex');
set(hc2, 'Position', [0.6337 0.7568 - 0.1287 - 0.0122 0.0095 0.1185]);
hTitle2_mm = title(hc2, '[mm]', 'Interpreter', 'latex', 'Position', [-4.8656 91.4996 0]);
clim(ax2, [0 20]);

hText2_label = text(ax2, 0, 0, '(B)', 'FontWeight', 'bold', 'FontSize', 12, 'Position', [-168.2678 82.8460 -216.8891]);

ax5 = subplot(3,3,5);
plot(reshape(z2(c1,c2,:) * 1e3, [Nz2 1]), reshape(parabolic_shift2(c1,c2,:) * 1e3, [Nz2 1]), 'LineWidth', 1, 'Color', color_order{1});
grid on; grid minor;
set(ax5, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax5, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax5, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
xlim(ax5, [-150 150]);
ylim(ax5, [0 20]);
subtitle(ax5, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(2), t_ramp(2), t_plateau(2)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(2), acceleration_factor(2), interleaves(2))}, 'Interpreter', 'latex', 'FontSize', 12);

hText5_label = text(ax5, -145, 20, '(E)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hTitle2_1 = text(ax5, 0, (31.96 - 0.6) * 2, {'Figure 3 (axial)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');
hTitle2_2 = text(ax5, 0, (31.96 - 2.0) * 2, {'Parabolic shift $$\delta_{\mathrm{pe}}(z)$$'}, 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'normal', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Title
%--------------------------------------------------------------------------
hText1 = text(ax5, 0, (+42.30) * 2, sprintf('Parabolic shift along the PE direction induced by concomitant fields'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

hText2 = text(ax5, 0, (+36.26 - 0.4) * 2, {'$$\lambda(x,y,z) = \frac{\gamma}{2B_{0}} \left\{ \left(G_x z - \frac{G_z x}{2} \right)^2 + \left(G_y z - \frac{G_z y}{2} \right)^2 \right\}  \left(\frac{2}{3}t_{\mathrm{ramp}} + t_{\mathrm{plateau}} \right)$$ [rad]', ...
                            '$$\delta_{\mathrm{pe}}(x,y,z) = \frac{\lambda(x,y,z)}{\mathrm{Number \; of \; in \textendash plane \; shots} \times R \times dk_{\mathrm{pe}}} $$, where $$dk_{\mathrm{pe}} = \frac{2 \pi}{FOV_{\mathrm{pe}}} $$'}, ...
                            'Color', color_order{1}, 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

hText3 = text(ax5, 0, (+33.60-0.85) * 2, {'$$G_x/G_y/G_z$$: readout gradient amplitude $$G_{\mathrm{u}}$$ in the physical axes [mT/m], $$R$$: acceleration factor', ...
                            '$$t_{\mathrm{ramp}}/t_{\mathrm{plateau}}$$: ramp-up time/plateau time of a trapezoidal readout gradient lobe [ms]'}, ...
                           'Color', 'k', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Create line (bottom)
annotation(gcf, 'line', [0.0700 0.9600], [0.8170 0.8170], 'LineWidth', 2);

% Create line (top)
annotation(gcf, 'line', [0.0700 0.9600], [0.8271 0.8271] + 0.138, 'LineWidth', 2);

%--------------------------------------------------------------------------
% Diffusion-weighted 2D SE-EPI (axial)
%--------------------------------------------------------------------------
s1 = Nx3;
s2 = Ny3;
s3 = 1;

c1 = floor(Nx3/2) + 1;
c2 = floor(Ny3/2) + 1;
c3 = floor(Nz3/2) + 1;

ax3 = subplot(3,3,3);
hold on;
surf(reshape(x3(:,:,s3) * 1e3, [Nx3 Ny3]), reshape(y3(:,:,s3) * 1e3, [Nx3 Ny3]), reshape(z3(:,:,s3) * 1e3, [Nx3 Ny3]), reshape(parabolic_shift3(:,:,s3) * 1e3, [Nx3 Ny3]), 'EdgeColor', 'none'); % axial
surf(reshape(x3(s1,:,:) * 1e3, [Ny3 Nz3]), reshape(y3(s1,:,:) * 1e3, [Ny3 Nz3]), reshape(z3(s1,:,:) * 1e3, [Ny3 Nz3]), reshape(parabolic_shift3(s1,:,:) * 1e3, [Ny3 Nz3]), 'EdgeColor', 'none'); % sagittal
surf(reshape(x3(:,s2,:) * 1e3, [Nx3 Nz3]), reshape(y3(:,s2,:) * 1e3, [Nx3 Nz3]), reshape(z3(:,s2,:) * 1e3, [Nx3 Nz3]), reshape(parabolic_shift3(:,s2,:) * 1e3, [Nx3 Nz3]), 'EdgeColor', 'none'); % coronal
axis image;
set(ax3, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'YDir', 'reverse', 'ZDir', 'reverse');
hXLabel3 = xlabel(ax3, 'x [mm]', 'Interpreter', 'latex', 'FontSize', 12);
hYLabel3 = ylabel(ax3, 'y [mm]', 'Interpreter', 'latex', 'FontSize', 12);
hZLabel3 = zlabel(ax3, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12, 'Position', [-168.0099 159.2428 3.9288]);
view(ax3, 45, 29);
colormap(ax3, cmap);
hc3 = colorbar('Location', 'eastoutside', 'TickLabelInterpreter', 'latex');
set(hc3, 'Position', [0.9326 0.7568 - 0.1287 - 0.0122 0.0095 0.1185]);
hTitle3_mm = title(hc3, '[mm]', 'Interpreter', 'latex', 'Position', [-4.8656 91.4996 0]);

hText3_label = text(ax3, 0, 0, '(C)', 'FontWeight', 'bold', 'FontSize', 12, 'Position', [-183.6170 86.3035 -232.0588]);

ax6 = subplot(3,3,6);
hold on;
plot(reshape(z3(1,1,:) * 1e3, [Nz3 1]), reshape(parabolic_shift3(1,1,:) * 1e3, [Nz3 1]), 'LineWidth', 1, 'Color', color_order{1});
plot(reshape(z3(c1,c2,:) * 1e3, [Nz3 1]), reshape(parabolic_shift3(c1,c2,:) * 1e3, [Nz3 1]), 'LineWidth', 1, 'Color', color_order{2});
grid on; grid minor;
xlim(ax6, [-150 150]);
ylim(ax6, [0 10]);
legend(ax6, '$$\delta_{\mathrm{pe}}(135,109,z)$$', '$$\delta_{\mathrm{pe}}(1,-11,z)$$', 'Interpreter', 'latex', 'Location', 'best');
set(ax6, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax6, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax6, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
%title(ax6, {'Figures 4 & 5 (tilted axial)'}, 'Interpreter', 'tex', 'Color', 'k');
subtitle(ax6, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(3), t_ramp(3), t_plateau(3)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(3), acceleration_factor(3), interleaves(3))}, 'Interpreter', 'latex');

hText6_label = text(ax6, -145, 10, '(F)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hTitle3 = text(ax6, 0, 31.96 - 0.6, {'Figures 4 & 5 (tilted axial)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');
hTitle4 = text(ax6, 0, 31.96 - 2.0, {'Parabolic shift $$\delta_{\mathrm{pe}}(x,y,z)$$'}, 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'normal', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Long-ETL 3D GRE-EPI (axial)
%--------------------------------------------------------------------------
c1 = floor(Nx4/2) + 1;
c2 = floor(Ny4/2) + 1;
c3 = floor(Nz4/2) + 1;

ax7 = subplot(3,3,7);
plot(reshape(z4(c1,c2,:) * 1e3, [Nz4 1]), reshape(parabolic_shift4(c1,c2,:) * 1e3, [Nz4 1]), 'LineWidth', 1, 'Color', color_order{1});
grid on; grid minor;
set(ax7, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax7, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax7, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
xlim(ax7, [-150 150]);
ylim(ax7, [0 10]);
title(ax7, {'Figure 6 (axial)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12);
subtitle(ax7, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(4), t_ramp(4), t_plateau(4)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(4), acceleration_factor(4), interleaves(4))}, 'Interpreter', 'latex');

hText7_label = text(ax7, -145, 10, '(G)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% High-resolution 3D GRE-EPI (axial)
%--------------------------------------------------------------------------
c1 = floor(Nx5/2) + 1;
c2 = floor(Ny5/2) + 1;
c3 = floor(Nz5/2) + 1;

ax8 = subplot(3,3,8);
plot(reshape(z5(c1,c2,:) * 1e3, [Nz5 1]), reshape(parabolic_shift5(c1,c2,:) * 1e3, [Nz5 1]), 'LineWidth', 1, 'Color', color_order{1});
grid on; grid minor;
set(ax8, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax8, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax8, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
xlim(ax8, [-150 150]);
ylim(ax8, [0 10]);
title(ax8, {'Figure 7 (axial)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12);
subtitle(ax8, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(5), t_ramp(5), t_plateau(5)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(5), acceleration_factor(5), interleaves(5))}, 'Interpreter', 'latex');

hText8_label = text(ax8, -145, 10, '(H)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% High-resolution 3D GRE-EPI in-vivo data (axial)
%--------------------------------------------------------------------------
c1 = floor(Nx6/2) + 1;
c2 = floor(Ny6/2) + 1;
c3 = floor(Nz6/2) + 1;

ax9 = subplot(3,3,9);
plot(reshape(z6(c1,c2,:) * 1e3, [Nz6 1]), reshape(parabolic_shift6(c1,c2,:) * 1e3, [Nz6 1]), 'LineWidth', 1, 'Color', color_order{1});
grid on; grid minor;
set(ax9, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax9, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax9, '[mm]', 'Interpreter', 'latex', 'FontSize', 12);
xlim(ax9, [-150 150]);
ylim(ax9, [0 10]);
title(ax9, {'Figure 8 (axial)'}, 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 12);
subtitle(ax9, {sprintf('$$G_{\\mathrm{u}}/t_{\\mathrm{ramp}}/t_{\\mathrm{plateau}}$$ = %5.2f/%4.2f/%4.2f', Gu(6), t_ramp(6), t_plateau(6)), ...
               sprintf('BW/R/shots = %3d/%d/%d', bandwidth(6), acceleration_factor(6), interleaves(6))}, 'Interpreter', 'latex');

hText9_label = text(ax9, -145, 10, '(I)', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

set(ax1, 'Position', [0.1300-0.02 0.5524 0.2134 0.2157]);
set(ax2, 'Position', [0.4108      0.5524 0.2134 0.2157]);
set(ax3, 'Position', [0.6916+0.02 0.5524 0.2134 0.2157]);

set(ax4, 'Position', [0.1300-0.02 0.3117 0.2326 0.1557]);
set(ax5, 'Position', [0.4108      0.3117 0.2326 0.1557]);
set(ax6, 'Position', [0.6916+0.02 0.3117 0.2326 0.1557]);

set(ax7, 'Position', [0.1300-0.02 0.0472 0.2326 0.1557]);
set(ax8, 'Position', [0.4108      0.0472 0.2326 0.1557]);
set(ax9, 'Position', [0.6916+0.02 0.0472 0.2326 0.1557]);

export_fig(sprintf('figure9'), '-r300', '-tif', '-c[10, 100, 0, 190]'); % [top,right,bottom,left]
