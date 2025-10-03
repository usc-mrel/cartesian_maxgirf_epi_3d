% demo_figure10.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 03/04/2025, Last modified: 03/11/2025

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
% ETL = 61, offset = 0 mm, Gridding/PHC
%--------------------------------------------------------------------------
output_path1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\meas_MID00279_FID18539_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% ETL = 61, offset = 0 mm, Gridding/PHC/CFC
%--------------------------------------------------------------------------
output_path2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\meas_MID00279_FID18539_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% ETL = 61, offset = 40 mm, Gridding/PHC
%--------------------------------------------------------------------------
output_path3 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% ETL = 61, offset = 40 mm, Gridding/PHC/CFC
%--------------------------------------------------------------------------
output_path4 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0';

%--------------------------------------------------------------------------
% ETL = 61, offset = 40 mm, Gridding/PHC/CFC/GNC
%--------------------------------------------------------------------------
output_path5 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc1_sfc0_gnc1_topup0_cutoff99_999';

%--------------------------------------------------------------------------
% ETL = 61, offset = 0 mm, Gridding/PHC
%--------------------------------------------------------------------------
dicom_path1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\dicom\S12_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1.66_0.8mm_offsetH0mm'; % ND

%--------------------------------------------------------------------------
% ETL = 61, offset = 40 mm, Gridding/PHC
%--------------------------------------------------------------------------
dicom_path2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\dicom\S14_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1.66_0.8mm_offsetH40mm'; % ND

%% Get directory information
dir_info = dir(fullfile(dicom_path1, '*IMA'));
nr_files = length(dir_info);

%% Get a DICOM header
dicom_info = dicominfo(fullfile(dir_info(1).folder, dir_info(1).name));

%% Parse the DICOM header
%--------------------------------------------------------------------------
% Rows Attribute
% Number of rows in the image
%--------------------------------------------------------------------------
N1 = double(dicom_info.Rows);

%--------------------------------------------------------------------------
% Columns  Attribute
% Number of columns in the image
%--------------------------------------------------------------------------
N2 = double(dicom_info.Columns);

%% Read DICOM files
img1_dicom = zeros(N1, N2, nr_files, 'single');
x_dicom1 = zeros(N1, N2, nr_files, 'single');
y_dicom1 = zeros(N1, N2, nr_files, 'single');
z_dicom1 = zeros(N1, N2, nr_files, 'single');

%% Get directory information
dir_info = dir(fullfile(dicom_path1, '*IMA'));
nr_files = length(dir_info);

for idx = 1:nr_files

    %% Set the fullpath of a DICOM file
    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('%s:(i=%3d/%3d) Reading a .dicom file: %s... ', datetime, idx, nr_files, dicom_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Get a DICOM header
    dicom_info = dicominfo(dicom_file);

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img1_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x_dicom1(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y_dicom1(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z_dicom1(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Read DICOM files
img2_dicom = zeros(N1, N2, nr_files, 'single');
x_dicom2 = zeros(N1, N2, nr_files, 'single');
y_dicom2 = zeros(N1, N2, nr_files, 'single');
z_dicom2 = zeros(N1, N2, nr_files, 'single');

%% Get directory information
dir_info = dir(fullfile(dicom_path2, '*IMA'));
nr_files = length(dir_info);

for idx = 1:nr_files

    %% Set the fullpath of a DICOM file
    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('%s:(i=%3d/%3d) Reading a .dicom file: %s... ', datetime, idx, nr_files, dicom_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Get a DICOM header
    dicom_info = dicominfo(dicom_file);

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img2_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x_dicom2(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y_dicom2(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z_dicom2(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Set variables
Nkx = 512;
Nky = 306;
Nkz = 176;

[Nx,Ny,Nz,nr_averages] = size(img1_dicom);

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, 'img_maxgirf_flat_gridding1_phc1_cfc0_sfc0_gnc0_topup0_i6_l0.00');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path2, 'img_maxgirf_flat_gridding1_phc1_cfc1_sfc0_gnc0_topup0_i6_l0.00');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, 'img_maxgirf_flat_gridding1_phc1_cfc0_sfc0_gnc0_topup0_i6_l0.00');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img3 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'img_maxgirf_flat_gridding1_phc1_cfc1_sfc0_gnc0_topup0_i6_l0.00');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path5, 'img_maxgirf_flat_gridding1_phc1_cfc1_sfc0_gnc1_topup0_i6_l0.00');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
slice_type = 'flat';

%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path3, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Adjust the image size of a custom reconstruction
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
idx3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(Nkz/2) + 1;

img1 = img1(idx1_range, idx2_range, idx3_range, :);
img2 = img2(idx1_range, idx2_range, idx3_range, :);
img3 = img3(idx1_range, idx2_range, idx3_range, :);
img4 = img4(idx1_range, idx2_range, idx3_range, :);
img5 = img5(idx1_range, idx2_range, idx3_range, :);

x1 = x1(idx1_range, idx2_range, idx3_range);
y1 = y1(idx1_range, idx2_range, idx3_range);
z1 = z1(idx1_range, idx2_range, idx3_range);

x2 = x2(idx1_range, idx2_range, idx3_range);
y2 = y2(idx1_range, idx2_range, idx3_range);
z2 = z2(idx1_range, idx2_range, idx3_range);

%% Flip variables
if read_sign1 < 0
    img1 = flip(img1,1);
    img2 = flip(img2,1);
    img3 = flip(img3,1);
    img4 = flip(img4,1);
    img5 = flip(img5,1);

    x1 = flip(x1,1);
    y1 = flip(y1,1);
    z1 = flip(z1,1);

    x2 = flip(x2,1);
    y2 = flip(y2,1);
    z2 = flip(z2,1);
end

if phase_sign1 < 0
    img1 = flip(img1,2);
    img2 = flip(img2,2);
    img3 = flip(img3,2);
    img4 = flip(img4,2);
    img5 = flip(img5,2);

    x1 = flip(x1,2);
    y1 = flip(y1,2);
    z1 = flip(z1,2);

    x2 = flip(x2,2);
    y2 = flip(y2,2);
    z2 = flip(z2,2);
end

%% Scale images
c1 = floor(Nx/2) + 1;
c2 = floor(Ny/2) + 1;
c3 = floor(Nz/2) + 1;

r1 = 140;
c1 = 25;
s1 = c3;

r2 = 140;
c2 = 28;
s2 = c3;

scale_factor_custom1 = abs(img1(r1,c1,s1)); % ETL = 61, offset = 0 mm, Gridding/PHC
scale_factor_custom2 = abs(img3(r2,c2,s2)); % ETL = 61, offset = 40 mm, Gridding/PHC

scale_factor_dicom1 = abs(img1_dicom(r1,c1,s1)); % ETL = 61, offset = 0 mm, Gridding/PHC
scale_factor_dicom2 = abs(img2_dicom(r2,c2,s2)); % ETL = 61, offset = 40 mm, Gridding/PHC

img1_scaled = img1 / scale_factor_custom1 * scale_factor_dicom1;
img2_scaled = img2 / scale_factor_custom1 * scale_factor_dicom1;

img3_scaled = img3 / scale_factor_custom2 * scale_factor_dicom2;
img4_scaled = img4 / scale_factor_custom2 * scale_factor_dicom2;
img5_scaled = img5 / scale_factor_custom2 * scale_factor_dicom2;

%% Calculate a spatial shift along the readout direction [m]
B0 = 0.55;
Gu = 19.94; % gradient amplitude [mT/m]

z_offsetH0mm = reshape(z1(1,1,:), [Nz 1]);
z_offsetH40mm = reshape(z2(1,1,:), [Nz 1]);

% [mT/m] * [T/1e3mT] / [T] * [m]^2 => [m]
delta_ro_offsetH0mm  = 1 / (2 * B0) * (Gu * 1e-3) * z_offsetH0mm.^2;
delta_ro_offsetH40mm = 1 / (2 * B0) * (Gu * 1e-3) * z_offsetH40mm.^2;

%% Display images
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

c1 = floor(Nx/2) + 1;
c2 = floor(Ny/2) + 1;
c3 = floor(Nz/2) + 1;

cmap = flip(brewermap([], "RdBu"), 1);

FontSize = 12;

climits = [0 1700];

axial_slice_number1 = 30;
axial_slice_number2 = c3;

s1 = 30;
s2 = 50;
s3 = 70;
s4 = 90;

img1_dicom_montage = reshape(flip(rot90(img1_dicom(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);
img2_dicom_montage = reshape(flip(rot90(img2_dicom(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);

img1_montage = reshape(flip(rot90(img1_scaled(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);
img2_montage = reshape(flip(rot90(img2_scaled(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);
img3_montage = reshape(flip(rot90(img3_scaled(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);
img4_montage = reshape(flip(rot90(img4_scaled(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);
img5_montage = reshape(flip(rot90(img5_scaled(:,:, [s1 s2 s3 s4]),-1),2), [Nx Ny * 4]);

figure('Color', 'w', 'Position', [3 5 950 987]);
%--------------------------------------------------------------------------
% DICOM: ETL = 61, offset = 0 mm, Gridding/PHC
%--------------------------------------------------------------------------
ax1 = subplot(5,1,1);
imagesc(abs(img1_dicom_montage));
axis image off;
colormap(gray(256));
clim(climits);
text(ax1, -30, Ny / 2, 'Traditional recon', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, 0, Ny / 2, 'Gridding/PHC', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, 2, 0, '(A)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

text(ax1, Nx / 2 + Nx * 0, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(2,1), color_order_rgb(2,2), color_order_rgb(2,3), s1)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 0, 0  , {sprintf('z = %4.2f mm', z1(1,1,s1) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH0mm(s1) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 1, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(3,1), color_order_rgb(3,2), color_order_rgb(3,3), s2)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 1, 0  , {sprintf('z = %4.2f mm', z1(1,1,s2) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH0mm(s2) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 2, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(4,1), color_order_rgb(4,2), color_order_rgb(4,3), s3)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 2, 0  , {sprintf('z = %4.2f mm', z1(1,1,s3) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH0mm(s3) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 3, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(5,1), color_order_rgb(5,2), color_order_rgb(5,3), s4)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, Nx / 2 + Nx * 3, 0  , {sprintf('z = %4.2f mm', z1(1,1,s4) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH0mm(s4) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax1, 2 * Nx, -105, {'\underline{\textbf{Table position = 0 mm}}'}, 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%--------------------------------------------------------------------------
% title
%--------------------------------------------------------------------------
text(ax1, 3 * Nx, -205, sprintf('Mitigation of a slice-dependent Nyquist ghost induced by concomitant fields'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax1, 3 * Nx, -155, sprintf('3D GRE-EPI: axial, 0.78 x 0.78 x 0.90 mm^3, 176 slices, R = 1, no PF, ETL = 61, 1 NSA'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: ETL = 61, offset = 0 mm, Gridding/PHC/CFC
%--------------------------------------------------------------------------
ax2 = subplot(5,1,2);
imagesc(abs(img2_montage));
axis image off;
colormap(gray(256));
clim(climits);
text(ax2, -30, Ny / 2, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, 0, Ny / 2, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, 2, 0, '(B)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% DICOM: ETL = 61, offset = 40 mm, Gridding/PHC
%--------------------------------------------------------------------------
ax3 = subplot(5,1,3);
imagesc(abs(img2_dicom_montage));
axis image off;
colormap(gray(256));
clim(climits);
text(ax3, -30, Ny / 2, 'Traditional recon', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, 0, Ny / 2, 'Gridding/PHC', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, 2, 0, '(C)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

text(ax3, Nx / 2 + Nx * 0, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(2,1), color_order_rgb(2,2), color_order_rgb(2,3), s1)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 0, 0  , {sprintf('z = %4.2f mm', z2(1,1,s1) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH40mm(s1) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 1, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(3,1), color_order_rgb(3,2), color_order_rgb(3,3), s2)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 1, 0  , {sprintf('z = %4.2f mm', z2(1,1,s2) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH40mm(s2) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 2, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(4,1), color_order_rgb(4,2), color_order_rgb(4,3), s3)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 2, 0  , {sprintf('z = %4.2f mm', z2(1,1,s3) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH40mm(s3) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 3, -70, {sprintf('{\\color[rgb]{%f %f %f}Slice %d}', color_order_rgb(5,1), color_order_rgb(5,2), color_order_rgb(5,3), s4)}, 'Interpreter', 'tex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, Nx / 2 + Nx * 3, 0  , {sprintf('z = %4.2f mm', z2(1,1,s4) * 1e3), sprintf('$$\\delta_{\\mathrm{ro}}(z)$$ = %4.3f mm', delta_ro_offsetH40mm(s4) * 1e3)}, 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax3, 2 * Nx, -105, {'\underline{\textbf{Table position = 40 mm}}'}, 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: ETL = 61, offset = 40 mm, Gridding/PHC/CFC
%--------------------------------------------------------------------------
ax4 = subplot(5,1,4);
imagesc(abs(img4_montage));
axis image off;
colormap(gray(256));
clim(climits);
text(ax4, -30, Ny, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, 0, Ny / 2, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, 2, 0, '(D)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: ETL = 61, offset = 40 mm, Gridding/PHC/CFC/GNC
%--------------------------------------------------------------------------
ax5 = subplot(5,1,5);
imagesc(abs(img5_montage));
axis image off;
colormap(gray(256));
clim(climits);
text(ax5, 0, Ny / 2, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Color', 'k', 'Interpreter', 'tex', 'Rotation', 90, 'FontSize', 12, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax5, 2, 0, '(E)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Spatial shift along the RO direction (table position = 0 mm)
%--------------------------------------------------------------------------
ax6 = axes;
hold on;
plot(z_offsetH0mm * 1e3, delta_ro_offsetH0mm * 1e3, 'LineWidth', 1);
plot(z_offsetH0mm(s1) * 1e3, delta_ro_offsetH0mm(s1) * 1e3, '.', 'Color', color_order{2}, 'MarkerSize', 20);
plot(z_offsetH0mm(s2) * 1e3, delta_ro_offsetH0mm(s2) * 1e3, '.', 'Color', color_order{3}, 'MarkerSize', 20);
plot(z_offsetH0mm(s3) * 1e3, delta_ro_offsetH0mm(s3) * 1e3, '.', 'Color', color_order{4}, 'MarkerSize', 20);
plot(z_offsetH0mm(s4) * 1e3, delta_ro_offsetH0mm(s4) * 1e3, '.', 'Color', color_order{5}, 'MarkerSize', 20);
grid on;
set(ax6, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 10);
legend(ax6, '', sprintf('Slice %d', s1), sprintf('Slice %d', s2), sprintf('Slice %d', s3), sprintf('Slice %d', s4), 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11);
hTitle6 = title(ax6, '\textbf{(F)} Spatial shift $$\delta_{\mathrm{ro}}(z)$$ along the RO direction', 'Interpreter', 'latex');
subtitle(ax6, {'For axial EPI: $$\delta_{\mathrm{ro}}(z) = \frac{1}{2 B_0} G_{\mathrm{ro}} z^2$$ [mm]', sprintf('$$B_0$$ = 0.55 T, $$G_{\\mathrm{ro}}$$ = %5.2f mT/m', Gu)}, 'Interpreter', 'latex', 'FontSize', 11);
xlabel(ax6, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);

%--------------------------------------------------------------------------
% Spatial shift along the RO direction (table position = 40 mm)
%--------------------------------------------------------------------------
ax7 = axes;
hold on;
plot(z_offsetH40mm * 1e3, delta_ro_offsetH40mm * 1e3, 'LineWidth', 1);
plot(z_offsetH40mm(s1) * 1e3, delta_ro_offsetH40mm(s1) * 1e3, '.', 'Color', color_order{2}, 'MarkerSize', 20);
plot(z_offsetH40mm(s2) * 1e3, delta_ro_offsetH40mm(s2) * 1e3, '.', 'Color', color_order{3}, 'MarkerSize', 20);
plot(z_offsetH40mm(s3) * 1e3, delta_ro_offsetH40mm(s3) * 1e3, '.', 'Color', color_order{4}, 'MarkerSize', 20);
plot(z_offsetH40mm(s4) * 1e3, delta_ro_offsetH40mm(s4) * 1e3, '.', 'Color', color_order{5}, 'MarkerSize', 20);
grid on;
set(ax7, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 10);
legend(ax7, '', sprintf('Slice %d', s1), sprintf('Slice %d', s2), sprintf('Slice %d', s3), sprintf('Slice %d', s4), 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11);
hTitle7 = title(ax7, '\textbf{(G)} Spatial shift $$\delta_{\mathrm{ro}}(z)$$ along the RO direction', 'Interpreter', 'latex');
subtitle(ax7, {'For axial EPI: $$\delta_{\mathrm{ro}}(z) = \frac{1}{2 B_0} G_{\mathrm{ro}} z^2$$ [mm]', sprintf('$$B_0$$ = 0.55 T, $$G_{\\mathrm{ro}}$$ = %5.2f mT/m', Gu)}, 'Interpreter', 'latex', 'FontSize', 11);
xlabel(ax7, 'z [mm]', 'Interpreter', 'latex', 'FontSize', 12);

set(ax1, 'Position', [0.0479 0.5029 + 0.1499 - 0.035       0.6153 0.3001]);
set(ax2, 'Position', [0.0479 0.3530 + 0.1499 - 0.035       0.6153 0.3001]);
set(ax3, 'Position', [0.0479 0.1651-0.022 + 0.1499 - 0.06 0.6153 0.3001]);
set(ax4, 'Position', [0.0479 0.0152-0.022 + 0.1499 - 0.06 0.6153 0.3001]);
set(ax5, 'Position', [0.0479 0.0152-0.022 - 0.06          0.6153 0.3001]);

set(ax6, 'Position', [0.7053-0.005 0.5765 0.2724 0.1879]);
set(ax7, 'Position', [0.7053-0.005 0.1116 0.2724 0.1879]);

set(hTitle6, 'Position', [-5.4170 0.1556 0]);
set(hTitle7, 'Position', [43.8262 0.3885 0]);

export_fig(sprintf('figure10'), '-r300', '-tif', '-c[20, 20, 20, 20]'); % [top,right,bottom,left]
