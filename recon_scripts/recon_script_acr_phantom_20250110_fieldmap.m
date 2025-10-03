% batch_acr_phantom_20250110_fieldmap.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/30/2025, Last modified: 07/18/2025

%% Clean slate
close all; clear all; clc;

%% Set source directories
if ispc
    package_path   = 'D:\cartesian_maxgirf_epi_3d';
    ismrmrd_path   = 'D:\ismrmrd';
    grad_file_path = 'D:\cartesian_maxgirf_epi_3d\GradientCoils';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    grad_file_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/GradientCoils';
end

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define a .json file
if ispc
    json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\json_files\meas_MID00262_FID07210_gre_field_mapping_gnc1.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/json_files/meas_MID00262_FID07210_gre_field_mapping_gnc1_server.json';
end

%% Calculate the number of json files
nr_json_files = length(json_files);

%% Perform image reconstruction
for json_number = 1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_type1_nufft_sense_2d_calculate_voxel_coordinates;

    %% Calculate voxel coordinates
    demo_type1_nufft_sense_2d_prepare_ksp_calibration;

    %% Prepare k-space data
    demo_type1_nufft_sense_2d_prepare_ksp_imaging;

    %% Estimate CSMs
    demo_type1_nufft_sense_2d_estimate_csm;

    %% Perform type-1 NUFFT based CG-SENSE reconstruction
    demo_type1_nufft_sense_2d_recon;
end
