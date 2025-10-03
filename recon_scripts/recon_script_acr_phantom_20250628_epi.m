% recon_script_acr_phantom_20250628_epi.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/29/2025, Last modified: 06/29/2025

%% Clean slate
close all; clear all; clc;

%% Set source directories
if ispc
    package_path = 'D:\cartesian_maxgirf_epi_3d';
    ismrmrd_path = 'D:\ismrmrd';
    grad_file_path = 'D:\cartesian_maxgirf_epi_3d\GradientCoils';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    grad_file_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/GradientCoils';
end

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define a list of .json files
if ispc
    %----------------------------------------------------------------------
    % ETL = 61, offset = 40 mm
    %----------------------------------------------------------------------
    json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{2} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';

    %----------------------------------------------------------------------
    % ETL = 101, offset = 40 mm
    %----------------------------------------------------------------------
    json_files{3} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00288_FID18548_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{4} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00288_FID18548_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH40mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';

    %----------------------------------------------------------------------
    % ETL = 61, offset = 0 mm
    %----------------------------------------------------------------------
    json_files{5} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00279_FID18539_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{6} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00279_FID18539_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';

    %----------------------------------------------------------------------
    % ETL = 61, offset = 0 mm, noise
    %----------------------------------------------------------------------
    json_files{7} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00296_FID18556_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_noise_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{8} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00296_FID18556_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH0mm_noise_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';

    %----------------------------------------------------------------------
    % ETL = 101, offset = 0 mm
    %----------------------------------------------------------------------
    json_files{9}  = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00277_FID18537_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{10} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00277_FID18537_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH0mm_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';

    %----------------------------------------------------------------------
    % ETL = 101, offset = 0 mm, noise
    %----------------------------------------------------------------------
    json_files{11} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00294_FID18554_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH0mm_noise_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{12} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250628\json_files\meas_MID00294_FID18554_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH0mm_noise_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';
else
    %----------------------------------------------------------------------
    % ETL = 61, offset = 40 mm
    %----------------------------------------------------------------------
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250628/json_files/meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_cutoff99_999_gridding1_phc1_cfc1_sfc0_gnc0_topup0_server.json';
    json_files{2} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250628/json_files/meas_MID00290_FID18550_ep3d_tra_AP_highres_TR226_TE88_etl61_esp1_66_0_8mm_offsetH40mm_cutoff99_999_gridding1_phc1_cfc1_sfc0_gnc1_topup0_server.json';

    %----------------------------------------------------------------------
    % ETL = 101, offset = 40 mm
    %----------------------------------------------------------------------
    json_files{3} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250628/json_files/meas_MID00288_FID18548_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH40mm_cutoff99_999_gridding1_phc1_cfc1_sfc0_gnc0_topup0_server.json';
    json_files{4} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250628/json_files/meas_MID00288_FID18548_ep3d_tra_AP_highres_TR226_TE88_etl101_esp1_66_0_8mm_offsetH40mm_cutoff99_999_gridding1_phc1_cfc1_sfc0_gnc1_topup0_server.json';
end

%% Define a list of a .json file (gre_field_mapping)
if ispc
    fieldmap_json_file = '';
else
    fieldmap_json_file = '';
end

%% Calculate the number of json files
nr_json_files = length(json_files);

%% Process per json file
for json_number = 1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_cartesian_maxgirf_3d_calculate_voxel_coordinates;

    %% Calculate an interpolated fieldmap
    demo_interpolate_fieldmap_gre_field_mapping;

    %% Prepare "imaging" k-space data
    demo_cartesian_maxgirf_3d_prepare_ksp_imaging;

    %% Estimate CSMs
    demo_cartesian_maxgirf_3d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_cartesian_maxgirf_3d_recon;
end
