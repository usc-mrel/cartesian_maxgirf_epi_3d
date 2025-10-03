% demo_convert_dicom_to_nifti.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/30/2025, Last modified: 04/30/2025

% From FSL FUGUE:
% "Fieldmap data from a SIEMENS scanner takes the form of one phase difference
% image and two magnitude images (one for each echo time)."

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Define the full path of FSL
fsl_path = '/home/image/fsl';

%% Define the full path of dcm2niix
dcm2niix_path = 'D:\dcm2niix_win';

%% Define the full path of data directory
magnitude_path = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S9_gre_field_mapping';
phase_path     = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S10_gre_field_mapping';

%% Define the full path of an output directory
output_path = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\fsl_fieldmap';

%% Make an output directory
mkdir(output_path);

%% Set up FSL commands
%--------------------------------------------------------------------------
% Define an FSL command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
fsl_command = sprintf('%s . %s/etc/fslconf/fsl.sh; %s/bin', command_prefix, fsl_path, fsl_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path
%--------------------------------------------------------------------------
if ispc
    fsl_output_path = strrep(output_path, '\', '/');
    fsl_output_path = sprintf('/mnt/%s/%s/', lower(fsl_output_path(1)), fsl_output_path(4:end));
else
    fsl_output_path = sprintf('%s/', output_path);
end

%% Convert DICOM images to a NIfTI file (magnitude)
command = sprintf('%s%sdcm2niix -f magnitude_gre_field_mapping -o %s -z y %s', dcm2niix_path, filesep, output_path, magnitude_path);
tstart = tic; fprintf('%s:[dcm2niix] Converting DICOM images to a NIfTI file:\n%s\n', datetime, command);
[status_dcm2niix1,result_dcm2niix1] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Convert DICOM images to a NIfTI file (phase)
command = sprintf('%s%sdcm2niix -f phase_gre_field_mapping -o %s -z y %s', dcm2niix_path, filesep, output_path, phase_path);
tstart = tic; fprintf('%s:[dcm2niix] Converting DICOM images to a NIfTI file:\n%s\n', datetime, command);
[status_dcm2niix2,result_dcm2niix2] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Get two echo times [msec]
TE = zeros(2, 1, 'double'); % [msec]

for echo_number = 1:2
    %----------------------------------------------------------------------
    % Define the full path of a .json file
    %----------------------------------------------------------------------
    json_file = sprintf('%s%smagnitude_gre_field_mapping_e%d.json', output_path, filesep, echo_number);

    %----------------------------------------------------------------------
    % Read a .json file
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
    fid = fopen(json_file);
    json_txt = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
    json = jsondecode(json_txt);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Get echo time [msec]
    %----------------------------------------------------------------------
    TE(echo_number) = json.EchoTime * 1e3;
end

%% Calculate deltaTE (the echo time difference of the fieldmap sequence)
deltaTE = TE(2) - TE(1);

%% Calculate a fieldmap [rad]
%--------------------------------------------------------------------------
% Usage: fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)> [--nocheck]
% 
%   Prepares a fieldmap suitable for FEAT from SIEMENS or GEHC data - saves output in rad/s format
%   <scanner> must be SIEMENS or GEHC_FIELDMAPHZ
%   <phase_image> should be the phase difference for SIEMENS and the fieldmap in HERTZ for GEHC_FIELDMAPHZ
%   <magnitude image> should be Brain Extracted (with BET or otherwise)
%   <deltaTE> is the echo time difference of the fieldmap sequence - find this out form the operator (defaults are *usually* 2.46ms on SIEMENS)
%             (defaults are *usually* 2.304ms for GEHC 2D-B0MAP at 3.0T and 2.272 ms GEHC 3D B0MAP at 3.0T)
%   --nocheck supresses automatic sanity checking of image size/range/dimensions
% 
%    e.g. fsl_prepare_fieldmap SIEMENS images_3_gre_field_mapping images_4_gre_field_mapping fmap_rads 2.65
%    e.g. fsl_prepare_fieldmap GEHC_FIELDMAPHZ 3dB0map_fieldmaphz mag_3dB0map fmap_rads 2.272
%--------------------------------------------------------------------------
scanner        = 'SIEMENS';
phase_file     = strcat(fsl_output_path, 'phase_gre_field_mapping_e2_ph');
magnitude_file = strcat(fsl_output_path, 'magnitude_gre_field_mapping_e1');
out_file       = strcat(fsl_output_path, 'fieldmap_rads');

command = sprintf('%s/fsl_prepare_fieldmap %s %s %s %s %f --nocheck', fsl_command, scanner, phase_file, magnitude_file, out_file, deltaTE);
tstart = tic; fprintf('%s:[FSL] Calculating a fieldmap using FSL topup:\n%s\n', datetime, command);
[status_fsl_prepare_fieldmap,result_fsl_prepare_fieldmap] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

status_fsl_prepare_fieldmap

result_fsl_prepare_fieldmap


%%
imain_file  = strcat(fsl_output_path, 'img');
datain_file = strcat(fsl_output_path, 'acq_param.txt');
%\\wsl.localhost\Ubuntu-20.04\home\image\fsl\pkgs\fsl-topup-2203.5-hb6de94e_0\etc\flirtsch
config_file = strcat(fsl_path, '/pkgs/fsl-topup-2203.5-hb6de94e_0/etc/flirtsch/b02b0.cnf');
out_file    = strcat(fsl_output_path, 'topup');
fout_file   = strcat(fsl_output_path, 'topup_fieldmap');
iout_file   = strcat(fsl_output_path, 'topup_img_cor');
command = sprintf('%s/topup --imain=%s --datain=%s --config=%s --out=%s --fout=%s --iout=%s --verbose', fsl_command, imain_file, datain_file, config_file, out_file, fout_file, iout_file);
tstart = tic; fprintf('%s:[FSL] Calculating a fieldmap using FSL topup:\n%s\n', datetime, command);
%[status_topup,result_topup] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));



%%
return

    %% Calculate a fieldmap using topup
    %----------------------------------------------------------------------
    % Part of FSL (ID: "")
    % topup
    %
    % Usage:
    % topup --imain=<some 4D image> --datain=<text file> --config=<text file with parameters> --out=my_topup_results
    %
    %
    % Compulsory arguments (You MUST set one or more of):
    %         --imain         name of 4D file with images
    %
    % Optional arguments (You may optionally specify one or more of):
    %         --datain        name of text file with PE directions/times
    %         --acqp          alternative way to specify text file with PE directions/times
    %         --out           base-name of output files (spline coefficients (Hz) and movement parameters)
    %         --fout          name of image file with field (Hz)
    %         --iout          name of 4D image file with unwarped images
    %         --featout       base-name of output for export to FEAT
    %         --logout        Name of log-file
    %         --warpres       (approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10
    %         --subsamp       sub-sampling scheme, default 1
    %         --fwhm          FWHM (in mm) of gaussian smoothing kernel, default 8
    %         --config        Name of config file specifying command line arguments
    %         --miter         Max # of non-linear iterations, default 5
    %         --lambda        Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documetation.
    %         --ssqlambda     If set (=1), lambda is weighted by current ssq, default 1
    %         --regmod        Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy
    %         --estmov        Estimate movements if set, default 1 (true)
    %         --minmet        Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)
    %         --splineorder   Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3
    %         --numprec       Precision for representing Hessian, double or float. Default double
    %         --interp        Image interpolation model, linear or spline. Default spline
    %         --scale         If set (=1), the images are individually scaled to a common mean, default 0 (false)
    %         --regrid        If set (=1), the calculations are done in a different grid, default 1 (true)
    %         --nthr          Number of threads to use (cannot be greater than numbers of hardware cores), default 1
    %         -h,--help       display help info
    %         -v,--verbose    Print diagonostic information while running
    %         -h,--help       display help info
    %----------------------------------------------------------------------
    imain_file  = strcat(fsl_output_path, 'img');
    datain_file = strcat(fsl_output_path, 'acq_param.txt');
    %\\wsl.localhost\Ubuntu-20.04\home\image\fsl\pkgs\fsl-topup-2203.5-hb6de94e_0\etc\flirtsch
    config_file = strcat(fsl_path, '/pkgs/fsl-topup-2203.5-hb6de94e_0/etc/flirtsch/b02b0.cnf');
    out_file    = strcat(fsl_output_path, 'topup');
    fout_file   = strcat(fsl_output_path, 'topup_fieldmap');
    iout_file   = strcat(fsl_output_path, 'topup_img_cor');
    command = sprintf('%s/topup --imain=%s --datain=%s --config=%s --out=%s --fout=%s --iout=%s --verbose', fsl_command, imain_file, datain_file, config_file, out_file, fout_file, iout_file);
    tstart = tic; fprintf('%s:[FSL] Calculating a fieldmap using FSL topup:\n%s\n', datetime, command);
    [status_topup,result_topup] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));



