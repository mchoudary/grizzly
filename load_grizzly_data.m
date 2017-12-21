%% Load data for template attacks from grizzly experiment
% Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none') % Remove TeX interpretation
tic

%% Load data
file_data = 'e2_bat_fb_beta_raw_s_0_3071.raw';
fprintf('Mapping data\n');
[m_data, metadata] = get_mmap(file_data);
nr_points = metadata.nr_points;
nr_groups = metadata.nr_groups;
nr_trials = metadata.nr_trials;
nr_blocks = floor(nr_trials / nr_groups);
toc

% The data is now available in m_data.
% See the file do_test_success_templates_bat_fb_dlinear.m for an example of
% how to use it to implement the template attacks and obtain some of
% the results that we shown in our paper.