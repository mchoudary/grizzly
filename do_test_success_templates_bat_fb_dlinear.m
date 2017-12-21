%% Test templates
% Author: Omar Choudary

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none') % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
fmap = 'e2_bat_fb_beta_raw_s_0_3071.raw';
data_title = 'Templates A2 BAT FB';
path_data = 'results/';
name_data = sprintf('a2_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat');
rand_iter = 10;
n_profile = 200; % ensure that: n_profile + n_attack < nr_blocks
nr_traces_vec = [1:10, 20:10:100, 200, 500, 1000];
bytes = 0:255;
atype = 'mvn';

%% Load files
fprintf('Mapping data\n');
[m_data, metadata] = get_mmap(fmap);
toc

%% Select idx for profile/attack
nr_blocks = 3072;
idx = 1:nr_blocks;
idx_profile = union(find(mod(idx,3) == 1), find(mod(idx,3) == 2));
idx_profile = idx_profile(randi([1, length(idx_profile)], 1, n_profile));
idx_attack = find(mod(idx,3) == 0);

%% Set up attack/result cells
results = cell(6, 1);

%% Run attack for LDA
cmethod = 'LDA';
cparams.lda_dimensions = 4;
discriminant = 'linearnocov';
results{1} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Run attack for PCA
cmethod = 'PCA';
cparams = [];
cparams.pca_threshold = 0.95;
cparams.pca_alternate = 0;
cparams.pca_dimensions = 4;
discriminant = 'linear';
results{2} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Run attack for 1ppc
cmethod = 'sample';
cparams = [];
cparams.curve = 'dom';
cparams.sel = '1ppc';
cparams.p1 = 240;
discriminant = 'linear';
results{3} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Run attack for 3ppc
cmethod = 'sample';
cparams = [];
cparams.curve = 'dom';
cparams.sel = '3ppc';
cparams.p1 = 240;
discriminant = 'linear';
results{4} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Run attack for 20ppc
cmethod = 'sample';
cparams = [];
cparams.curve = 'dom';
cparams.sel = '20ppc';
cparams.p1 = 240;
discriminant = 'linear';
results{5} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Run attack for allap
cmethod = 'sample';
cparams = [];
cparams.curve = 'dom';
cparams.sel = 'allap';
cparams.p1 = 0.95;
discriminant = 'linear';
results{6} = run_template_attack(...
    m_data, metadata, idx_profile, m_data, metadata, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, rand_iter, nr_traces_vec);

%% Save all variables and clean up
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results');
toc

%% Exit, only use for condor. Matlab will quit.
% exit

