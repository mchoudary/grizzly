%% Test templates
% Author: Omar Choudary

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none') % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
fmap_profile = 'e2_bat_fb_alpha_raw_s_0_3071.raw';
fmap_attack = 'e2_bat_fb_beta_raw_s_0_3071.raw';
data_title = 'Templates A2D';
path_data = 'results/';
name_data = sprintf('a2d_ab_bat_fb_templates_adapt_boffset_dlinear_n1000r_slr_g1000_r10.mat');
rand_iter = 10;
n_profile = 1000; % ensure that: n_profile + n_attack < nr_blocks
nr_traces_vec = [1:10, 20:10:100, 200, 500, 1000];
bytes = 0:255;
atype = 'boffset';
rng('default'); % To use always same random values

%% Load files
fprintf('Mapping data for profile\n');
[m_data_profile, metadata_profile] = get_mmap(fmap_profile);
toc

fprintf('Mapping data for attack\n');
[m_data_attack, metadata_attack] = get_mmap(fmap_attack);
toc

%% Select idx for profile/attack
nr_blocks = 3072;
idx = 1:nr_blocks;
idx_profile = union(find(mod(idx,3) == 1), find(mod(idx,3) == 2));
idx_profile = idx_profile(randi([1, length(idx_profile)], 1, n_profile));
idx_attack = find(mod(idx,3) == 0);

%% Set up attack/result cells
results = cell(4, 1);

%% Run attack for LDA, m=4
cmethod = 'LDA';
cparams.lda_dimensions = 4;
cparams.lda_threshold = 0.95;
discriminant = 'linearnocov';
eparams = [];
s_profile = [];
s_profile.nr_sets = 1;
s_profile.mmap_data{1} = m_data_profile;
s_profile.metadata{1} = metadata_profile;
s_profile.idx{1} = idx_profile;
s_helper = [];
results{1} = run_template_attack_adapt(...
    s_profile, s_helper, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams);

%% Run attack for LDA, m=5
cmethod = 'LDA';
cparams.lda_dimensions = 5;
cparams.lda_threshold = 0.95;
discriminant = 'linearnocov';
eparams = [];
s_profile = [];
s_profile.nr_sets = 1;
s_profile.mmap_data{1} = m_data_profile;
s_profile.metadata{1} = metadata_profile;
s_profile.idx{1} = idx_profile;
s_helper = [];
results{2} = run_template_attack_adapt(...
    s_profile, s_helper, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams);

%% Run attack for PCA, m=4
cmethod = 'PCA';
cparams = [];
cparams.pca_threshold = 0.95;
cparams.pca_alternate = 0;
cparams.pca_dimensions = 4;
discriminant = 'linear';
eparams = [];
s_profile = [];
s_profile.nr_sets = 1;
s_profile.mmap_data{1} = m_data_profile;
s_profile.metadata{1} = metadata_profile;
s_profile.idx{1} = idx_profile;
s_helper = [];
results{3} = run_template_attack_adapt(...
    s_profile, s_helper, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams);

%% Run attack for PCA, m=5
cmethod = 'PCA';
cparams = [];
cparams.pca_threshold = 0.95;
cparams.pca_alternate = 0;
cparams.pca_dimensions = 5;
discriminant = 'linear';
eparams = [];
s_profile = [];
s_profile.nr_sets = 1;
s_profile.mmap_data{1} = m_data_profile;
s_profile.metadata{1} = metadata_profile;
s_profile.idx{1} = idx_profile;
s_helper = [];
results{4} = run_template_attack_adapt(...
    s_profile, s_helper, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams);
    
%% Save results
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results', '-v7.3');
toc

%% Exit, only use for condor. Matlab will quit.
% exit

