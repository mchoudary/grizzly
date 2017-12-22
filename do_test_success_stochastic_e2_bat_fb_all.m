% Test stochastic model attacks 
% Author: Omar Choudary

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none'); % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
fmap = 'e2_bat_fb_beta_raw_s_0_3071.raw';
data_title = 'Stochastic A2 BAT FB';
path_data = 'results/';
name_data = sprintf('a2_bat_fb_stochastic_f9_all_dlinear_n200r_slr_g1000_r10.mat');
rand_iter = 10;
n_profile = 200; % ensure that: n_profile + n_attack < nr_blocks
nr_traces_vec = [1:10, 20:10:100, 200, 500, 1000];
rng('default'); % use same randomisation to get consistent results

%% Load files
fprintf('Mapping data\n');
[m_data, metadata] = get_mmap(fmap);
toc

%% Select idx for profile/attack
nr_groups = metadata.nr_groups;
nr_blocks = floor(metadata.nr_trials / nr_groups);
idx = 1:nr_blocks;
idx_profile = union(find(mod(idx,3) == 1), find(mod(idx,3) == 2));
idx_attack = find(mod(idx,3) == 0);
idx_profile_all = [];
for k=1:length(idx_profile)
    v = (idx_profile(k)-1)*nr_groups+1:idx_profile(k)*nr_groups;
    idx_profile_all = [idx_profile_all; v(:)];
end

%% Select basis
base = 'F9';

%% Select different factors for n_profile
% This is used to get results compatible with my PCA and LDA methods,
% where the number of traces used for computing the stochastic parameters
% is either the number of traces per byte (n_profile) or the number of
% traces per byte times the number of bytes
%nvec = {1, nr_groups, ceil(nr_groups^(1/2)), ceil(nr_groups^(1/3)), ceil(nr_groups^(1/4))};
nvec = {nr_groups};

%% Set up attack/result cells
results = cell(length(nvec), 6);

%% Run stochastic attacks for each combination of N and compression
for k=1:length(nvec)
    %% Select traces for stochastic model
    idx_profile_sel = idx_profile_all(randi([1, length(idx_profile_all)], 1, n_profile*nvec{k}));
    
    %% Run attack for LDA
    atype = 'templatelda';
    aparams = [];
    aparams.lda_threshold = 0.95;
    aparams.lda_dimensions = 4;
    discriminant = 'linearnocov';
    eparams = [];
    results{k, 1} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
    
    %% Run attack for PCA
    atype = 'templatepca';
    aparams = [];
    aparams.pca_threshold = 0.95;
    aparams.pca_alternate = 0;
    aparams.pca_dimensions = 4;
    discriminant = 'linear';
    eparams = [];
    results{k, 2} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
    
    %% Run attack for classic, 1ppc
    atype = 'same_profile';
    aparams = [];
    aparams.signal = 'dom';
    aparams.sel = '1ppc';
    aparams.p1 = 240;
    aparams.p2 = 0.95;
    discriminant = 'linear';
    eparams = [];
    results{k, 3} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
    
    %% Run attack for classic, 3ppc
    atype = 'same_profile';
    aparams = [];
    aparams.signal = 'dom';
    aparams.sel = '3ppc';
    aparams.p1 = 240;
    aparams.p2 = 0.95;
    discriminant = 'linear';
    eparams = [];
    results{k, 4} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
    
    %% Run attack for classic, 20ppc
    atype = 'same_profile';
    aparams = [];
    aparams.signal = 'dom';
    aparams.sel = '20ppc';
    aparams.p1 = 240;
    aparams.p2 = 0.95;
    discriminant = 'linear';
    eparams = [];
    results{k, 5} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
    
    %% Run attack for classic, allap
    atype = 'same_profile';
    aparams = [];
    aparams.signal = 'dom';
    aparams.sel = 'allap';
    aparams.p1 = 0.95;
    aparams.p2 = [];
    discriminant = 'linear';
    eparams = [];
    results{k, 6} = run_stochastic_attack(...
        m_data, metadata, idx_profile_sel, base, ...
        m_data, metadata, idx_attack, ...
        atype, aparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
end

%% Save all variables and clean up
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results', '-v7.3');
toc

%% Exit when running in script mode
% exit
