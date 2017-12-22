% Test stochastic model attacks 
% Author: Omar Choudary

%% Reset environment
% close all;
clear;
set(0, 'DefaulttextInterpreter', 'none'); % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
fmap = 'e5_all_bat_fb_2mhz_n500_beta.raw';
data_title = 'Stochastic A6 BAT FB';
path_data = 'results/';
name_data = sprintf('a6_bat_fb_stochastic_f17_all_dlinear_n1000r_p1_slr_g1000_r10.mat');
rand_iter = 10;
total_values = 2^16;
pmul = 2^4;
n_profile = 1000*pmul; % ensure that: n_profile + n_attack*nr_attack_groups < nr_trials
n_attack = 150;
n_profile_per_group = 100;
nr_traces_vec = [1, 2, 5, 10, 20, 50, 100];
rng('default'); % use same randomisation to get consistent results

%% Load files
fprintf('Mapping data\n');
[m_data, metadata] = get_mmap(fmap);
toc

%% Extract values in mapped data
fprintf('Extracting values corresponding to mapped data\n');
D_all_bytes = m_data.data(1).B(2:3,:)';
D_all = bitxor(bitshift(uint16(D_all_bytes(:,1)), 8), uint16(D_all_bytes(:,2)));
toc

%% Select traces for profile/attack
fprintf('Selecting traces for profile and attack...\n');
nr_trials = metadata.nr_trials;
total_trials_per_value = nr_trials / total_values;
idx_group = 1:total_trials_per_value;
idx_profile_group = find(mod(idx_group,4) == 0);
total_profile_trials_per_group = length(idx_profile_group);
idx_attack_group = union(union(find(mod(idx_group,4) == 1), find(mod(idx_group,4) == 2)), ...
                   find(mod(idx_group,4) == 3));
total_profile = total_values*total_profile_trials_per_group;
idx_profile_all = zeros(total_profile, 1);
[~, si] = sort(D_all, 1, 'ascend');
for k=1:total_values
    idx_val = si(total_trials_per_value*(k-1)+1:total_trials_per_value*k);
    idx_profile_all(total_profile_trials_per_group*(k-1)+1:total_profile_trials_per_group*k) = idx_val(idx_profile_group);
end
idx_profile_sel = idx_profile_all(randi([1, total_profile], n_profile, 1));
toc

%% Select base for coefficients
base = 'F17';

%% Select values for discriminant evaluation
V_discriminant = 0:(total_values-1);

%% Select number of values for V_profile and V_attack
ng_profile = {2^16}; % groups for profile
ng_attack = {2^8}; % groups for attack - note for 2^16 takes time...
nr_eval = 3; % How many times to evaluate with each selection

%% Set up attack/result cells
nr_attacks = 5;
results = cell(length(ng_attack), length(ng_profile), nr_eval, nr_attacks);

for r=1:nr_eval
    for k=1:length(ng_attack)    
        % Select groups for V_attack
        if ng_attack{k} == total_values
            V_attack = V_discriminant; % full guessing entropy
        else
            V_attack = randi([0, total_values-1], ng_attack{k}, 1); % partial guessing entropy
        end
    
        for i=1:length(ng_profile)
            fprintf('Running attacks for (k,i,r)=(%d,%d,%d)\n', k, i, r);
            % Select groups for V_profile
            if ng_profile{i} == total_values
                V_profile = 0:total_values-1;
            else
                V_profile = randi([0, total_values-1], ng_profile{i}, 1);
            end
            
           %% Run attack for selection, 20ppc
            atype = 'selection';
            aparams = [];
            aparams.signal = 'bnorm';
            aparams.sel = '20ppc';
            aparams.p1 = 40;
            aparams.p2 = 0.90;
            discriminant = 'linearfast';
            eparams = [];
            eparams.save_signals = 1;
            results{k,i,r,1} = run_stochastic_attack_generic_mmap(...
                m_data, D_all, idx_profile_sel, ...
                m_data, D_all, idx_attack_group, ...
                V_profile, V_attack, V_discriminant, ...
                base, atype, aparams, discriminant, ...
                rand_iter, nr_traces_vec, eparams);
            
           %% Run attack for S-PCA (aka Master PCA)
            atype = 'templatepca';
            aparams = [];            
            aparams.pca_threshold = 0.95;
            aparams.pca_alternate = 0;
            aparams.pca_dimensions = 10; % see semilogy of eigenvalues - clear stuff
            discriminant = 'linearfast';
            eparams = [];
            results{k,i,r,2} = run_stochastic_attack_generic_mmap(...
                m_data, D_all, idx_profile_sel, ...
                m_data, D_all, idx_attack_group, ...
                V_profile, V_attack, V_discriminant, ...
                base, atype, aparams, discriminant, ...
                rand_iter, nr_traces_vec, eparams);
            
           %% Run attack for T-PCA (aka Slave PCA)
            atype = 'pca';
            aparams = [];
            aparams.idx_traces = idx_profile_group;
            aparams.pca_threshold = 0.95;
            aparams.pca_alternate = 0;
            aparams.pca_dimensions = 10; % see semilogy of eigenvalues - clear stuff
            aparams.cov_from_sel = 0;
            discriminant = 'linearfast';
            eparams = [];
            results{k,i,r,3} = run_stochastic_attack_generic_mmap(...
                m_data, D_all, idx_profile_sel, ...
                m_data, D_all, idx_attack_group, ...
                V_profile, V_attack, V_discriminant, ...
                base, atype, aparams, discriminant, ...
                rand_iter, nr_traces_vec, eparams);
            
           %% Run attack for S-LDA (aka Master LDA)
            atype = 'templatelda';
            aparams = [];            
            aparams.lda_threshold = 0.95;
            aparams.lda_dimensions = 10; % see semilogy of eigenvalues
            discriminant = 'linearnocov';
            eparams = [];
            results{k,i,r,4} = run_stochastic_attack_generic_mmap(...
                m_data, D_all, idx_profile_sel, ...
                m_data, D_all, idx_attack_group, ...
                V_profile, V_attack, V_discriminant, ...
                base, atype, aparams, discriminant, ...
                rand_iter, nr_traces_vec, eparams);
            
           %% Run attack for T-LDA (aka Slave LDA)
            atype = 'lda';
            aparams = [];           
            aparams.idx_traces = idx_profile_group;
            aparams.lda_threshold = 0.95;
            aparams.lda_dimensions = 10; % see semilogy of eigenvalues
            aparams.cov_from_sel = 0;
            discriminant = 'linearfastnocov';
            eparams = [];
            results{k,i,r,5} = run_stochastic_attack_generic_mmap(...
                m_data, D_all, idx_profile_sel, ...
                m_data, D_all, idx_attack_group, ...
                V_profile, V_attack, V_discriminant, ...
                base, atype, aparams, discriminant, ...
                rand_iter, nr_traces_vec, eparams);            
        end
    end
end

%% Save all variables and clean up
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results', '-v7.3');
toc

%% Exit when running in script mode
% exit
