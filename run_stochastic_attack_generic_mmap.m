function [results] = run_stochastic_attack_generic_mmap(...
    m_data_profile, D_profile_all, idx_profile_all, ...
    m_data_attack, D_attack_all, idx_attack_group, ...
    V_profile, V_attack, V_discriminant, ...
    base, atype, aparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams)
%RUN_STOCHASTIC_ATTACK_GENERIC_MMAP Runs a generic stochastic attack
%   [results] = RUN_STOCHASTIC_ATTACK_GENERIC_MMAP(...
%       m_data_profile, D_profile_all, idx_profile_all, ...
%       m_data_attack, D_attack_all, idx_attack_group, ...
%       V_profile, V_attack, V_discriminant, ...
%       base, atype, aparams, discriminant, ...
%       rand_iter, nr_traces_vec, eparams)
%   runs a stochastic attack on generic data given via a memory mapped object
%   and returns a results structure that is defined below.
%
%   This method allows attacks on arbitrary data values and number of bits.
%
%   m_data_profile should be a memory mapped object containing the data for
%   profile in the field "m_data_profile.data(1).X". Here X should have
%   size nr_samples x nr_trials, where nr_trials is the number of traces
%   each having nr_samples.
%
%   D_profile_all should be a vector of length nr_trials containing the
%   data values corresponding to all the traces in m_data_profile.data(1).X.
%
%   idx_profile_all should be a vector of indices specifying which traces
%   should be used for the profile data. Note: this vector should provide
%   the complete list of all the traces to be used for profiling, not the
%   traces per group as in the case of idx_attack_group (see below).
%
%   m_data_attack should be a memory mapped object containing the data for
%   attack in the field "m_data_attack.data(1).X". Here X should have
%   size nr_samples x nr_trials, where nr_trials is the number of traces
%   each having nr_samples.
%
%   D_attack_all should be a vector of length nr_trials containing the data
%   values corresponding to all the traces in m_data_attack.data(1).X.
%
%   idx_attack_group should be a vector of indices specifying which traces
%   should be used for the attack data. This vector should specify only the
%   traces to be used for each attack value/group specified in V_attack.
%
%   V_profile should be a vector with indices specifying for which values
%   (out of those in D_profile) to estimate profiling parameters. For each
%   attack type (see atype below), this vector is used as follows:
%   -> 'selection' and 'classic: this vector defines the values that
%   are used for the computation of signal strength estimates.
%   -> 'pca' and 'lda': this vector defines the subset of values used to
%   estimate the PCA or LDA parameters.
%   -> 'templatepca' and 'templatelda': this vector defines the values that
%   are used for the estimation of PCA and LDA parameters.
%
%   V_attack should be a vector with indices specifying which values (out
%   of those in D_attack_all) will be used to obtain the attack data. This
%   data will then be used for the estimation of the attack result,
%   e.g. the guessing entropy.
%
%   V_discriminant should be a vector with indices specifying for which
%   values (of similar class as those in D_attack_all) to compute
%   the discriminant of the attack, i.e. which values are used for the
%   estimation of the guessing entropy or other success information.
%   For example, it is possible to  provide data from only some values
%   through X_profile but to estimate the mean vectors for many more values.
%
%   The difference between V_attack and V_discriminant is that V_attack
%   specifies the real data that is used to estimate the attack result,
%   while V_discriminant specifies the data that is only estimated (e.g.
%   the means via stochastic models) in order to compare traces from the
%   attack data with estimated data on all possible values V_discriminant.
%   This is useful for the "Partial Guessing Entropy".
%
%   base should be a string specifying the base vectors to be used in the
%   computation of stochastic parameters and attack. See 'get_map_base' for
%   the possible options (including 'F9', 'F17').
%
%   atype should be a string specifying the type of stochastic attack to be
%   used. Currently supported are:
%   - 'classic': which uses the classic approach described in the CHES paper.
%   - 'selection': uses the same profiling traces to compute both the
%   base coeficients and the covariance matrix and a signal strength
%   estimate to determine which points to select from a trace for
%   compression.
%   - 'pca': first compute the sum of squares and corss product matrices
%   and the PCA parameters from a subset of traces, then
%   projects traces and then computes stochastic coefficients.
%   - 'templatepca': method very similar to the PCA used with template
%   attacks, only adapted to compute the "mean vectors" using the
%   stochastic model.
%   - 'lda': similar to 'pca' but uses Fisher's LDA instead.
%   - 'templatelda': method very similar to the LDA used with template
%   attacks, only adapted to compute the "mean vectors" using the
%   stochastic model.
%
%   aparams should be a structure of params specific to the given attack
%   type. For each attack type the parameters are as follows:
%   - 'classic':
%       -> aparams.signal: is a string specifying the algorithm used to
%       to compute the signal strength estimate. Possible choices are:
%       'dom', 'snr', 'sosd', 'bnorm' or 'bnorm_std'.
%       -> aparams.sel is a string specifying the class of selection.
%       -> aparams.p1 is a parameter for the class of selection.
%       -> aparams.p2 is an additional optional parameter.
%   - 'selection': same as for 'classic'.
%   - 'pca':
%       -> aparams.m_data_subset: memory mapped object, similar to
%       m_data_profile, but containing data for subsets from which PCA
%       params will be obtained.
%       -> aparams.D_subset_all: similar to D_profile_all but for m_data_subset.
%       -> aparams.idx_traces: vector with index of traces per group to be
%       used with the data in m_data_subset.
%       -> aparams.pca_threshold: used for PCA
%       -> aparams.pca_alternate: used for PCA
%       -> aparams.pca_dimensions: used for PCA
%       -> aparams.cov_from_sel: set this to 1 or true to compute the
%       covariance from the selection used to compute the PCA parameters.
%       Otherwise, if this is 0 or false, the covariance will be computed
%       from the traces used to compute the stochastic model coefficients.
%   - 'templatepca:
%       -> aparams.pca_threshold: used for PCA
%       -> aparams.pca_alternate: used for PCA
%       -> aparams.pca_dimensions: used for PCA
%   - 'lda':
%       -> aparams.m_data_subset: memory mapped object, similar to
%       m_data_profile, but containing data for subsets from which PCA
%       params will be obtained.
%       -> aparams.D_subset_all: similar to D_profile_all but for m_data_subset.
%       -> aparams.idx_traces: vector with index of traces per group to be
%       used with the data in m_data_subset.
%       -> aparams.lda_threshold: used for LDA
%       -> aparams.lda_dimensions: used for LDA
%       -> aparams.cov_from_sel: set this to 1 or true to compute the
%       covariance from the selection used to compute the LDA parameters.
%       Otherwise, if this is 0 or false, the covariance will be computed
%       from the traces used to compute the stochastic model coefficients.
%   - 'templatelda':
%       -> aparams.lda_threshold: used for LDA
%       -> aparams.lda_dimensions: used for LDA
%
%   discriminant should be a string specifying the type of discriminant to
%   be used. The possible options are:
%   - 'linear': uses a pooled common covariance matrix with a linear
%      discriminant.
%   - 'linearnocov': does not use a covariance matrix. Might be useful in
%      particular with LDA, where the covariance should be the
%      identity if the eigenvectors are chosen carefully.
%   - 'linearfast': uses the linear discriminant with some precomputed
%      values to obtain a faster evaluation of success information.
%   Note that there is no 'linearfastnocov' since it is faster to use
%   'linearnocov' than using 'linearfast' with an identity covariance (i.e
%   nocov).
%
%   rand_iter should be a positive integer specifying the number of
%   iterations to run the evaluation (guessing_entropy) computation. The
%   returned results may contain either the individual or the average
%   results. Check below for details.
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used for each element.
%
%   eparams is a structure of extra parameters that may be needed for some
%   options. Some of these extra parameters include:
%   - 'save_xdata': if this field exists and is non zero then the x_profile
%   and x_attack data will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'save_eval': if this field exists and is non zero then the
%   handle_eval pointer and the data for evaluation (generally the
%   templates) will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'save_signals': set this to 1 to save signal strength estimates
%   (for classic and selection only).
%   - 'save_ssp': use this to save the SSP matrices M,B,W (only for PCA, LDA).
%   - 'v_attack': if this field is given then it should contain a vector of
%   the same length this method's parameter V_attack. This vector will
%   replace V_attack. This may be useful to provide V_attack as values
%   covering a wide range (e.g. 16-bit data) but then running the
%   evaluation on only a part of the data (e.g. the most significant 8
%   bits). I used it to test attacks on 8-bit data influenced by pipeline
%   of other 8-bit data. Hence the original V_attack had 16-bit data but
%   then I run the attack on only one byte, assuming the remaining byte is
%   noise that the attack has to deal with.
%
%   The 'results' structure contains the following:
%   -> results.atype: the atype string.
%   -> results.aparams: the aparams structure.
%   -> results.discriminant: the discriminant string.
%   -> results.rand_iter: the number of iterations.
%   -> results.nr_traces_vec: the vector with number of attack traces.
%   -> results.coef: the matrix of coefficients for the stochastic model
%   -> results.signals: vectors of signal strength estimates (optional).
%   -> results.M: the matrix of group means (optional).
%   -> results.B: the between-groups matrix (optional).
%   -> results.W: the matrix of variances and covariances across all data (optional).
%   -> results.x_profile: the profiling data, after compression. (optional)
%   -> results.x_attack: the attack data, after compression. (optional)
%   -> results.handle_prepare: function used to extract features for templates.
%   -> results.pp1 ... results.pp5: parameters of results.handle_prepare.
%   -> results.handle_eval: function used to evaluate templates. (optional)
%   -> results.pe3 ... results.pe6: parameters for results.handle_eval. (optional)
%   -> results.success_info: guessing entropy information, as returned by
%   the get_success_info_like method.
%   -> results.error: an optional messsage in case of an error. This may
%   avoid the abortion of a larger set of runs while allowing to detect the
%   cases in which an error occurred.
%
%   See the papers "A Stochastic Model for Differential Side Channel
%   Cryptanalysis", Schindler et al. 2005, CHES 2005,
%   and "Efficient Template Attacks", Choudary and Kuhn, CARDIS 2013.
%
%   See also get_mmap.

%% Check and initialise parameters
np = length(idx_profile_all);
results = [];
results.atype = atype;
results.aparams = aparams;
results.discriminant = discriminant;
results.rand_iter = rand_iter;
results.nr_traces_vec = nr_traces_vec;
if nargin < 11
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_stochastic_attack() ...\n');

%% Start selected attack type
if strcmp(atype, 'classic')
    %% Select data for computing the coefficients
    fprintf('Obtaining data for stochastic coefficients...\n');
    N1 = floor(np/2);
    idx_n1 = idx_profile_all(1:N1);
    X1 = double(m_data_profile.data(1).X(:,idx_n1)');
    D1 = D_profile_all(idx_n1);
    N2 = np-N1;
    idx_n2 = idx_profile(N1+1:N1+N2);
    X2 = double(m_data_profile.data(1).X(:,idx_n2)');
    D2 = D_profile_all(idx_n2);
    toc
    
    %% Select basis for coefficients
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients on raw data...\n');
    [coef] = compute_coef_stochastic(X1, D1, map_base);
    toc
    
    %% Compute covariance
    fprintf('Computing raw covariance matrix...\n');
    Z2 = X2 - get_leakage_from_map_base(D2, coef, map_base);
    C = (Z2'*Z2) / (N2-1);
    toc
    
    %% Compute signal strengths
    fprintf('Computing signal strength estimate...\n');
    if strcmp(aparams.signal, 'bnorm') || strcmp(aparams.signal, 'bnorm_std')
        signals = get_signal_strength_coef(X1, coef, base);
    else
        smean_r = get_leakage_from_map_base(V_profile, coef, map_base);
        signals = get_signal_strength_ssp(smean_r, [], C*(N2-1), N2);
    end
    toc
    if isfield(eparams, 'save_signals') && (eparams.save_signals ~= 0)
        results.signals = signals;
    end

    %% Select points (samples)
    spoints = get_selection(signals.(aparams.signal), aparams.sel, ...
                            aparams.p1, aparams.p2);
                        
    %% Restrict coefficiencts only to selected points
    coef = coef(:, spoints);
    results.coef = coef;
    
    %% Obtain mean and covariance on selected points
    
    %% Compute compressed mean and covariance for attack step
    fprintf('Obtaining smean/scov for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    scov = C(spoints, spoints);
    toc
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template;
    pp1 = spoints;
    pp2 = [];
    pp3 = [];
    pp4 = [];
    pp5 = [];    
elseif strcmp(atype, 'selection')
    %% Obtain data for computing the coefficients
    fprintf('Obtaining data for stochastic coefficients...\n');
    X = double(m_data_profile.data(1).X(:,idx_profile_all)');
    D = D_profile_all(idx_profile_all);
    toc
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients on raw data...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Compute covariance
    fprintf('Computing raw covariance matrix...\n');
    Z = X - get_leakage_from_map_base(D, coef, map_base);
    C = (Z'*Z) / (np-1);
    toc
    
    %% Compute signal strengths
    fprintf('Computing signal strength estimate...\n');
    if strcmp(aparams.signal, 'bnorm') || strcmp(aparams.signal, 'bnorm_std')
        signals = get_signal_strength_coef(X, coef, base);
    else
        smean_r = get_leakage_from_map_base(V_profile, coef, map_base);
        signals = get_signal_strength_ssp(smean_r, [], C*(np-1), np);
    end    
    toc
    if isfield(eparams, 'save_signals') && (eparams.save_signals ~= 0)
        results.signals = signals;
    end

    %% Select points (samples)
    spoints = get_selection(signals.(aparams.signal), aparams.sel, ...
                            aparams.p1, aparams.p2);
                        
    %% Restrict coefficiencts only to selected points
    coef = coef(:, spoints);
    results.coef = coef;
    
    %% Compute compressed mean and covariance for attack step
    fprintf('Obtaining smean/scov for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    scov = C(spoints, spoints);
    toc
                        
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template;
    pp1 = spoints;
    pp2 = [];
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(atype, 'pca')    
    %% Obtain the SSP matrices (M, B, W) on selected bytes
    fprintf('Obtaining sums of squares and cross products on selected bytes...\n');
    [M, B, W] = compute_ssp_generic_mmap(aparams.m_data_subset, ...
                                         aparams.D_subset_all, ...
                                         V_profile, aparams.idx_traces);
    if isfield(eparams, 'save_ssp') && (eparams.save_ssp ~= 0)
        results.M = M;
        results.B = B;
        results.W = W;
    end
    xmm = mean(M, 1);
    toc
                                
    %% Compute PCA params
    fprintf('Computing PCA parameters...\n');
    [U, ~, ~, K] = compute_params_pca(M, aparams.pca_threshold);
    if isfield(aparams, 'pca_dimensions') && (aparams.pca_dimensions > 0)
        U = U(:,1:aparams.pca_dimensions);
    else
        U = U(:,1:K);
    end
    toc
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = U;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
    %% Obtain data for computing the coefficients   
    fprintf('Obtaining data for stochastic coefficients...\n');
    L = double(m_data_profile.data(1).X(:,idx_profile_all)');
    X = handle_prepare(L, pp1, pp2, pp3, pp4, pp5);
    D = D_profile_all(idx_profile_all);
    toc
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors from stochastic model
    fprintf('Aproximating mean vectors for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    toc
    
    %% Compute covariance
    fprintf('Computing covariance...\n');
    if aparams.cov_from_sel
        fprintf('Computing data for covariance from selection...\n');
        x_cov = compute_features_generic_mmap(aparams.m_data_subset, ...
                                    aparams.D_subset_all, ...
                                    V_profile, aparams.idx_traces, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5);
        toc
        [~, C] = compute_template(x_cov);
        scov = mean(C, 3);
    else
        Z = X - get_leakage_from_map_base(D, coef, map_base);
        scov = (Z'*Z) / (np-1);
    end
    toc
elseif strcmp(atype, 'templatepca')
    %% Obtain data for computing the coefficients
    fprintf('Obtaining data for stochastic coefficients...\n');
    X = double(m_data_profile.data(1).X(:,idx_profile_all)');
    D = D_profile_all(idx_profile_all);
    toc
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients for raw data...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Approximate mean vectors from stochastic model
    fprintf('Aproximating raw mean vectors from stochastic model...\n');
    smean_r = get_leakage_from_map_base(V_profile, coef, map_base);
    toc
    
    %% Compute PCA params
    fprintf('Computing PCA parameters...\n');
    [U, ~, xmm, K] = compute_params_pca(smean_r, aparams.pca_threshold, ...
                                        aparams.pca_alternate);
    if isfield(aparams, 'pca_dimensions') && (aparams.pca_dimensions > 0)
        U = U(:,1:aparams.pca_dimensions);
    else
        U = U(:,1:K);
    end
    toc
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = U;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
    %% Project data using PCA
    Y = handle_prepare(X, pp1, pp2, pp3, pp4, pp5);
    
    %% Compute Stochastic coefficients in PCA space
    fprintf('Computing Stochastic coefficients for compressed data...\n');
    [coef] = compute_coef_stochastic(Y, D, map_base);
    results.coef = coef;
    toc
    
    %% Approximate mean vectors in PCA space
    fprintf('Aproximating mean vectors for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    toc
    
    %% Compute covariance in PCA space
    fprintf('Computing covariance...\n');
    Z = Y - get_leakage_from_map_base(D, coef, map_base);
    scov = (Z'*Z) / (np-1);
    toc
elseif strcmp(atype, 'lda')
    %% Obtain the SSP matrices (M, B, W) on selected bytes
    fprintf('Obtaining sums of squares and cross products on selected bytes...\n');
    [M, B, W] = compute_ssp_generic_mmap(aparams.m_data_subset, ...
                                         aparams.D_subset_all, ...
                                         V_profile, aparams.idx_traces);
    if isfield(eparams, 'save_ssp') && (eparams.save_ssp ~= 0)
        results.M = M;
        results.B = B;
        results.W = W;
    end
    xmm = mean(M, 1);
    toc
    
    %% Compute LDA params
    fprintf('Computing Fishers LDA parameters...\n');
    nr_values_profile = length(V_profile);
    nr_traces_per_value = length(aparams.idx_traces);
    Spool = W / (nr_values_profile*(nr_traces_per_value-1));
    [A ,~, K] = compute_params_lda(B, Spool, nr_values_profile, aparams.lda_threshold);
    if isfield(aparams, 'lda_dimensions') && (aparams.lda_dimensions > 0)
        FW = A(:,1:aparams.lda_dimensions);        
    else
        FW = A(:,1:K);
    end
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = FW;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
    %% Obtain data for computing the coefficients    
    fprintf('Obtaining data for stochastic coefficients...\n');
    L = double(m_data_profile.data(1).X(:,idx_profile_all)');
    X = handle_prepare(L, pp1, pp2, pp3, pp4, pp5);
    D = D_profile_all(idx_profile_all);
    toc
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors from stochastic model
    fprintf('Aproximating mean vectors for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    toc
    
    %% Compute covariance
    fprintf('Computing covariance...\n');
    if aparams.cov_from_sel
        fprintf('Computing data for covariance from selection...\n');
        x_cov = compute_features_generic_mmap(aparams.m_data_subset, ...
                                    aparams.D_subset_all, ...
                                    V_profile, aparams.idx_traces, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5);
        toc
        [~, C] = compute_template(x_cov);
        scov = mean(C, 3);
    else
        Z = X - get_leakage_from_map_base(D, coef, map_base);
        scov = (Z'*Z) / (np-1);
    end
    toc
elseif strcmp(atype, 'templatelda')
    %% Obtain data for computing the coefficients
    fprintf('Obtaining data for stochastic coefficients...\n');
    X = double(m_data_profile.data(1).X(:,idx_profile_all)');
    D = D_profile_all(idx_profile_all);
    toc
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients for raw data...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Approximate raw mean vectors from stochastic model
    fprintf('Aproximating raw mean vectors from stochastic model...\n');
    smean_r = get_leakage_from_map_base(V_profile, coef, map_base);
    toc
    
    %% Compute raw covariance matrix
    fprintf('Computing raw covariance...\n');
    Z = X - get_leakage_from_map_base(D, coef, map_base);
    C = (Z'*Z) / (np-1);
    toc
    
    %% Compute raw between-groups matrix B
    fprintf('Computing between-groups matrix B...\n');
    nr_values_profile = length(V_profile);
    xmm = mean(smean_r, 1);
    T = smean_r - ones(nr_values_profile, 1)*xmm;
    B = T'*T;
    toc
    
    %% Compute LDA params
    fprintf('Computing Fishers LDA parameters...\n');
    [A ,~, K] = compute_params_lda(B, C, nr_values_profile, aparams.lda_threshold);
    if isfield(aparams, 'lda_dimensions') && (aparams.lda_dimensions > 0)
        FW = A(:,1:aparams.lda_dimensions);
    else
        FW = A(:,1:K);
    end
    toc
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = FW;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
    %% Project data using LDA
    Y = handle_prepare(X, pp1, pp2, pp3, pp4, pp5);
    
    %% Compute Stochastic coefficients in LDA space
    fprintf('Computing Stochastic coefficients for compressed data...\n');
    [coef] = compute_coef_stochastic(Y, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors in LDA space
    fprintf('Aproximating mean vectors for attack...\n');
    smean = get_leakage_from_map_base(V_discriminant, coef, map_base);
    toc
    
    %% Compute covariance in LDA space
    fprintf('Computing covariance...\n');
    Z = Y - get_leakage_from_map_base(D, coef, map_base);
    scov = (Z'*Z) / (np-1);
    toc
else
    error('Unknown atype: %s', atype);
end

%% Store handle_prepare data
results.handle_prepare = handle_prepare;
results.pp1 = pp1;
results.pp2 = pp2;
results.pp3 = pp3;
results.pp4 = pp4;
results.pp5 = pp5;

%% Load data for attack
fprintf('Computing attack data...\n');
X_attack = compute_features_generic_mmap(m_data_attack, D_attack_all, ...
                                    V_attack, idx_attack_group, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_attack = X_attack;
end
toc

%% Set evaluation parameters
fprintf('Computing evaluation parameters...\n');
tmiu = smean;
ic0 = inv(scov);
if strcmp(discriminant, 'linear')
    handle_discriminant = @compute_discriminant;
    pe3 = tmiu;
    pe4 = ic0;
    pe5 = [];
    pe6 = [];
elseif strcmp(discriminant, 'linearnocov')
    handle_discriminant = @compute_discriminant;
    pe3 = tmiu;
    pe4 = [];
    pe5 = [];
    pe6 = [];
elseif strcmp(discriminant, 'linearfast')
    [ng_miu, m] = size(tmiu);
    Y = zeros(ng_miu, m);
    Z = zeros(m, 1);
    for k=1:ng_miu
        Y(k,:) = tmiu(k,:)*ic0;
        Z(k) = Y(k,:)*tmiu(k,:)';
    end
    handle_discriminant = @compute_dlinear_fast;
    pe3 = Y;
    pe4 = Z;
    pe5 = [];
    pe6 = [];
else
    error('Unsupported discriminant type: %s', discriminant);
end
toc

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.handle_discriminant = handle_discriminant;
    results.pe3 = pe3;
    results.pe4 = pe4;
    results.pe5 = pe5;
    results.pe6 = pe6;
end

%% Replace V_attack if requested
% This might be needed to run attacks on only part of the target data while
% allowing to select attack data from the full range.
if isfield(eparams, 'v_attack')
    if length(eparams.v_attack) ~= length(V_attack)
        error('eparams.v_attack given but has incompatible length with V_attack');
    end
    V_attack = eparams.v_attack;
end

%% Compute the success information
fprintf('Computing success info...\n');
[results.success_info] = get_success_info_generic(...
                               X_attack, V_attack, V_discriminant,...
                               rand_iter, ...
                               nr_traces_vec, ...
                               handle_discriminant, pe3, pe4, pe5, pe6);
toc

end
