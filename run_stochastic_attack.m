function [results] = run_stochastic_attack(...
    m_data_profile, metadata_profile, ...
    idx_profile, base, ...
    m_data_attack, metadata_attack, idx_attack, ...
    atype, aparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams)
%RUN_STOCHASTIC_ATTACK Runs a stochastic model attack
%   [results] = RUN_STOCHASTIC_ATTACK(...
%       m_data_profile, metadata_profile, ...
%       idx_profile, base, ...
%       m_data_attack, metadata_attack, idx_attack, ...
%       atype, aparams, discriminant, ...
%       rand_iter, nr_traces_vec, eparams)
%   runs a stochastic attack with the given parameters and returns a results
%   structure that is defined below.
%
%   m_data_profile and metadata_profile should be the memory mapped object
%   and associated metadata info for the profiling data. Use get_mmap 
%   on the selected data to obtain these objects.
%
%   idx_profile should be a vector of indices specifying which traces
%   should be used for the profile data. Note: this vector should provide
%   the complete list of all the traces to be used for profiling, not the
%   traces per group as in the case of idx_attack (see below).
%
%   base should be a string specifying the base vectors to be used in the
%   computation of stochastic parameters and attack. The possible options
%   are:
%   - 'F9': use a constant factor plus the individual bits (1 to 8).
%
%   m_data_attack and metadata_attack should be the memory mapped object
%   and associated metadata info for the attack data. Use get_mmap 
%   on the selected data to obtain these objects. You can pass the same
%   objects for profiling and attack, which should make the attack faster.
%
%   idx_attack should be a vector of indices specifying which traces
%   should be used for the attack data. Unlike for idx_profile, this vector
%   should specify the traces to be used per group for the attack data.
%
%   atype should be a string specifying the type of stochastic attack to be
%   used. Currently supported are:
%   - 'classic': which uses the classic approach described in the CHES paper.
%   - 'same_profile': uses the same profiling traces to compute both the
%   base coeficients and the covariance matrix.
%   - 'pca': first compute the sum of squares and corss product matrices
%   and the PCA parameters from a subset of traces, then
%   projects traces and then computes stochastic coefficients.
%   - 'templatepca': method very similar to the PCA used with template
%   attacks, only adapted to compute the "mean vectors" using the
%   stochastic model.
%   - 'badpca': method proposed in "A New Difference Method for
%   Side-Channel Analysis with High-Dimensional Leakage Models". As the
%   name suggests, this is a very bad way of implementing PCA and I don't
%   recommended. Used here only to compare results.
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
%       'dom', 'snr', 'sosd', 'bnorm'.
%       -> aparams.sel is a string specifying the class of selection.
%       -> aparams.p1 is a parameter for the class of selection.
%       -> aparams.p2 is an additional optional parameter.
%   - 'same_profile': same as for 'classic'.
%   - 'pca':
%       -> aparams.byte_sel: vector having a selection of bytes to be used
%       for the computation of the SSP matrices.
%       -> aparams.idx_traces: vector with index of traces per group to be
%       used when computing the SSP matrices.
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
%   - 'badpca': same params as for 'templatepca'.
%   - 'lda':
%       -> aparams.byte_sel: vector having a selection of bytes to be used
%       for the computation of the SSP matrices.
%       -> aparams.idx_traces: vector with index of traces per group to be
%       used when computing the SSP matrices.
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
%   - 'save_ssp': use this to save the SSP matrices M,B,W (only for PCA, LDA).
%
%   The 'results' structure contains the following:
%   -> results.metadata_profile: the metadata structure for profile.
%   -> results.idx_profile: the idx_profile vector.
%   -> results.metadata_attack: the metadata structure for attack.
%   -> results.idx_attack: the idx_attack vector.
%   -> results.atype: the atype string.
%   -> results.aparams: the aparams structure.
%   -> results.discriminant: the discriminant string.
%   -> results.rand_iter: the number of iterations.
%   -> results.nr_traces_vec: the vector with number of attack traces.
%   -> results.coef: the matrix of coefficients for the stochastic model
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

%% Check and initialise parameters
nr_groups = metadata_profile.nr_groups;
np = length(idx_profile);
results = [];
results.metadata_profile = metadata_profile;
results.idx_profile = idx_profile;
results.metadata_attack = metadata_attack;
results.idx_attack = idx_attack;
results.atype = atype;
results.aparams = aparams;
results.discriminant = discriminant;
results.rand_iter = rand_iter;
results.nr_traces_vec = nr_traces_vec;
if nargin < 13
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_stochastic_attack() ...\n');

%% Start selected attack type
if strcmp(atype, 'classic')
    %% Obtain data for computing the coefficients
    N1 = floor(np/2);
    idx_n1 = idx_profile(1:N1);
    X1 = double(m_data_profile.data(1).X(:,idx_n1)');
    D1 = double(m_data_profile.data(1).B(2,idx_n1)');
    N2 = np-N1;
    idx_n2 = idx_profile(N1+1:end);
    X2 = double(m_data_profile.data(1).X(:,idx_n2)');
    D2 = double(m_data_profile.data(1).B(2,idx_n2)');
    
    %% Select basis for coefficients
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X1, D1, map_base);
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean_r = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance
    % (note I've observed classic and same_profiling give same/similar
    % covariance as I expected. Therefore I really don't see the point in
    % the stochastic approach paper about using different sets for this.
    Z2 = X2 - get_leakage_from_map_base(D2, coef, map_base);
    C = (Z2'*Z2) / (N2-1);
    
    %% Compute signal strengths
    if strcmp(aparams.signal, 'bnorm')
        signals = get_signal_strength_coef(X1, coef, base);
    else
        signals = get_signal_strength_ssp(smean_r, [], C*(nr_groups*(N2-1)), N2);
    end

    %% Select points (samples)
    spoints = get_selection(signals.(aparams.signal), aparams.sel, ...
                            aparams.p1, aparams.p2);
                        
    %% Restrict coefficiencts only to selected points
    coef = coef(:, spoints);
    results.coef = coef;
    
    %% Obtain mean and covariance on selected points
    smean = smean_r(:,spoints);
    scov = C(spoints, spoints);   
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template;
    pp1 = spoints;
    pp2 = [];
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
elseif strcmp(atype, 'same_profile')
    %% Obtain data for computing the coefficients    
    X = double(m_data_profile.data(1).X(:,idx_profile)');
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean_r = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance
    % (note I've observed classic and same_profiling give same/similar
    % covariance as I expected. Therefore I really don't see the point in
    % the stochastic approach paper about using different sets for this.
    Z = X - get_leakage_from_map_base(D, coef, map_base);
    C = (Z'*Z) / (np-1);
    
    %% Compute signal strengths
    if strcmp(aparams.signal, 'bnorm') || strcmp(aparams.signal, 'bnorm_std')
        signals = get_signal_strength_coef(X, coef, base);
    else
        signals = get_signal_strength_ssp(smean_r, [], C*(nr_groups*(np-1)), np);
    end    

    %% Select points (samples)
    spoints = get_selection(signals.(aparams.signal), aparams.sel, ...
                            aparams.p1, aparams.p2);
                        
    %% Restrict coefficiencts only to selected points
    coef = coef(:, spoints);
    results.coef = coef;
    
    %% Obtain mean and covariance on selected points
    smean = smean_r(:,spoints);
    scov = C(spoints, spoints);
                        
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
    [M, B, W] = compute_ssp_e2_mmap(m_data_profile, metadata_profile, ...
                                    aparams.idx_traces, aparams.byte_sel);
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
    
    %% Set compression/selection parameters
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = U;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
    
    %% Obtain data for computing the coefficients    
    L = double(m_data_profile.data(1).X(:,idx_profile)');
    X = handle_prepare(L, pp1, pp2, pp3, pp4, pp5);
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance
    fprintf('Computing covariance...\n');
    if aparams.cov_from_sel
        fprintf('Computing data for covariance from selection...\n');
        x_cov = compute_features_e2_mmap(...
                    m_data_profile, metadata_profile, aparams.idx_traces, ...
                    handle_prepare, pp1, pp2, pp3, pp4, pp5, aparams.byte_sel);
        toc
        [~, C] = compute_template(x_cov);
        scov = mean(C, 3);
    else
        Z = X - get_leakage_from_map_base(D, coef, map_base);
        scov = (Z'*Z) / (np-1);
    end
elseif strcmp(atype, 'templatepca')
    %% Obtain data for computing the coefficients
    X = double(m_data_profile.data(1).X(:,idx_profile)');
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean_r = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute PCA params
    fprintf('Computing PCA parameters...\n');
    [U, ~, xmm, K] = compute_params_pca(smean_r, aparams.pca_threshold, ...
                                        aparams.pca_alternate);
    if isfield(aparams, 'pca_dimensions') && (aparams.pca_dimensions > 0)
        U = U(:,1:aparams.pca_dimensions);
    else
        U = U(:,1:K);
    end
    
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
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(Y, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors in PCA space
    smean = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance in PCA space
    fprintf('Computing covariance...\n');
    Z = Y - get_leakage_from_map_base(D, coef, map_base);
    scov = (Z'*Z) / (np-1);
elseif strcmp(atype, 'badpca')
    %% Obtain data for computing the coefficients
    X = double(m_data_profile.data(1).X(:,idx_profile)');
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean_r = get_leakage_from_map_base(data, coef, map_base);
    xmm = mean(smean_r, 1);
    
    %% Compute full covariance estimate
    fprintf('Computing full covariance...\n');
    Z = X - get_leakage_from_map_base(D, coef, map_base);
    C = (Z'*Z) / (np-1);
    
    %% Compute "bad" PCA params
    fprintf('Computing bad PCA parameters...\n');
    [U, ~, ~] = svd(C);    
    if isfield(aparams, 'pca_dimensions') && (aparams.pca_dimensions > 0)
        U = U(:,1:aparams.pca_dimensions);
    else
        U = U(:,1);
    end
    
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
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(Y, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors in PCA space
    smean = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance in PCA space
    fprintf('Computing covariance...\n');
    Z = Y - get_leakage_from_map_base(D, coef, map_base);
    scov = (Z'*Z) / (np-1);
elseif strcmp(atype, 'lda')
    %% Obtain the SSP matrices (M, B, W) on selected bytes
    fprintf('Obtaining sums of squares and cross products on selected bytes...\n');
    [M, B, W] = compute_ssp_e2_mmap(m_data_profile, metadata_profile, ...
                                    aparams.idx_traces, aparams.byte_sel);
    if isfield(eparams, 'save_ssp') && (eparams.save_ssp ~= 0)
        results.M = M;
        results.B = B;
        results.W = W;
    end
    xmm = mean(M, 1);
    toc
    
    %% Compute LDA params
    fprintf('Computing Fishers LDA parameters...\n');
    len_byte_sel = length(aparams.byte_sel);
    nr_traces_per_byte = length(aparams.idx_traces);
    Spool = W / (len_byte_sel*(nr_traces_per_byte-1));
    [A ,~, K] = compute_params_lda(B, Spool, len_byte_sel, aparams.lda_threshold);
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
    L = double(m_data_profile.data(1).X(:,idx_profile)');
    X = handle_prepare(L, pp1, pp2, pp3, pp4, pp5);
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance
    fprintf('Computing covariance...\n');
    if aparams.cov_from_sel
        fprintf('Computing data for covariance from selection...\n');
        x_cov = compute_features_e2_mmap(...
                    m_data_profile, metadata_profile, aparams.idx_traces, ...
                    handle_prepare, pp1, pp2, pp3, pp4, pp5, aparams.byte_sel);
        toc
        [~, C] = compute_template(x_cov);
        scov = mean(C, 3);
    else
        Z = X - get_leakage_from_map_base(D, coef, map_base);
        scov = (Z'*Z) / (np-1);
    end
elseif strcmp(atype, 'templatelda')
    %% Obtain data for computing the coefficients
    X = double(m_data_profile.data(1).X(:,idx_profile)');
    D = double(m_data_profile.data(1).B(2,idx_profile)');
    
    %% Select basis
    map_base = get_map_base(base);
    
    %% Compute Stochastic coefficients
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(X, D, map_base);
    toc
    
    %% Apprximate raw mean vectors from stochastic model
    data = 0:(metadata_profile.nr_groups-1);
    smean_r = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute raw covariance matrix
    fprintf('Computing raw covariance...\n');
    Z = X - get_leakage_from_map_base(D, coef, map_base);
    C = (Z'*Z) / (np-1);
    
    %% Compute raw between-groups matrix B
    fprintf('Computing between-groups matrix B...\n');
    xmm = mean(smean_r, 1);
    T = smean_r - ones(nr_groups, 1)*xmm;
    B = T'*T;
    
    %% Compute LDA params
    fprintf('Computing Fishers LDA parameters...\n');
    [A ,~, K] = compute_params_lda(B, C, nr_groups, aparams.lda_threshold);
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
    
    %% Project data using LDA
    Y = handle_prepare(X, pp1, pp2, pp3, pp4, pp5);
    
    %% Compute Stochastic coefficients in LDA space
    fprintf('Computing Stochastic coefficients...\n');
    [coef] = compute_coef_stochastic(Y, D, map_base);
    results.coef = coef;
    toc
    
    %% Apprximate mean vectors in LDA space
    smean = get_leakage_from_map_base(data, coef, map_base);
    
    %% Compute covariance in LDA space
    fprintf('Computing covariance...\n');
    Z = Y - get_leakage_from_map_base(D, coef, map_base);
    scov = (Z'*Z) / (np-1);
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
x_attack = compute_features_e2_mmap(...
                m_data_attack, metadata_attack, idx_attack, ...
                handle_prepare, pp1, pp2, pp3, pp4, pp5);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_attack = x_attack;
end
toc

%% Set evaluation parameters
tmiu = smean;
ic0 = inv(scov);
handle_eval = @evaluate_discriminant;
if strcmp(discriminant, 'linear')
    pe3 = tmiu;
    pe4 = ic0;
    pe5 = [];
    pe6 = [];
elseif strcmp(discriminant, 'linearnocov')
    pe3 = tmiu;
    pe4 = [];
    pe5 = [];
    pe6 = [];
else
    error('Unsupported discriminant type: %s', discriminant);
end

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.handle_eval = handle_eval;
    results.pe3 = pe3;
    results.pe4 = pe4;
    results.pe5 = pe5;
    results.pe6 = pe6;
end

%% Compute the success information
fprintf('Computing success info...\n');
[results.success_info] = get_success_info_like(x_attack, rand_iter, ...
                                       nr_traces_vec, ...
                                       handle_eval, pe3, pe4, pe5, pe6);

end
