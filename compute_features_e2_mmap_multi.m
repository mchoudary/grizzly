function [xdata] = compute_features_e2_mmap_multi(s_data, ...
                                    func_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5, bytes, roffset)
%COMPUTE_FEATURES_E2_MMAP_MULTI Computes features from sample vectors
%   [xdata] = COMPUTE_FEATURES_E2_MMAP_MULTI(s_data, ...
%                              func_prepare, ...
%                              pp1, pp2, pp3, pp4, pp5, bytes, roffset)
%   computes some features from sample vectors. This method is targeted
%   to the E2-related experiments, where a series of bytes are loaded to
%   some registers, and generally the value of a particular byte is being
%   analysed while the values of the other bytes are either fixed or
%   random.
%
%   This method is very similar to compute_features_e2_mmap but extracts the
%   features from multiple sets (i.e. memory mapped objects).
%
%   This method will apply the func_prepare method (details below) on the
%   given data, in order to compute some
%   features from each sample (trace) in the input data. In addition, this
%   method will split the data into the different group categories
%   available, such that the output of this method will be a matrix xdata
%   of size nr_rows x nr_features x nr_groups.
%
%   s_data should be a structure containing the data from which to compute
%   the sums of squares and cross-product matrices. This structure should
%   contain the following fields:
%   - 'nr_sets': the number of sets (i.e. memory mapped objects) available
%   - 'mmap_data': a cell of length nr_sets, having the memmapfile objects, which
%   should contain at least:
%       - transposed data matrix X of size nr_points x nr_trials
%       - transposed input bytes matrix B of size nr_bytes x nr_trials
%   - 'metadata': a cell of length nr_sets, having the metadata structures.
%   These should provide at least the number of groups
%   (nr_groups) for which different data was acquired. The assumption here
%   is that the data in the matrix X is from random input byte values, as
%   specified by the matrix B, and that the number of different byte values
%   is given by nr_groups. Note that this method expects that all the data
%   sets contain data for all the groups, even if the sample size per group
%   is different across sets.
%   - 'idx': a cell of length nr_sets, having the vectors of indeces
%   used to select which block samples will be used to compute the SSP matrices.
%   The data in X is assumed to
%   be taken over nr_trials trials, where each consecutive nr_groups trials
%   form a block. In each block all the possible values of the target byte
%   should be covered. Therefore, when selecting samples corresponding to a
%   particular byte value, this method will select one sample (trial) from
%   each block, for a total of nr_blocks samples. The bindex vector should
%   provide a filter to select only the specified samples, allowing to test
%   different selection sizes. The values in bindex should be between 1 and
%   nr_blocks inclusive. nr_blocks can be computed as nr_trials/nr_groups
%   from the metadata information.
%
%   The func_prepare method should transform the data (e.g. compress it)
%   for the further processing. The output data can then be used to test
%   leakage, classification, etc. An exmaple of func_prepare method is
%   prepare_data_pca which uses some PCA parameters to retain only some
%   principal directions of the data. The prototype of func_prepare is:
%   data_out = func_prepare(data_in, pp1, pp2, pp3, ...) where data_in is
%   the raw input trace and pp1, pp2, pp3 are parameters needed by the
%   specific function.
%
%   bytes is an optional vector parameter which specifies the indices of
%   the bytes for which to compute the features, e.g. only a subset of all
%   the target bytes for which data is available.
%
%   roffset is an optional parameter that can be used to add a random
%   offset to each trace. This should be a cell of length nr_sets, where
%   each element is a matrix of size np_set x nr_groups, containing the
%   random offset to be added to each trace and group per set. If not
%   necessary pass [] (empty) or ignore this parameter.
%
%   The output xdata is a matrix of size nr_rows x nr_features x nr_groups,
%   where nr_rows is the total (across all sets) number of traces per group,
%   nr_features is the number of features remaining per trace after
%   applying the func_prepare method and nr_groups is the number of groups
%   specified in the metadata.
%
%   See also compute_features, compute_features_e2_mmap, get_mmap.

%% Check and initialize data
nr_sets = s_data.nr_sets;
if nargin > 7 && ~isempty(bytes)
    nr_groups = length(bytes);
else
    nr_groups = s_data.metadata{1}.nr_groups;
    bytes = 0:(nr_groups-1);
end
if min(bytes) < 0 || max(bytes) > (nr_groups-1)
    error('Some index in bytes is outside available groups');
end
np = 0;
for k=1:nr_sets
    np = np + length(s_data.idx{k});
end
nr_points = s_data.metadata{1}.nr_points;
if nargin < 9
    roffset = [];
end
fprintf('Running compute_features_e2_mmap_multi()...\n');

%% Extract data from each group leakage matrix
for k=1:nr_groups
    idx = 1;
    L = zeros(np, nr_points);
    for i=1:nr_sets
        kindex = find(s_data.mmap_data{i}.data(1).B(2,:)==bytes(k));
        lindex = kindex(s_data.idx{i});
        % Below note transpose and conversion in case we had integer class
        L(idx:idx+length(lindex)-1,:) = double(s_data.mmap_data{i}.data(1).X(:,lindex)');
        if ~isempty(roffset)
            L(idx:idx+length(lindex)-1,:) = ...
                L(idx:idx+length(lindex)-1,:) + roffset{i}(:,bytes(k)+1)*ones(1,nr_points);
        end
        idx = idx + length(lindex);
    end
    
    data_out = func_prepare(L, pp1, pp2, pp3, pp4, pp5);    
    if k==1        
        [nr_traces, nr_features] = size(data_out);
        xdata = zeros(nr_traces, nr_features, nr_groups);
    end    
    xdata(:,:,k) = data_out;
end
                              
end
    
