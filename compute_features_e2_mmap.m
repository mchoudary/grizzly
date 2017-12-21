function [xdata] = compute_features_e2_mmap(map, metadata, bindex, ...
                                    func_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5, bytes)
%COMPUTE_FEATURES_E2_MMAP Computes features from sample vectors
%   [xdata] = COMPUTE_FEATURES_E2_MMAP(map, metadata, bindex, ...
%                              func_prepare, ...
%                              pp1, pp2, pp3, pp4, pp5, bytes)
%   computes some features from sample vectors. This method is targeted
%   to the E2-related experiments, where a series of bytes are loaded to
%   some registers, and generally the value of a particular byte is being
%   analysed while the values of the other bytes are either fixed or
%   random.
%
%   This method will apply the func_prepare method (details below) on the
%   data represented by the memmapfile object map, in order to compute some
%   features from each sample (trace) in the input data. In addition, this
%   method will split the data into the different group categories
%   available, such that the output of this method will be a matrix xdata
%   of size nr_rows x nr_features x nr_groups.
%
%   map should be a memmapfile object with 1 data entry containing at least:
%   - transposed data matrix X of size nr_points x nr_trials
%   - transposed input bytes matrix B of size nr_bytes x nr_trials
%
%   metadata should provide in addition at least the number of groups
%   (nr_groups) for which different data was acquired. The assumption here
%   is that the data in the matrix X is from random input byte values, as
%   specified by the matrix B, and that the number of different byte values
%   is given by nr_groups.
%
%   bindex should be a vector of indeces used to select which block samples
%   will be used to compute the signal curves. The data in X is assumed to
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
%   The output xdata is a matrix of size nr_rows x nr_features x nr_groups,
%   where nr_rows is the number of traces per group,
%   nr_features is the number of features remaining per trace after
%   applying the func_prepare method and nr_groups is the number of groups
%   specified in the metadata.
%
%   See also compute_features, get_mmap.

%% Check and initialize data
if nargin > 9 && ~isempty(bytes)
    nr_groups = length(bytes);
else
    nr_groups = metadata.nr_groups;
    bytes = 0:(nr_groups-1);
end

%% Extract data from each group leakage matrix
for k=1:nr_groups
    fprintf('Computing features for data group %d\n', bytes(k));
    kindex = find(map.data(1).B(2,:)==bytes(k));
    lindex = kindex(bindex);
    % Below note transpose and conversion in case we had integer class
    L = double(map.data(1).X(:,lindex)');
    
    data_out = func_prepare(L, pp1, pp2, pp3, pp4, pp5);    
    if k==1        
        [nr_traces, nr_features] = size(data_out);
        xdata = zeros(nr_traces, nr_features, nr_groups);
    end    
    xdata(:,:,k) = data_out;
end
                              
end
    
