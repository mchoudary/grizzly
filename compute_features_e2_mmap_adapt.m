function [xdata] = compute_features_e2_mmap_adapt(map, metadata, bindex, ...
                                    s_adapt, ...
                                    func_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5, bytes)
%COMPUTE_FEATURES_E2_MMAP_ADAPT Computes features from sample vectors
%   [xdata] = COMPUTE_FEATURES_E2_MMAP_ADAPT(map, metadata, bindex, ...
%                              s_adapt, ...
%                              func_prepare, ...
%                              pp1, pp2, pp3, pp4, pp5, bytes)
%   computes some features from sample vectors. This method is targeted
%   to the E2-related experiments, where a series of bytes are loaded to
%   some registers, and generally the value of a particular byte is being
%   analysed while the values of the other bytes are either fixed or
%   random.
%
%   This method is very similar to compute_features_e2_mmap, but it first
%   adapts the data for each group using the data in s_adapt and then
%   applies func_prepare.
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
%   s_adapt should be a structure with adaptation data. The mandatory field
%   is 'type', which specifies the type of adaptation. Based on 'type'
%   different parameters are expected.
%   The current types and parameters supported are:
%   - 'none': no adaptation (has no parameters).
%   - 'offset_median': This type of adaptation will simply compute an
%   offset between the median of each data trace and the median of the
%   xmm parameter and substract that from the data trace.
%   This type has the following parameters:
%       - 'xmm': this should be the overall mean trace from the profiling
%       data.
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
fprintf('Running compute_features_e2_mmap_adapt()...\n');

%% Extract data from each group leakage matrix
for k=1:nr_groups
    fprintf('Computing features for data group %d\n', bytes(k));
    kindex = find(map.data(1).B(2,:)==bytes(k));
    lindex = kindex(bindex);
    % Below note transpose and conversion in case we had integer class
    L = double(map.data(1).X(:,lindex)');
    
    if strcmp(s_adapt.type, 'none')
        % Do nothing
    elseif strcmp(s_adapt.type, 'offset_median')
        xmm = s_adapt.xmm;
        nr_traces = size(L, 1);
        for i=1:nr_traces
            offset = median(L(i,:)) - median(xmm);
            L(i,:) = L(i,:) - offset;
        end
    else
        error('Adaptation type not recognized: %s', s_adapt.type);
    end
    
    data_out = func_prepare(L, pp1, pp2, pp3, pp4, pp5);    
    if k==1        
        [nr_traces, nr_features] = size(data_out);
        xdata = zeros(nr_traces, nr_features, nr_groups);
    end    
    xdata(:,:,k) = data_out;
end
                              
end
    
