function [M, B, W, np] = compute_ssp_e2_mmap_multi(s_data, bytes, pts, roffset)
%COMPUTE_SSP_E2_MMAP_MULTI Computes sums of squares and cross products
%   [M, B, W, np] = COMPUTE_SSP_E2_MMAP_MULTI(s_data, bytes, pts, roffset)
%   computes the sums of squares and cross products matrices from memory
%   mapped data, allowing to work with large data sets.
%
%   The sums of squares and cross product matrices are those matrices used
%   to compute the multivariate analysis of variance (MANOVA). They
%   represent the "treatment" matrix (B) and the "residual" or "error"
%   matrix (W). W can be seen as a scaled pooled covariance matrix, as it
%   is in fact the sum of the covariances from each data group.
%
%   Given B and W it is possible to compute the total sum of squares and
%   products corrected for the mean as:
%   B+W = sum_l*sum_j[(x_lj - x_bar)(x_lj-x_bar)']
%   where x_lj is the sample value of variable j (1<=j<=nr_points) for
%   group l (1<=l<=nr_groups), and x_bar is the overall mean vector.
%   The number of degrees of freedom of (B+W) is
%   nr_groups*np - 1,
%   where np is the total number of samples (from all sets) per each group.
%   See below for more details.
%   
%   With B, W and B+W it is then possible to compute Wilks' Lambda:
%   L = |W| / |B+W|
%   in order to test if the hypothesis H0 (there is no statistical
%   difference between the groups/treatments) can be rejected or not.
%   See the "Statistical Multivariate Analysis" book, p.302.
%
%   This method is targeted to the E2-related experiments, where a series
%   of bytes are loaded to some registers, and generally the value of a
%   particular byte is being analysed while the values of the other bytes
%   are either fixed or random.
%
%   This method is very similar to compute_ssp_e2_mmap, but allows to
%   combine data from different sets (i.e. different memory maped files).
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
%   bytes should be a vector of byte values that will select which groups
%   to consider in order to compute the covariance matrix. Each index in
%   the bytes vector should be between 0 and (nr_groups-1).
%   Note: group 1 corresponds to byte value 0, group 2 to byte value
%   1, ..., group 256 to byte value 255.
%
%   pts should be a vector of indeces that specifies which points from each
%   trace will be used to create SSP matrices. If not given or empty then
%   the default from metadata (1:metadata.nr_points) will be used.
%
%   roffset is an optional parameter that can be used to add a random
%   offset to each trace. This should be a cell of length nr_sets, where
%   each element is a matrix of size np_set x nr_groups, containing the
%   random offset to be added to each trace and group per set. If not
%   necessary pass [] (empty) or ignore this parameter.
%
%   The outputs of this method are:
%   - M: a matrix of size nr_groups x nr_points, having the mean vector
%   for each group specified by the "bytes" parameter.
%   - B: the "treatment" matrix, containing the variance caused only by the
%   mean vectors of each group. The number of degrees of freedom of B is
%   df_B = nr_groups - 1
%   where nr_groups is given by the first dimension of M.
%   - W: the "residual" matrix, containing the pooled (sumed) covariance
%   of the data within each group. The number of degrees of freedom of W is
%   df_W = nr_groups x (np - 1)
%   Both B and W are matrices of size nr_points x nr_points.
%   
%   Note: W can be used as a "pooled" covariance matrix, i.e. a better
%   estimate of a covariance matrix, when the covariance of each group is
%   similar. However, in this case, W should be divided by the number of
%   degrees of freedom df_W.
%
%   See also compute_ssp_mmap, compute_features_e2_mmap, get_mmap.

%% Check and initialize data
nr_sets = s_data.nr_sets;
nr_groups = s_data.metadata{1}.nr_groups;
for k=2:nr_sets
    if s_data.metadata{2}.nr_groups ~= nr_groups
        error('Set %d has different number of groups than set 1: ', k);
    end
end
nr_bytes = length(bytes);
if min(bytes) < 0 || max(bytes) > (nr_groups-1)
    error('Some index in bytes is outside available groups');
end
if nargin < 3 || isempty(pts)
    nr_points = s_data.metadata{1}.nr_points;
    pts = 1:nr_points;
else
    pts = pts(:);
    nr_points = length(pts);
end
if nargin < 4
    roffset = [];
end
M = zeros(nr_bytes, nr_points);
B = zeros(nr_points, nr_points);
W = zeros(nr_points, nr_points);
np = 0;
for k=1:nr_sets
    np = np + length(s_data.idx{k});
end
fprintf('Running compute_ssp_e2_mmap_multi()...\n');

%% Compute the group means and the residual matrix W
fprintf('Computing the mean vectors M and the residual matrix W...\n');
for k=1:nr_bytes
    idx = 1;
    L = zeros(np, nr_points);
    for i=1:nr_sets
        kindex = find(s_data.mmap_data{i}.data(1).B(2,:)==bytes(k));
        lindex = kindex(s_data.idx{i});
        % Below note transpose and conversion in case we had integer class
        L(idx:idx+length(lindex)-1,:) = double(s_data.mmap_data{i}.data(1).X(pts,lindex)');
        if ~isempty(roffset)
            L(idx:idx+length(lindex)-1,:) = ...
                L(idx:idx+length(lindex)-1,:) + roffset{i}(:,bytes(k)+1)*ones(1,nr_points);
        end
        idx = idx + length(lindex);
    end    
    M(k,:) = mean(L, 1);
    X = L - ones(np, 1)*M(k,:);
    W = W + X'*X;
end
mvec = mean(M, 1);

%% Compute the treatment matrix B
fprintf('Computing the treatment matrix B...\n');
for k=1:nr_bytes
    xm = M(k,:) - mvec;
    B = B + xm'*xm;
end
B = B * np;
                              
end
    
