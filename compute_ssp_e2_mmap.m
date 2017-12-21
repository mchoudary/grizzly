function [M, B, W] = compute_ssp_e2_mmap(map, metadata, bindex, bytes, pts)
%COMPUTE_SSP_E2_MMAP Computes sums of squares and cross products
%   [M, B, W] = COMPUTE_SSP_E2_MMAP(map, metadata, bindex, bytes, pts)
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
%   nr_groups*nr_samples_per_groups - 1. See below for more details.
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
%   will be used to compute the covariance matrix. The data in X is assumed to
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
%   The outputs of this method are:
%   - M: a matrix of size nr_groups x nr_points, having the mean vector
%   for each group specified by the "bytes" parameter.
%   - B: the "treatment" matrix, containing the variance caused only by the
%   mean vectors of each group. The number of degrees of freedom of B is
%   df_B = nr_groups - 1
%   where nr_groups is given by the first dimension of M.
%   - W: the "residual" matrix, containing the pooled (sumed) covariance
%   of the data within each group. The number of degrees of freedom of W is
%   df_W = nr_groups x (nr_samples_per_group - 1)
%   where nr_samples_per_group is given by the length of the supplied
%   parameter "bindex".
%   Both B and W are matrices of size nr_points x nr_points.
%   
%   Note: W can be used as a "pooled" covariance matrix, i.e. a better
%   estimate of a covariance matrix, when the covariance of each group is
%   similar. However, in this case, W should be divided by the number of
%   degrees of freedom df_W = nr_groups x (nr_samples_per_group - 1).
%
%   See also compute_features_e2_mmap, get_mmap.

%% Check and initialize data
nr_groups = metadata.nr_groups;
nr_bytes = length(bytes);
if min(bytes) < 0 || max(bytes) > (nr_groups-1)
    error('Some index in bytes is outside available groups');
end
if nargin < 5 || isempty(pts)
    nr_points = metadata.nr_points;
    pts = 1:nr_points;
else
    pts = pts(:);
    nr_points = length(pts);
end
M = zeros(nr_groups, nr_points);
B = zeros(nr_points, nr_points);
W = zeros(nr_points, nr_points);
nr_samples_per_group = length(bindex);

%% Compute the group means and the residual matrix W
fprintf('Computing the mean vectors M and the residual matrix W...\n');
for k=1:nr_bytes
    kindex = find(map.data(1).B(2,:)==bytes(k));
    lindex = kindex(bindex);
    % Below note transpose and conversion in case we had integer class
    L = double(map.data(1).X(pts,lindex)');
    M(k,:) = mean(L, 1);
    X = L - ones(nr_samples_per_group, 1)*M(k,:);
    W = W + X'*X;
end
mvec = mean(M, 1);

%% Compute the treatment matrix B
fprintf('Computing the treatment matrix B...\n');
for k=1:nr_bytes
    xm = M(k,:) - mvec;
    B = B + xm'*xm;
end
B = B * nr_samples_per_group;
                              
end
    
