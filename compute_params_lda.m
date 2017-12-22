function [A, D, K] = compute_params_lda(B, S, nr_groups, threshold, use_eig)
%COMPUTE_PARAMS_LDA Compute Fisher's LDA parameters
%   [A, D, K] = COMPUTE_PARAMS_LDA(B, S, nr_groups, threshold)
%   computes the parameters A (matrix of coefficients/eigenvectors) and
%   D (vector of eigenvalues of S^(-1)*B) that are needed to
%   obtain Fisher's linear discriminants. Optionally the method returns
%   the number K of coefficients needed to retain the given threshold of
%   the total variance.
%
%   B should be "between means" sum of squares and crossproducts matrix:
%   B = Sum((x_i_mean - x_mean)(x_i_mean - x_mean), where x_i_mean is the
%   mean vector of group i and x_mean is the overall mean vector.
%   Note that it does not matter much if B is scaled in some manner, the
%   only difference will be that the returned eigenvalues will be scaled
%   accordingly.
%
%   S should be the "pooled" pooled covariance matrix, i.e. the scaled
%   "within groups" sum of squares and crossproducts matrix:
%   S = (Sum_i Sum_j (x_i_j - x_i_mean)(x_i_j - x_i_mean)) /  ...
%       (nr_groups * (nr_samples_per_group - 1))
%   where i goes over all groups and j over all observations in a given
%   group.
%
%   Both B and S should be of size nr_points x nr_points.
%
%   The following 2 arguments are optional and will be used only to compute
%   K, the number of components needed to reach a certain threshold of the
%   total variance:
%   - nr_groups should specify how many groups/classes where used to
%     compute B and S.
%   - threshold is an optional real value (between 0 and 1) that should
%     specify the minimum threshold of the total variance.
%
%   use_eig is an optional argument that specifies if the actual
%   eigenvalues should be used instead of the singular values (default).
%   Pass a non-zero value for using the eigenvalues. Otherwise this method
%   uses in fact a singular value decomposition, which is more stable even
%   if the result does not represent precisely the eigenvectors.
%
%   It is important that the given matrix S is already scaled according to
%   the degrees of freedom as shown above, since the eigenvectors e
%   returned in matrix A are such that e'*S*e=1, and this in turn will
%   help in classification applications, since the covariance of
%   discriminants of the form y1=e1'*x, y2=e2'*x, ... will be the identity
%   (I).
%
%   Please see the book "Applied Multivariate Statistical Analysis",
%   Section 11.6 for more details.
%
%   This method returns:
%   - the eigenvector matrix A having the scaled eigenvectors e of S^(-1)*B
%   such that e'*S*e = 1. These eigenvectors, also known as coefficients,
%   can be used to compute the linear discriminants y = e'*x.
%   This matrix can be used with methods such as prepare_data_pca_v2.
%   A is of size nr_points x nr_points, having the scaled eigenvectors as
%   columns in A. Note however that only the first
%   s <= min(nr_points, nr_groups-1) eigenvectors are non-zero, the rest
%   can be ignored.
%   - the diagonal matrix D also of size nr_points x nr_points,
%   containing the eigenvalues of S^(-1)*B.
%   - an optional number of components K that represent the number of
%     eigenvalues needed to reach the specified threshold of the total
%     variance. This is the cummulative variance method of determining the
%     number of components to use.
%
%   See also prepare_data_pca_v1, prepare_data_pca_v2.

%% Check and initialize stuff
nr_points = size(B, 1);
if size(B,2) ~= nr_points
    error('B not square');
end
if size(B) ~= size(S)
    error('size of B different than size of S');
end
if nargin < 5
    use_eig = 0;
end

%% Compute the eigenvalues
if use_eig
    [A,D] = eig(S\B);
else
    [A, D, ~] = svd(S\B);
end
D = diag(D);

%% Scale eigenvalues to have e'Se = 1 for each e
Q = A'*S*A;
Z = diag(1./sqrt(diag(Q)));
A = A*Z;

%% Return K if needed
if nargout > 2
    for k=1:nr_groups
        f = sum(D(1:k)) / sum(D);
        if f >= threshold
            K = k;
            break
        end
    end
end

end