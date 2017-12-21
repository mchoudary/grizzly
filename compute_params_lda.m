function [A, D] = compute_params_lda(B, S)
%COMPUTE_PARAMS_LDA Compute Fisher's LDA parameters
%   [A, D] = COMPUTE_PARAMS_LDA(B, S)
%   computes the parameters A (matrix of coefficients/eigenvectors) and
%   D (diagonal matrix of eigenvalues of S^(-1)*B) that are needed to
%   obtain Fisher's linear discriminants.
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

%% Check and initialize stuff
nr_points = size(B, 1);
if size(B,2) ~= nr_points
    error('B not square');
end
if size(B) ~= size(S)
    error('size of B different than size of S');
end

%% Compute the eigenvalues
[A, D, ~] = svd(S\B);

%% Scale eigenvalues to have e'Se = 1 for each e
Q = A'*S*A;
Z = diag(1./sqrt(diag(Q)));
A = A*Z;

end