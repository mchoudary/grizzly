function [coef] = compute_coef_stochastic(X, D, map_base)
%COMPUTE_COEF_STOCHASTIC Computes the stochastic coeficients beta
%   [coef] = COMPUTE_COEF_STOCHASTIC(X, D, map_base)
%   computes beta coeficients for a stochastic model based on the given
%   leakage, input data and mappings.
%
%   X should be a matrix of size nr_traces x nr_points containing the
%   leakage data on some set of points.
%   
%   D should be a vector of length nr_traces, containing the data
%   corresponding to X.
%
%   map_base should be either a cell of size u x 1 if the same mapping
%   bases should be used on all the nr_points from tdata, or a cell of size
%   u x nr_points if different mapping bases should be used for each point
%   in each trace from tdata. Each element of map_base should be a function
%   handle of the following form:
%       map_base{j, t} = func(v)
%   returning the real value of the mapping for the value v
%   at index j and time t.
%
%   This method returns a matrix of coefficients of size u x nr_points
%   that should be used with the given mapping bases in order to model the
%   given leakage using the stochastic appraoch. Each column of
%   coefficients should be used with the particular point to which it
%   corresponds.
%
%   See the paper "A Stochastic Model for Differential Side Channel
%   Cryptanalysis", by Schindler et al.

%% Check and initialize data
nr_traces = size(X, 1);
nr_points = size(X, 2);
D = D(:);
if length(D) ~= nr_traces
    error('Wrong size of xdata');
end
u = size(map_base, 1);
coef = zeros(u, nr_points);

if size(map_base, 2) == 1
    %% Compute coefficients using same base for all points (fast)
    A = zeros(nr_traces, u);
    for k = 1:nr_traces
        for i=1:u
            A(k,i) = map_base{i}(D(k));
        end
    end
    M = (A'*A)\A';
    coef = M * X;    
else
    %% Compute coeficients for each point using individual bases
    if size(map_base, 2) ~= nr_points
        error('Wrong size of map_base');
    end
    
    for j = 1:nr_points
        fprintf('Computing coeficients for j=%d\n', j);
        x = X(:,j);
        A = zeros(nr_traces, u);
        for k = 1:nr_traces
            for i=1:u
                A(k,i) = map_base{i,j}(D(k));
            end
        end
        coef(:,j) = ((A'*A)\A') * x;
    end
end

end



