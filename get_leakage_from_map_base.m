function [X] = get_leakage_from_map_base(D, coef, map_base)
%GET_LEAKAGE_FROM_MAP_BASE Returns leakage traces using base mappings.
%   [X] = GET_LEAKAGE_FROM_MAP_BASE(D, coef, map_base)
%   returns leakage traces X from the input values D, stochastic
%   coefficients coef and the base mappings map_base.
%
%   D should be a vector of length nr_traces, containing input data, for
%   which to return corresponding leakage traces.
%
%   coef should be a matrix of size nr_bases x nr_points that will be used
%   together with map_base in order to produce traces of length nr_points.
%
%   map_base should be either a cell of size u x 1 if the same mapping
%   bases should be used on all the nr_points from tdata, or a cell of size
%   u x nr_points if different mapping bases should be used for each point
%   in each trace from tdata. Each element of map_base should be a function
%   handle of the following form:
%       map_base{k, t} = func(v)
%   returning the real value of the mapping for the value v
%   at index k and time t.
%
%   This method returns the matrix X of size nr_traces x nr_points, having
%   nr_traces of nr_points each.

%% Check and initialise parameters
nr_traces = length(D);
[nr_bases, nr_points] = size(coef);
if size(map_base, 1) ~= nr_bases
    error('Incompatible sizes between coef and map_base');
end

if size(map_base, 2) == 1
    %% Compute using same base (faster)
    V = zeros(nr_traces, nr_bases);
    for i=1:nr_traces
        for k=1:nr_bases
            V(i,k) = map_base{k}(D(i));
        end
    end

    X = V * coef;
else
    %% Compute using individual bases (slower)
    if size(map_base, 2) ~= nr_bases
        error('wrong size of map_base');
    end
    
    X = zeros(nr_traces, nr_points);
    for i=1:nr_traces
        for j=1:nr_points
            v = zeros(1, nr_bases);
            for k=1:nr_bases
                v(k) = map_base{k,j}(D(i));
            end
            X(i,j) = v * coef(:,j);
        end
    end
    
end

end

