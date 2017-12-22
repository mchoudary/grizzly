function [points] = get_selection(curve, selection, p1, p2)
%GET_SELECTION Returns a vector of interesting points for template attacks
%   [points] = GET_SELECTION(curve, selection, p1, p2)
%   returns a vector of interesting points for template attacks based on
%   the signal curve and selection parameters.
%
%   curve should be a vector of length nr_points containing the signal
%   curve from which to select points.
%
%   selection should be a string defining the selection method to be used.
%
%   p1 should be a parameter for the given selection method, e.g. the
%   number of samples distance for 1ppc.
%
%   p2 is an additional but optional parameter which might be used to
%   change some default values, such as the minimum percentile for 1ppc.
%   
%   Choose the selection method from:
%   - '1ppc': 1 peak per clock at most, and must be above some percentile.
%   p1 selects minimum peak distance; this should be slightly less than
%   a clock period. p2 selects percentile, should be in [0,1] (default 0.95).
%   - '3ppc': same as '1ppc' but taking also the left and right neighbor
%   points of each peak.
%   - '20ppc': 20 points per clock at most, and must be above 95%
%   percentile.
%   - 'allam': all points above mean curve. p1 has no effect.
%   - 'allap': all points above given percentile. p1 selects percentile,
%   e.g. give p1=0.95 to select all points above the 95% percentile.
%   - 'all: this returns all points, basically ignoring the curve.
%
%   This method returns a vector of points for the given selection.

%% Check and initialize stuff
curve = curve(:);
if nargin < 4
    p2 = [];
end

%% Compute selection
if strcmp(selection, '1ppc')
    scurve = sort(curve);
    if isempty(p2)
        p2 = 0.95;
    end
    idx = floor(p2*length(curve));
    [~, points] = findpeaks(curve, 'minpeakdistance', p1, ...
                         'minpeakheight', scurve(idx));
elseif strcmp(selection, '3ppc')
    scurve = sort(curve);
    if isempty(p2)
        p2 = 0.95;
    end
    idx = floor(p2*length(curve));
    [~, pstart] = findpeaks(curve, 'minpeakdistance', p1, ...
                         'minpeakheight', scurve(idx));
    points = [];
    for p=1:length(pstart)
        if pstart(p) > 1
            points = [points; pstart(p)-1];
        end;
        points = [points; pstart(p)];
        if pstart(p) < length(curve)
            points = [points; pstart(p)+1];
        end
    end
elseif strcmp(selection, '20ppc')
    scurve = sort(curve);
    if isempty(p2)
        p2 = 0.95;
    end
    idx = floor(p2*length(curve));
    [~, pks] = findpeaks(curve, 'minpeakdistance', p1, ...
                         'minpeakheight', scurve(idx));
    points = [];
    offset = 1;
    if (pks(1) - floor(p1/2)) <= 0
        offset = offset + abs(pks(1) - floor(p1/2));
    end
    for k=1:length(pks)
        start = pks(k) - floor(p1/2) + offset;
        ival = start:start+p1;
        ival = ival(ival <= length(curve));
        [~, spoints] = sort(curve(ival), 1, 'descend');
        pset = [];
        for i=1:min(length(ival), 20)
            if curve(spoints(i) + start - 1) > scurve(idx)
                pset = [pset; spoints(i) + start - 1];
            end
        end
        pset = sort(pset);
        points = [points; pset(:)];
    end        
elseif strcmp(selection, 'allam')
    mcurve = mean(curve);
    points = find(curve > mcurve);
elseif strcmp(selection, 'allap')
    if p1 < 0 || p1 >= 1
        error('Wrong percentile number, must be between 0 and 1');
    end
    scurve = sort(curve);
    idx = floor(p1*length(curve));
    points = find(curve > scurve(idx));    
elseif strcmp(selection, 'all')
    points = 1:length(curve);
end

end

