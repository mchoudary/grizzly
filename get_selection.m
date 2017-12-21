function [points] = get_selection(curve, selection, p1)
%GET_SELECTION Returns a vector of interesting points for template attacks
%   [points] = GET_SELECTION(curve, selection)
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
%   Choose the selection method from:
%   - '1ppc': 1 peak per clock at most, and must be above 95% percentile.
%   p1 selects minimum peak distance; this should be slightly less than
%   a clock period.
%   - '3ppc': same as '1ppc' but taking also the left and right neighbor
%   points of each peak. p1 selects minimum peak distance.
%   - '20ppc': 20 points per clock at most, and must be above 95%
%   percentile. p1 selects minimum peak distance; this should be slightly
%   less than a clock period.
%   - 'allam': all points above mean curve. p1 has no effect.
%   - 'allap': all points above given percentile. p1 selects percentile,
%   e.g. give p1=0.95 to select all points above the 95% percentile.
%
%   This method returns a vector of points for the given selection.

%% Check and initialize stuff
curve = curve(:);

%% Compute selection
if strcmp(selection, '1ppc')
    scurve = sort(curve);
    idx = floor(0.95*length(curve));
    [~, points] = findpeaks(curve, 'minpeakdistance', p1, ...
                         'minpeakheight', scurve(idx));
elseif strcmp(selection, '3ppc')
    scurve = sort(curve);
    idx = floor(0.95*length(curve));
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
    idx = floor(0.95*length(curve));
    [~, pks] = findpeaks(curve, 'minpeakdistance', p1, ...
                         'minpeakheight', scurve(idx));
    points = [];
    for k=1:length(pks)
        start = pks(k) - floor(p1/2);
        if start <= 0
            ival = 1:start+p1;
            start = 1;
        else
            ival = start:start+p1;
        end
        [~, spoints] = sort(curve(ival), 1, 'descend');
        for i=1:20
            if curve(spoints(i) + start - 1) > scurve(idx)
                points = [points; spoints(i) + start - 1];
            end
        end
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
end

end

