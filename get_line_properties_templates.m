function [slines] = get_line_properties_templates(uid, style)
%GET_LINE_PROPERTIES_TEMPLATES Returns line properties for plots
%   [slines] = GET_LINE_PROPERTIES_TEMPLATES(uid, style)
%   returns line properties for plots of template attack results.
%
%   This method may be useful to get consitent line properties for each
%   compression method and parameters so that plots with same parameters
%   may use the same color, width, etc.
%
%   uid should be an integer that uniquely identifies the desired line
%   properties. This could be used to make sure that the same properties
%   are used across different figures. This method currently supports only
%   the integers from 1 to 6.
%
%   slines is a structure containing the following information:
%   - 'Color'
%   - 'LineStyle'
%   - 'LineWidth'
%   - 'Marker'
%
%   style is a string that specifies a kind of style, that may change the
%   properties returned. Currently supported styles are:
%   - 'normal'
%   - 'fancy'
%
%   The slines structure can be used with make_figures_ge to use the
%   selected properties.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

if strcmp(style, 'normal')
    if uid == 1
        slines.Color = 'm';
        slines.LineStyle = '-';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    elseif uid == 2
        slines.Color = 'c';
        slines.LineStyle = '-.';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    elseif uid == 3        
        slines.Color = 'g';
        slines.LineStyle = '-';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    elseif uid == 4
        slines.Color = 'k';
        slines.LineStyle = '-.';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    elseif uid == 5
        slines.Color = 'b';
        slines.LineStyle = '-';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    elseif uid == 6
        slines.Color = 'r';
        slines.LineStyle = '-.';
        slines.LineWidth = 4;
        slines.Marker = 'none';
    else
        error('uid not supported');
    end
elseif strcmp(style, 'fancy')
    if uid == 1
        slines.Color = 'm';
        slines.LineStyle = '-';
        slines.LineWidth = 2;
        slines.Marker = 'o';
    elseif uid == 2
        slines.Color = 'c';
        slines.LineStyle = '--';
        slines.LineWidth = 2;
        slines.Marker = '+';
    elseif uid == 3
        slines.Color = 'g';
        slines.LineStyle = '-.';
        slines.LineWidth = 2;
        slines.Marker = '*';
    elseif uid == 4
        slines.Color = 'k';
        slines.LineStyle = '-';
        slines.LineWidth = 2;
        slines.Marker = '.';
    elseif uid == 5
        slines.Color = 'b';
        slines.LineStyle = '--';
        slines.LineWidth = 2;
        slines.Marker = 'x';
    elseif uid == 6
        slines.Color = 'r';
        slines.LineStyle = '-.';
        slines.LineWidth = 2;
        slines.Marker = 's';
    else
        error('uid not supported');
    end
else
    error('Unknown style');
end

end