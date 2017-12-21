function make_figures_ge(G, nr_traces_vec, ...
                         rpath, rprefix, ...
                         title_ge, slegend, font_size, ...
                         slines, options, yrange)
%MAKE_FIGURES_GE creates figures for guessing entropy data
%   MAKE_FIGURES_GE(G, nr_traces_vec, ...
%                        rpath, rprefix, ...
%                        title_ge, slegend, font_size, ...
%                        slines, options, yrange)
%   creates figures for guessing entropy data.
%
%   G must be a matrix of size nr_exp x nr_test_groups, where nr_exp
%   specifies the number of experiments and nr_test_groups is the
%   number of different iterations (x axis) for which the guessing entropy
%   was computed at each experiment.
%
%   As an example, the success_info structure returned by
%   get_success_info_like can be used to produce guessing entropy values
%   using the function get_ge_from_success_info.
%
%   nr_traces_vec is a vector containing the number of test traces used for
%   the guessing entropy in each experiment. This will be used for the
%   x-axis in figures.
%
%   rpath is a string containing the path where the figures should be
%   saved.
%
%   rprefix is a string that will be used to prefix all saved figures.
%
%   title_result specifies an overall title for the figures. Pass []
%   (empty) to ignore.
%
%   slegend should be a cell of strings containing the legend to be shown
%   for each experiment. Pass [] to omit legends.
%
%   font_size specifies the font size of the text. Pass [] to use the
%   default (24).
%
%   slines should be a cell of structures of length nr_exp, containing line
%   properties to be used with each experiment. Each structure in slines
%   should contain the following fields:
%   - 'Color'
%   - 'LineStyle'
%   - 'LineWidth'
%   - 'Marker'
%   See "doc plot" for details on these parameters.
%
%   options is a string that can specify options for the plot.
%   Pass one or more of the following values if desired:
%       'y': use the default Matlab ylim values for the given data (to
%       maximize the space over which the data is shown). If this is not
%       given then the range specified by yrange (see below) is used.
%       'g': use grid on. Default is off.
%       'L': use large figure size (useful with 'Publish' mode).
%
%   yrange is a 2-element vector that specifies the ylim values for the
%   figures.

%% Initialize and check stuff
nr_exp = size(G, 1);
if isempty(font_size)
    font_size = 24;
end
nr_test_groups = length(nr_traces_vec);
if size(G, 2) ~= nr_test_groups
    error('Incompatible nr_test_groups');
end
xl_str = '$n_a$ (log axis)';
font_small = 14;
if strfind(options, 'L')
    fig_rect = [1 1 1024 768];
else
    fig_rect = [1 1 640 480];
end

%% Plot guessing entropy data
h = figure('Position', fig_rect);
for i=1:nr_exp
    semilogx(nr_traces_vec, G(i,:), ...
        'Color', slines{i}.Color, ...
        'LineStyle', slines{i}.LineStyle, ...
        'LineWidth', slines{i}.LineWidth, ...
        'Marker', slines{i}.Marker);        
    hold on;
end
hold off;
if isempty(strfind(options, 'y'))
    ylim(yrange);
end

set(gca,'FontSize', font_size);
title_str = sprintf('Guessing entropy\n%s', ...
                     title_ge);
set(0, 'DefaulttextInterpreter', 'latex');
lh = xlabel(xl_str);
set(lh, 'FontSize', font_size);
if ~isempty(title_ge)
    title(title_str);
end
set(0, 'DefaulttextInterpreter', 'none');
ylabel('Guessing entropy (bits)');
if ~isempty(slegend)
    set(gca,'FontSize', font_small);
    lh = legend(slegend);
    set(lh, 'interpreter', 'none');
end
if ~isempty(strfind(options, 'g'))
    grid on;
end

orient landscape;
print(h, '-dpdf', [rpath, rprefix, 'guess_entropy.pdf']);

end
