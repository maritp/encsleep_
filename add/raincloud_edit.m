% Use like: h = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
% Where 'data' is an M x N cell array of N data series and M measurements
% And 'colours' is an N x 3 array defining the colour to plot each data series
% plot_top_to_bottom: Default plots left-to-right, set to 1 to rotate.
% density_type: 'ks' (default) or 'RASH'. 'ks' uses matlab's inbuilt 'ksdensity' to
% determine the shape of the rainclouds. 'RASH' will use the 'rst_RASH'
% method from Cyril Pernet's Robust Stats toolbox, if that function is on
% your matlab path.
% bandwidth: If density_type == 'ks', determines bandwidth of density estimate
% h is a cell array containing handles of the various figure parts:
% h.p{i,j} is the handle to the density plot from data{i,j}
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.m(i,j) is the handle to the single, large dot that represents nanmean(data{i,j})
% h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)

%% TO-DO:
% Patch can create colour gradients using the 'interp' option to 'FaceColor'. Allow this?

function h = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
%% check dimensions of data

[n_plots_per_series, n_series] = size(data);

% make sure we have correct number of colours
assert(all(size(colours) == [n_series 3]), 'number of colors does not match number of plot series');

%% default arguments

if nargin < 3
    plot_top_to_bottom  = 0;    % left-to-right plotting by default
end

if nargin < 4
    density_type        = 'ks'; % use 'ksdensity' to create cloud shapes
end

if nargin < 5
    bandwidth           = [];   % let the function specify the bandwidth
end

%% Calculate properties of density plots

% Probably okay to hard-code this as it just determines the granularity of
% the density estimate
density_granularity = 200;

n_bins = repmat(density_granularity, n_plots_per_series, n_series);

% calculate kernel densities
for i = 1:n_plots_per_series
    for j = 1:n_series
       
        switch density_type
            
            case 'ks'
                
                % compute density using 'ksdensity'
                [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', bandwidth);
                
            case 'rash'
                
                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found 
                assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
                
                % compute density using RASH
                [x{i, j}, ks{i, j}] = rst_RASH(data{i, j});
                
                % override default 'n_bins' as rst_RASH determines number of bins
                n_bins(i, j) = size(ks{i, j}, 2);
        end
        
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i, j}     = (1:n_bins(i, j) - 1)';
        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)];
        
    end
end

% determine spacing between plots
spacing     = 3 * nanmean(nanmean(cellfun(@max, ks)));
ks_offsets  = [0:n_plots_per_series-1] .* spacing;

% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

% calculate patch vertices from kernel density
for i = 1:n_plots_per_series
    for j = 1:n_series
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
    end
end

%% jitter for the raindrops

jit_width = spacing/20;
raindrop_size = 25;

for i = 1:n_plots_per_series
    for j = 1:n_series
        jit{i,j} = jit_width + rand(1, length(data{i,j})) * jit_width +j*spacing*1/5 -jit_width;
    end
end

%% descriptives of data

cell_nanmeans = cellfun(@nanmean, data);
cell_std = cellfun(@nanstd, data);
cell_n = cellfun(@numel, data);
cell_95cis = 1.96*(cell_std./sqrt(cell_n));

%% plot
% note - we *could* plot everything here in one big loop, but then
% different figure parts would overlay each other in a silly way.

hold on

% patches
for i = 1:n_plots_per_series
    for j = 1:n_series
        
        % plot patches
        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        
        wdth = jit_width;
        % scatter rainclouds
        h.s{i, j} = scatter(data{i, j}, -jit{i, j} + ks_offsets(i) + wdth, 'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5, 'SizeData', raindrop_size);
        
        % 95% ci lines
        h.q = line([cell_nanmeans(i,j) - cell_95cis(i,j)  cell_nanmeans(i,j) + cell_95cis(i,j)],...
            [(ks_offsets(i) + 0.5*wdth -nanmean(jit{i,j}-spacing*1/9)) (ks_offsets(i) + 0.5*wdth -nanmean(jit{i,j}-spacing*1/9))],...
            'col',[0 0 0],'LineWidth',1.5);

        % mean dots
        h.k = scatter(cell_nanmeans(i,j), ks_offsets(i) + 0.5*wdth -nanmean(jit{i,j}-spacing*1/9),...
            'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1, 'SizeData', raindrop_size * 3.5, 'LineWidth', 1.5);
        
    end
end


%% clear up axis labels

% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
ax(1) = gca;
set(ax(1), 'YTick', fliplr(ks_offsets));
set(ax(1), 'YTickLabel', n_plots_per_series:-1:1);

%% set xlabel
set(ax(1), 'XTick', [-2 -1 0 1 2])

%% determine plot rotation
% default option is left-to-right
% plot_top_to_bottom can be set to 1 
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here. That's why the
% logical is constructed the way it is.

% rotate and flip
if ~plot_top_to_bottom
    view([90 -90]);
    axis ij
end

%% add y axis for correlation values
ax(1) = gca;
ax(2) = axes('Position',ax(1).Position,...
    'YLim', ax(1).XLim,...
    'XAxisLocation','top', 'YAxisLocation', 'right', 'color','none') 
%% clear up axis labels

% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
% set(ax(1), 'YTick', fliplr(ks_offsets));
% set(ax(1), 'YTickLabel', n_plots_per_series:-1:1);

%% set xlabel
% set(ax(1), 'XTick', [-2 -1 0 1 2])
set(ax(2), 'YTick', [-2 -1 0 1 2], 'YTickLabel', [-0.96 -0.76 0 0.76 0.96])
set(ax(2), 'XTickLabel', [], 'XTick', []) 



