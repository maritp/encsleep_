function [x, y] = getting_layout(layout_name, neighbours_name, plot)

% getting_layout('Bham-64CH-Lay.mat', 'Bham-64CH-Neighbours.mat', 0)
% plot = 1 if yes and 0 if no


cfg = [];
cfg.layout = layout_name;
x = ft_prepare_layout(cfg);


cfg = [];
cfg.method = 'template';
cfg.template = neighbours_name;
y = ft_prepare_neighbours(cfg);


if plot
    figure;
    ft_plot_lay(x);
    
    cfg = [];
    cfg.layout = x;
    cfg.neighbours = y;
    ft_neighbourplot(cfg);
end

