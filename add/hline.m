function hline(y_pos,spec)
if nargin < 2
    spec = 'k-';
end
range=axis;
hold on
plot([range(1) range(2)],[y_pos y_pos],spec)
hold off
