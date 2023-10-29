function ind = forceDirectionPlot(X, Y, Fx, Fy, M, p, ind)

minX = min(X, [], 'all');
maxX = max(X, [], 'all');
minY = min(Y, [], 'all');
maxY = max(Y, [], 'all');

Dx = Fx ./ M / 2.2;
Dy = Fy ./ M / 2.2;

if ~exist('ind', 'var')
    ind = M > max(M, [], 'all') * p / 100;
end

quiver(X(ind) - Dx(ind) / 2, Y(ind) - Dy(ind) / 2, Dx(ind), Dy(ind), 'AutoScale', 'off', 'LineWidth', 1, 'ShowArrowHead', 'off');
xlim([minX maxX]);
ylim([minY maxY]);
daspect([1 1 1]);
xticks([]);
yticks([]);
end