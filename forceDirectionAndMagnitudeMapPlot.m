function forceDirectionAndMagnitudeMapPlot(X, Y, Fx, Fy, M, ax1, ax2, maxM)

[m, n] = size(X);
minX = min(X, [], 'all');
maxX = max(X, [], 'all');
minY = min(Y, [], 'all');
maxY = max(Y, [], 'all');
x = linspace(minX, maxX, 10 * n);
y = linspace(minY, maxY, 10 * m);
[Xq, Yq] = meshgrid(x, y);
Mq = interp2(X, Y, M, Xq, Yq, 'makima');
Mq = Mq * max(M, [], 'all') / max(Mq, [], 'all');
Mq(Mq < 0) = 0;

if exist("maxM", 'var')
    Mq = Mq / maxM;
else
    maxM = max(Mq, [], 'all');
end

if ~exist('ax1', 'var')
    figure('WindowStyle', 'docked');
    ax1 = gca;
end

if ~exist("ax2", 'var')
    figure('WindowStyle', 'docked');
    ax2 = gca;
end

pq = pcolor(ax1, padarray(Xq, [1, 1], -0.1, 'pre'), padarray(Yq, [1, 1], -0.1, 'pre'), padarray(Mq, [1, 1], maxM, 'pre'));
pq.LineStyle = 'none';
xlim(ax1, [minX maxX]);
ylim(ax1, [minY maxY]);
daspect(ax1, [1 1 1]);
xticks(ax1, []);
yticks(ax1, []);

colormap(ax1, viridis(100));
colorbar(ax1);

forceDirectionPlot(ax2, X, Y, Fx, Fy, M, 10);
xlim(ax2, [minX maxX]);
ylim(ax2, [minY maxY]);
daspect(ax2, [1 1 1]);
xticks(ax2, []);
yticks(ax2, []);

end

function forceDirectionPlot(ax, X, Y, Fx, Fy, M, p)

Dx = Fx ./ M;
Dy = Fy ./ M;
ind = M < max(M, [], 'all') * p / 100;

quiver(ax, X(~ind), Y(~ind), Dx(~ind), Dy(~ind), 'AutoScale', 'off', 'LineWidth', 1, 'ShowArrowHead', 'off');

end
