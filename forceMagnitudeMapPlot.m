function forceMagnitudeMapPlot(X, Y, M, maxM)

[m, n] = size(X);
minX = min(X, [], 'all');
maxX = max(X, [], 'all');
minY = min(Y, [], 'all');
maxY = max(Y, [], 'all');
x = linspace(minX, maxX, 10 * n);
y = linspace(minY, maxY, 10 * m);
[Xq, Yq] = meshgrid(x, y);
Mq = interp2(X, Y, M, Xq, Yq, 'linear');
% Mq = Mq * max(M, [], 'all') / max(Mq, [], 'all');
Mq(Mq < 0) = 0;

if ~exist("maxM", 'var')
    maxM = max(M, [], 'all');
end

pq = pcolor(padarray(Xq, [1, 1], -0.1, 'pre'), padarray(Yq, [1, 1], -0.1, 'pre'), padarray(Mq, [1, 1], maxM, 'pre'));
pq.LineStyle = 'none';
xlim([minX maxX]);
ylim([minY maxY]);
daspect([1 1 1]);
xticks([]);
yticks([]);

% colormap(viridis(100));
colorbar;

end