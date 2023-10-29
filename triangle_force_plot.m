

% h(1) = figure('WindowStyle', 'docked');
% % f = estimated_force_log;
% % [f, ~] = cell_force_field(1, force_number, force_scale);
% % [~, ~, f, ~] = initial_estimation(cell_number, force_scale);
% f = efb;
% ef = efb;
% % tempf = f;
% % f(1:2:end) = tempf(1:2:end);
% % f(2:2:end) = tempf(2:2:end);
% % tempf = ef;
% % ef(1:2:end) = tempf(1:2:end);
% % ef(2:2:end) = tempf(2:2:end);
% % f = tf;
% % [~, ~, ~, TFin] = cell_shape(division_number, 1);
% % f(TFin) = 0;
% [m, ~] = mag_dir_vector(ef, 2);
% % maxmag = max(m);
% % minmag = min(m);
% maxmag = max(ef_img_Lasso, [], 'all');
% minmag = 0;
% 
% % maxmag = 200 / 75000;
% minmag = 0;
% 
% % [M, ~] = mag_dir_map(division_number, f, true, maxmag, minmag);
% M = ef_img_Bayes / maxmag;
% M = M';
% x = linspace(0, width, 2 * division_number + 1);
% x = x(2:2:2 * division_number + 1);
% [X, Y] = meshgrid(x, x);
% x = linspace(0, width, 20 * division_number + 1);
% [Xq, Yq] = meshgrid(x, x);
% Mq = interp2(X, Y, M, Xq, Yq, 'makima');
% Mq = Mq * max(M, [], 'all') / max(Mq, [], 'all');
% Mq(Mq < 0) = 0;
% Mq = sqrt(Mq);
% pq = pcolor(padarray(Xq, [1, 1], -0.1, 'pre'), padarray(Yq, [1, 1], -0.1, 'pre'), padarray(Mq', [1, 1], 1, 'pre'));
% pq.LineStyle = 'none';
% xlim([0 width]);
% ylim([0 height]);
% daspect([1 1 1]);
% xticks([]);
% yticks([]);
% % line([10 - 0.5, 15 - 0.5], [0.5, 0.5], 'Color', 'white', 'LineWidth', 3);
% 
% colormap viridis(100);
% 
% h(2) = figure('WindowStyle', 'docked');
% % c = cell_shape(division_number, cell_number);
% % c = polyshape(c);
% % p = plot(c);
% % p.FaceAlpha = 0;
% 
% hold on;
% 
% % triangle_plot(division_number, f, grid_points(division_number), maxmag, minmag);
% 
% triangle_plot(division_number, f, grid_points(division_number), max(sqrt(f(1:2:end) .^ 2 + f(2:2:end) .^ 2)), 0);
% % 
% % figure;
% % maxmag = max(sqrt(gradx .^ 2 + grady .^ 2));
% % maxmag = max(max(sqrt(gradx .^ 2 + grady .^ 2)));
% % minmag = min(min(sqrt(gradx .^ 2 + grady .^ 2)));
% % gradx = reshape(gradx, 1, []);
% % grady = reshape(grady, 1, []);
% % gradf = reshape([gradx; grady], 1, []);
% % gX = g.xs{1};
% % gX = reshape(gX, 1, []);
% % gY = g.xs{2};
% % gY = reshape(gY, 1, []);
% % ggX = reshape([gX; gY], 1, []);
% % gradf = -(ggX - 0.5);
% % triangle_plot(division_number, gradf, ggX, maxmag, minmag);
% 
% xlim([0 1]);
% ylim([0 1]);
% daspect([1 1 1]);
% 
% xticks([]);
% yticks([]);
% 
% colormap viridis(100);
% % colormap gray(100);
% c = colorbar;
% c.Ticks = c.Limits;
% c.TickLabels = {minmag * 0.4 * 300, maxmag * 0.4 * 300};
% 
% savefig(h, 'triangle_force_plot.fig');
% 
% Numcolor = size(viridis(100), 1);
% Imgmin = min(Mq(:));
% Imgmax = max(Mq(:));
% mappedImage = uint8( (Mq-Imgmin)./(Imgmax-Imgmin).* (Numcolor-1) );
% imwrite(mappedImage, viridis(100), 'efb.tif');
% 
% function triangle_plot(division_number, F, gX, maxmag, minmag)
% 
% X = gX(1:2:end);
% Y = gX(2:2:end);
% 
% [M, D] = mag_dir_vector(F, 2);
% M = (M - minmag) ./ (maxmag - minmag);
% M(M > 1) = 1;
% M(M < 0) = 0;
% M = ceil(M * 100);
% M(M <= 0) = 1;
% 
% DX = D(1:2:end);
% DY = D(2:2:end);
% 
% x = zeros(1, length(F) / 2);
% y = zeros(size(x));
% u = zeros(size(x));
% v = zeros(size(x));
% 
% hold on;
% for i=1:length(F) / 2
%     
%     if M(i) < 50
%         continue;
%     end
%     
% %     plot([0 0.9 * DX(i) / division_number] + X(i), [0 0.9 * DY(i) / division_number] + Y(i), ...
% %         'Marker', 'none', 'Color', 'k', ...
% %         'MarkerIndices', 1, 'MarkerEdgeColor', 'none', ...
% %         'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
% %         'LineWidth', 2);
%     
%     x(i) = X(i);
%     y(i) = Y(i);
%     u(i) = 0.9 * DX(i) / division_number;
%     v(i) = 0.9 * DY(i) / division_number;
% end
% 
% q = quiver(x, y, u, v, 'AutoScale', 'off', 'LineWidth', 1, 'ShowArrowHead', 'off');
% q.MaxHeadSize = 0.5;
% end
% 
% function R = rotate_matrix(theta)
%     R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% end