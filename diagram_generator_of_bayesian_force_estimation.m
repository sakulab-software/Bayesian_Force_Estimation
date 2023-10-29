division_number = 1;
n = division_number;
m = division_number;
N = n * m;

B = 100;

width = 15;
height = width;

force_scale = 1e-1;

X = width / 2;
Y = height / 2;

Fx = zeros(n, m);
Fy = zeros(n, m);

rng(49);
for i=1
    Fx(i) = randn;
    Fy(i) = randn;
end

Fx = Fx * force_scale;
Fy = Fy * force_scale;

rng(60);
IBLx = rand(B, 1) * width;
IBLy = rand(B, 1) * height;

tfmc = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy);
tfmc = tfmc.simulate();
[BDx, BDy] = tfmc.observe_without_noise();
BD = sqrt(BDx .^ 2 + BDy .^ 2);
% [dBDx, dBDy] = tfmc.observe(mean(BD) * 10 ^ (- (25) / 10), 49);
[dBDx, dBDy] = tfmc.observe(1e-10, 49);
% [dBDx, dBDy] = tfmc.observe_without_noise();

figure, quiver(X, Y, Fx, Fy, 'DisplayName', 'force');
hold on; quiver(IBLx, IBLy, BDx, BDy, 'DisplayName', 'disp');
hold on; quiver(IBLx, IBLy, dBDx, dBDy, 'DisplayName', 'disp(noise)');
title("Bead");
xlim([0 width]);
ylim([0 height]);
daspect([1 1 1]);
legend();

%%
BD = reshape([dBDx, dBDy]', 1, [])';
[FxData, FyData] = meshgrid(-0.6:0.01:0.6, -0.6:0.01:0.6);
Li = zeros(size(FxData));
for i=1:numel(Li)
    Li(i) = exp(-1e-1 * norm(BD - tfmc.G * [FxData(i); FyData(i)]) ^ 2);
end

% figure;
% pc = pcolor(FxData, FyData, Li);
% pc.LineStyle = 'none';
% colormap gray;
% colorbar;
% daspect([1 1 1]);
% hold on;
% line([-1, 1], [0, 0], 'Color', 'white', 'LineWidth', 2);
% line([0, 0], [-1, 1], 'Color', 'white', 'LineWidth', 2);

Pr = zeros(size(FxData));
for i=1:numel(Pr)
    Pr(i) = exp(-norm([2 * Fx; 2 * Fy] - [FxData(i); FyData(i)]));
end
% 
% figure;
% pc = pcolor(FxData, FyData, Pr);
% pc.LineStyle = 'none';
% colormap gray;
% colorbar;
% daspect([1 1 1]);
% hold on;
% line([-1, 1], [0, 0], 'Color', 'white', 'LineWidth', 2);
% line([0, 0], [-1, 1], 'Color', 'white', 'LineWidth', 2);
% 
% figure;
% pc = pcolor(FxData, FyData, Pr * Li);
% pc.LineStyle = 'none';
% colormap gray;
% colorbar;
% daspect([1 1 1]);
% hold on;
% line([-1, 1], [0, 0], 'Color', 'white', 'LineWidth', 2);
% line([0, 0], [-1, 1], 'Color', 'white', 'LineWidth', 2);

shape = zeros(2, 100);
for i=1:size(shape, 2)
    theta = 2 * pi / size(shape, 2) * i;
    alpha = deg2rad(45);
    R = [cos(alpha), -sin(alpha);
         sin(alpha), cos(alpha)];
    shape(:, i) = ...
        R * ([cos(theta); sin(theta)] .* ...
        [1; 0.5] * ...
        0.7) + ...
        [-0.45; -0.53];
end

IBLx2 = (IBLx - width / 2) / 3;
IBLy2 = (IBLy - height / 2) / 3;

figure; hold on;
ax = gca;
ax.FontSize = 15;
pc = pcolor(FxData, FyData, Pr);
[c, h] = contour(FxData, FyData, Li, 'LineColor', 'k', 'ShowText', 'off');
clabel(c, h);
plot(shape(1, :), shape(2, :), 'k--', 'LineWidth', 2);
pc.LineStyle = 'none';
colormap gray;
c = colorbar;
c.Label.String = "Probability density of Prior";
c.Label.FontSize = 18;
daspect([1 1 1]);
% quiver(IBLx2, IBLy2, dBDx, dBDy, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1 / 100);
line([-1, 1], [0, 0], 'Color', 'white', 'LineWidth', 2);
line([0, 0], [-1, 1], 'Color', 'white', 'LineWidth', 2);
quiver(0, 0, -0.32, -0.1, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 3);
xlabel('X component [Pa]');
xlim([-0.6, 0.6]);
xticks(0);
ylabel('Y component [Pa]');
ylim([-0.6, 0.6]);
yticks(0);

