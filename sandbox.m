%%
max_f = 0.0603;
min_f = 0;
ef_img = tf_img;
ef_img = (ef_img - min_f) ./ (max_f - min_f);
Xq = linspace(0, width, 3 * division_number);
Yq = linspace(0, height, 3 * division_number);
[Xq, Yq] = meshgrid(Xq, Yq);
q = interp2(X, Y, ef_img, Xq, Yq);
q(:, end) = 1;
figure;
p = pcolor(Xq, Yq, q);
p.LineStyle = 'none';
colormap viridis;
colorbar;
xticks([]);
yticks([]);
title('Synthetic force');
pbaspect([1 1 1]);
line([10 - 0.5, 15 - 0.5], [0.5, 0.5], 'Color', 'white', 'LineWidth', 3);

%%

load('BayesROC.mat');
tpfb = tpfs;
fpfb = fpfs;
clear tpfs fpfs;
load('RidgeROC.mat');
tpfr = tpfs;
fpfr = fpfs;
clear tpfs fpfs;
load('LassoROC.mat');
tpfl = tpfs;
fpfl = fpfs;
clear tpfs fpfs;
load('ENROC.mat');
tpfe = tpfs;
fpfe = fpfs;
clear tpfs fpfs;

figure; hold on;
ax = gca;
ax.FontSize = 16;
plot(fpfb, tpfb, 'Marker', 'o', 'MarkerFaceColor', 'white', 'DisplayName', 'Bayes');
plot(fpfr, tpfr, 'Marker', 'o', 'MarkerFaceColor', 'white', 'DisplayName', 'Ridge');
plot(fpfl, tpfl, 'Marker', 'o', 'MarkerFaceColor', 'white', 'DisplayName', 'Lasso');
plot(fpfe, tpfe, 'Marker', 'o', 'MarkerFaceColor', 'white', 'DisplayName', 'EN');
title('ROC');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside', 'AutoUpdate', 'off');
line([0, 1], [0, 1], 'Color', 'k');
pbaspect([1 1 1]);

%%
clear;
rng default % For reproducibility
n = 50;
X = randn(n, 5);
weights = [0;2;0;-3;0]; % Only two nonzero coefficients
y = X*weights + randn(n, 1)*0.1; % Small added noise

[B, FitInfo] = lasso(X,y, 'CV', 10);
lassoPlot(ef,FitInfo,'PlotType','CV');
legend('show') % Show legend

%%
clear;

n = 30;

cell_id = 5;

[CSPs, II, EI, TFin] = cell_shape(n, cell_id);
II = union(II, EI);

shape = false(n);
shape(II) = true;

II = repmat(reshape(shape, 1, []), [n ^ 2, 1]); % inner index
m = length(find(shape));

width = 15;
height = width;

cs = polyshape(CSPs .* width);

x = linspace(0, width, n);
y = linspace(0, height, n);

[X, Y] = meshgrid(x, y);
X2 = repmat(reshape(X, 1, []), [n ^ 2, 1]);
Y2 = repmat(reshape(Y, 1, []), [n ^ 2, 1]);
DX = X2 - repmat(X2(1, :)', [1, n ^ 2]);
DX(DX < 0) = -1;
DX(DX > 0) = 1;
DY = Y2 - repmat(Y2(1, :)', [1, n ^ 2]);
DY(DY < 0) = -1;
DY(DY > 0) = 1;
D = sqrt(DX .^ 2 + DY .^ 2);
D(~II) = 0;
U = reshape(sum(D .* DX, 2), [n, n]);
U(~II(1, :)) = 0;
U = U ./ max(U, [], 'all') ./ sqrt(2);
V = reshape(sum(D .* DY, 2), [n, n]);
V(~II(1, :)) = 0;
V = V ./ max(V, [], 'all') ./ sqrt(2);

figure; hold on;
plot(cs, 'DisplayName', 'Cell');
scatter(reshape(X, 1, []), reshape(Y, 1, []), 'Marker', 'o', 'MarkerFaceColor', 'b', 'SizeData', 3, 'DisplayName', 'Grid');
quiver(X, Y, U, V, 'Color', 'k', 'DisplayName', 'Traction');
xticks([]);
yticks([]);
legend('Location', 'eastoutside');
daspect([1 1 1]);

load("cell_shape\force_field_" + cell_id + ".mat")
U = reshape(U, 1, []);
V = reshape(V, 1, []);
f = reshape([U; V], [], 1)';

%%
save("cell_shape\force_field_" + cell_id + ".mat", 'EI', 'II', 'cs', 'd', 'f');

%% 
clear;
load('cell_shape\force_field_1.mat')
B = f;
max(sqrt(B(1:2:length(B)) .^ 2 + B(2:2:length(B))), [], 'all')
load('cell_shape\force_field_1 - ƒRƒs[.mat')
B = f;
max(sqrt(B(1:2:length(B)) .^ 2 + B(2:2:length(B))), [], 'all')

%%
clear;
load('201022I.mat')
load("201018beadsdata.mat");

%%
GSarea = beadsdata.GSarea(end:-1:1, :);
c = [300, 150];
yrange = c(2) + (-100:100);
ybound = [yrange(1), yrange(end)];
xrange = c(1) + (-100:100);
xbound = [xrange(1), xrange(end)];

GSarea = GSarea(yrange, xrange);

% figure; hold on;
% [width, height] = size(GSarea);
% plot([0, width], [0, height]);
% imagesc(GSarea); axis equal; axis tight;

division_number = 30;
n = division_number;
m = division_number;
[X, Y, U, V] = prior_force_from_GSarea(GSarea, 30);

figure; hold on;
imagesc(GSarea); axis equal; axis tight;
scatter(reshape(X, 1, []), reshape(Y, 1, []), 'Marker', 'o', 'MarkerFaceColor', 'b', 'SizeData', 3, 'DisplayName', 'Grid');
quiver(X, Y, U, V, 'Color', 'k', 'DisplayName', 'Traction');
xticks([]);
yticks([]);
legend('Location', 'eastoutside');
daspect([1 1 1]);

%%

frame = 2;
coords1 = beadsdata.coords(:, :, 1);
coords2 = beadsdata.coords(:, :, frame);
IBLx = coords1(:, 1);
IBLy = coords1(:, 2);
IBLy = IBLy(:);
IBLx = IBLx(:);
I = (IBLx > xbound(1)) & (xbound(2) > IBLx) & (IBLy > ybound(1)) & (ybound(2) > IBLy);
IBLx = IBLx(I) - xbound(1);
IBLy = IBLy(I) - ybound(1);
BD = coords2 - coords1;
BDx = BD(:, 1);
BDy = BD(:, 2);
BDx = BDx(I);
BDy = BDy(I);
Bx = U(:) * 1e-1;
By = V(:) * 1e-1;
B = [Bx, By]';
B = B(:);

% [IBLx, IBLy, BDx, BDy] = noise_cut(IBLx, IBLy, BDx, BDy);
BD = [BDx, BDy]';
BD = BD(:);

tfmc = TFM_computation(X, Y, Bx, By, IBLx, IBLy);
tfmc = tfmc.simulate();
G = tfmc.G;

bayestfm = BayesTFM2D(1, 1e-1, 1);

II = find(I);
III = reshape([2 * II - 1; 2 * II], 1, []);

tic;
bayestfm = bayestfm.fit(BD, tfmc.G(:, III), B(III));
ef = zeros(1, 2 * n * m);
ef(III) = bayestfm.force;
% disp(mean((tf - ef) .^ 2, 'all') * 2);
toc;

efx = ef(1:2:length(ef));
efy = ef(2:2:length(ef));

% figure; hold on;
% quiver(IBLx, IBLy, BD(1:2:length(BD)), BD(2:2:length(BD)));
% quiver(IBLx, IBLy, tfmc.BDx, tfmc.BDy);

[width, height] = size(GSarea);
figure; hold on;
imagesc(GSarea); colormap gray; axis equal; axis tight;
scatter(reshape(X, 1, []), reshape(Y, 1, []), 'Marker', 'o', 'MarkerFaceColor', 'b', 'SizeData', 3, 'DisplayName', 'Grid');
hold on; quiver(X, Y, reshape(B(1:2:length(B)), n, m), reshape(B(2:2:length(B)), n, m), 'DisplayName', 'prior');
quiver(IBLx, IBLy, BD(1:2:length(BD)), BD(2:2:length(BD)), 'DisplayName', 'Bead displacements');
quiver(X, Y, reshape(efx, n, m), reshape(efy, n, m), 'DisplayName', 'estimated', 'LineWidth', 2);
title("Bayes");
xlim([0 width]);
ylim([0 height]);
daspect([1 1 1]);
legend();


function [X, Y, U, V] = prior_force_from_GSarea(GSarea, division_number)
[shape, X, Y] = shape_from_GSarea(GSarea, division_number);
n = size(shape, 1);

II = repmat(reshape(shape, 1, []), [n ^ 2, 1]); % inner index

X2 = repmat(reshape(X, 1, []), [n ^ 2, 1]);
Y2 = repmat(reshape(Y, 1, []), [n ^ 2, 1]);
DX = X2 - repmat(X2(1, :)', [1, n ^ 2]);
DX(DX < 0) = -1;
DX(DX > 0) = 1;
DY = Y2 - repmat(Y2(1, :)', [1, n ^ 2]);
DY(DY < 0) = -1;
DY(DY > 0) = 1;
D = sqrt(DX .^ 2 + DY .^ 2);
D(~II) = 0;
U = reshape(sum(D .* DX, 2), [n, n]);
U(~II(1, :)) = 0;
U = U ./ max(U, [], 'all') ./ sqrt(2);
V = reshape(sum(D .* DY, 2), [n, n]);
V(~II(1, :)) = 0;
V = V ./ max(V, [], 'all') ./ sqrt(2);
end

function [shape, X, Y] = shape_from_GSarea(GSarea, division_number)

[width, height] = size(GSarea);
X = linspace(1, width, division_number);
Y = linspace(1, height, division_number);
[X, Y] = meshgrid(X, Y);

roundX = round(X);
roundY = round(Y);
shape = false(size(roundX));

for i=1:size(shape, 1)
    for j=1:size(shape, 2)
        shape(j, i) = (GSarea(roundX(i, j), roundY(i, j)) == 1);
    end
end

end

function [cIBLx, cIBLy, cBDx, cBDy] = noise_cut(IBLx, IBLy, BDx, BDy)
BD = sqrt((BDx .^ 2 + BDy .^ 2));
[counts, x] = hist(BD, 10);

fun = @(coef, xdata) coef(1) * exp(coef(2) * xdata);
coef0 = [1, -1];
coef = lsqcurvefit(fun, coef0, x, counts);

Thr = - log(2) / coef(2);
TI = BD > Thr;

cIBLx = IBLx(TI);
cIBLy = IBLy(TI);
cBDx = BDx(TI);
cBDy = BDy(TI);
end