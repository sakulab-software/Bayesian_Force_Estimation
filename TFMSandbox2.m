%
% Simulation and BATELAB Analysis as of Dec. 28, 2020
%

%% 

visible = true;

divisionNumber = 30;
BeadNumber = 50;
w = 15;
h = w;
forceScale = 1e-1;
cellID = 5;
forceID = 1;
randomState = 51;

[X, Y, SFx, SFy, IBLx, IBLy, ~, BDx, BDy, dBD, dBDx, dBDy, tfmc] = ... 
    generateSyntheticData(divisionNumber, BeadNumber, w, h, forceScale, cellID, forceID, randomState);

if visible
    figure, quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(IBLx, IBLy, BDx, BDy, 'DisplayName', 'disp');
    hold on; quiver(IBLx, IBLy, dBDx, dBDy, 'DisplayName', 'disp(noise)');
    title("Bead");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end

sf_img = sqrt(SFx .^ 2 + SFy .^ 2);
% 1 x (2 * divisionNumber ^ 2)
sf = reshape([reshape(SFx, 1, []); reshape(SFy, 1, [])], 1, []);

%% 

visible = true;

if ~visible
    th = beadDispThreshold(dBD, dBDx, dBDy);
else
    [th, dBDx2, dBDy2] = beadDispThreshold(dBD, dBDx, dBDy);
    figure;
    hold on; quiver(IBLx, IBLy, BDx, BDy, 'DisplayName', 'disp');
    hold on; quiver(IBLx, IBLy, dBDx2, dBDy2, 'DisplayName', 'disp(noise)');
    title("Bead");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end

%%

visible = true;

gm = beadDensityMap(X, Y, BeadNumber, IBLx, IBLy, divisionNumber, cellID);

if visible
    figure; hold on;
    p = pcolor(X, Y, gm);
    p.LineStyle = 'none';
    colorbar;
    colormap viridis;
    cell = cell_shape(divisionNumber, cellID);
    cell = polyshape({cell(:, 1) * w}, {cell(:, 2) * h});
    plot(cell);
    title("Local bead density");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
end

%% 

%
% Previouse Ridge Code
%
% tic;
% fold = 10;
% kfolds = kfoldsplitter(tfmc.G, reshape([dBDx, dBDy]', 1, [])', B, fold);
% bestparam = RidgeCV(kfolds, logspace(-10, 3, 14));
% efr = ridge(reshape([dBDx, dBDy]', 1, [])', tfmc.G, bestparam, 0);
% efr = efr(2:end);% + ef(1);
% efr = efr';
% disp(mean((tf - efr) .^ 2, 'all') * 2);
% toc;
% 
% efx = efr(1:2:end);
% efy = efr(2:2:end);
% 
% figure, quiver(X, Y, Fx, Fy, 'DisplayName', 'force');
% hold on; quiver(X, Y, reshape(efx, n, m), reshape(efy, n, m), 'DisplayName', 'estimated');
% title("Ridge");
% xlim([0 width]);
% ylim([0 height]);
% daspect([1 1 1]);
% legend();

visible = true;

if visible
    tic;
end

[REFx, REFy, REFMap] = TFMWithRidge(divisionNumber, tfmc.G, dBDx, dBDy);

if visible
    disp(mean((SFx - REFx) .^ 2 + (SFy - REFy) .^ 2, 'all') * 2);
    toc;

    figure, quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, REFx, REFy, 'DisplayName', 'estimated');
    title("Ridge");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end


%%

visible = true;

if visible
    tic;
end

[LEFx, LEFy, LEFMap] = TFMWithLASSO(divisionNumber, tfmc.G, dBDx, dBDy);

if visible
    disp(mean((SFx - LEFx) .^ 2 + (SFy - LEFy) .^ 2, 'all') * 2);
    toc;

    figure, quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, LEFx, LEFy, 'DisplayName', 'estimated');
    title("LASSO");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end

%%

visible = true;

if visible
    tic;
end

[EEFx, EEFy, EEFMap] = TFMWithElasticNet(divisionNumber, tfmc.G, dBDx, dBDy, 1e-3);

if visible
    disp(mean((SFx - EEFx) .^ 2 + (SFy - EEFy) .^ 2, 'all') * 2);
    toc;

    figure, quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, EEFx, EEFy, 'DisplayName', 'estimated');
    title("Elastic Net");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end

%%

visible = true;

if visible
    tic;
end

[BEFx, BEFy, BEFMap] = TFMWithBayesianEstimation(divisionNumber, tfmc.G, dBDx, dBDy, cellID, forceScale, 1, 1e0, gm);

if visible
    disp(mean((SFx - BEFx) .^ 2 + (SFy - BEFy) .^ 2, 'all') * 2);
    toc;

    [~, ~, Pr, ~] = initial_estimation(cellID, forceScale);
    Pr = Pr';
    Prx = reshape(Pr(1:2:end), divisionNumber, divisionNumber);
    Pry = reshape(Pr(2:2:end), divisionNumber, divisionNumber);
    
    figure, quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, BEFx, BEFy, 'DisplayName', 'estimated');
%     hold on; quiver(X, Y, Prx, Pry, 'DisplayName', 'prior');
    title("Bayes");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
end

rng('default');

%%
load('BATELABResult.mat');
fdata = data;
load('BATLABResult.mat');
bdata = data;
load('roi.mat');

forceScale = 1e-1;
width = 15;
height = 15;

n = 30;
m = n;

X = reshape(fdata.x, 30, 30) * width;
Y = reshape(fdata.y, 30, 30) * height;

II = find(II);
II = reshape([2 * II - 1, 2 * II]', 1, []);

IBLx = (bdata.x - roiPosition(1)) / roiPosition(3) * width;
IBLy = (bdata.y - roiPosition(2)) / roiPosition(4) * height;
BD = reshape([bdata.u / roiPosition(3) * width, bdata.v / roiPosition(4) * height]', 1, [])';

X = linspace(0, width, n);
Y = linspace(0, height, m);
[X, Y] = meshgrid(X, Y);
gm = local_bead_density(length(IBLx), IBLx, IBLy, 1, n, X, Y);
gmmax_inside = max(gm, [], 'all');
gmmin = min(reshape(gm, 1, []));
gm = (gm - gmmin) ./ (gmmax_inside - gmmin);
gm = 1 ./ (1 + exp(- 10 * (gm - 0.5)));

Fx = reshape(fdata.u, 30, 30) * forceScale .* gm;
Fy = reshape(fdata.v, 30, 30) * forceScale .* gm;
PFx = reshape(fdata.u, 30, 30) * forceScale;
PFy = reshape(fdata.v, 30, 30) * forceScale;

tfmc = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy);
tfmc = tfmc.simulate();

BD2 = reshape([tfmc.BDx, tfmc.BDy]', 1, [])';

B = zeros(2 * n * m, 1);
B(1:2:end) = reshape(Fx, 1, []);
B(2:2:end) = reshape(Fy, 1, []);
B2 = zeros(2 * n * m, 1);
B2(1:2:end) = reshape(PFx, 1, []);
B2(2:2:end) = reshape(PFy, 1, []);

bayestfm = BayesTFM2D(1, 1e0, 10000);

tic;
bayestfm = bayestfm.fit(BD, tfmc.G(:, II), B(II));
ef = zeros(1, 2 * n * m);
ef(II) = bayestfm.force;
toc;

% ef = B2;
% efy = B2(1:2:end);
% efx = B2(2:2:end);
efx = ef(1:2:end);
efy = ef(2:2:end);

figure;
% hold on; quiver(X, Y, reshape(efx, n, m), reshape(efy, n, m), 'DisplayName', 'Estimated Force');
hold on; scatter(IBLx, IBLy);
hold on; quiver(IBLx, IBLy, BD(1:2:end), BD(2:2:end), 'DisplayName', 'Bead Displacements', 'AutoScale', 'off', 'ShowArrowHead', 'off', 'LineWidth', 2);
title("");
xlim([0 width]);
xticks([]);
ylim([0 height]);
yticks([]);
daspect([1 1 1]);

figure;
hold on; quiver(X, Y, PFx, PFy, 'DisplayName', 'Prior Force');
title("Bayes");
xlim([0 width]);
ylim([0 height]);
daspect([1 1 1]);
legend();

ef_img_Bayes = reshape(efx .^ 2 + efy .^ 2, n, m);

%% 
% figure;
% magef = sqrt(efx .^ 2 + efy .^ 2);
% magmax = max(magef);
% magmin = min(magef);
% [counts, x] = hist(magef, 10);
% stem(x, counts);
% T = otsuthresh(counts);
% T = T * (magmax - magmin) + magmin;
% l = line(ones(1, 2) * T, [0, max(counts)]);

%% 

% efx2 = efx;
% efx2(magef < T) = 0;
% efy2 = efy;
% efy2(magef < T) = 0;
% 
% figure, quiver(X, Y, Fx, Fy, 'DisplayName', 'force');
% hold on; quiver(X, Y, reshape(efx2, n, m), reshape(efy2, n, m), 'DisplayName', 'estimated');
% title("Bayes");
% xlim([0 width]);
% ylim([0 height]);
% daspect([1 1 1]);
% legend();

%%
function [X, Y, Fx, Fy, IBLx, IBLy, BD, BDx, BDy, dBD, dBDx, dBDy, Sim] = ... 
    generateSyntheticData(divisionNumber, BeadNumber, w, h, forceScale, cellID, forceID, randomState)

n = divisionNumber;
m = divisionNumber;
N = n * m;

X = linspace(0, w, n);
Y = linspace(0, h, m);
[X, Y] = meshgrid(X, Y);

[tf, ~] = cell_force_field(cellID, forceID, forceScale);
Fx = tf(1:2:length(tf))';
Fy = tf(2:2:length(tf))';
Fx = reshape(Fx, n, m);
Fy = reshape(Fy, n, m);

rng(randomState);
IBLx = rand(BeadNumber, 1) * w;
IBLy = rand(BeadNumber, 1) * h;

Sim = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy);
Sim = Sim.simulate();
[BDx, BDy] = Sim.observe_without_noise();
BD = sqrt(BDx .^ 2 + BDy .^ 2);
[dBDx, dBDy] = Sim.observe(1e-3, 49);
dBD = sqrt(dBDx .^ 2 + dBDy .^ 2);

end

%% 

function [th, thBDx, thBDy] = beadDispThreshold(BD, BDx, BDy, t, visible)

if ~exist("t", 'var')
    t = 2;
end

if ~exist("visible", 'var')
    visible = false;
end

[counts, edges] = histcounts(BD, 10);
x = mean([edges(1:end-1); edges(2:end)]);

fun = @(coef, xdata) coef(1) * exp(coef(2) * xdata);
coef0 = [1, -1];
coef = lsqcurvefit(fun, coef0, x, counts);

th = - log(t) / coef(2);
ind = BD < th;

thBDx = BDx;
thBDy = BDy;
thBDx(ind) = 1e-12;
thBDy(ind) = 1e-12;

if visible
    figure; hold on;
    stem(x, counts);
    line(ones(1, 2) * th, [0, max(counts)]);
    plot(x, fun(coef, x));
end

end

%%
function gm = beadDensityMap(X, Y, BeadNumber, IBLx, IBLy, divisionNumber, cellID)

gm = local_bead_density(BeadNumber, IBLx, IBLy, 1, divisionNumber, X, Y);
[~, II, ~, ~] = cell_shape(divisionNumber, cellID);
gmmax_inside = max(gm(II), [], 'all');
gmmin = min(reshape(gm, 1, []));
gm = (gm - gmmin) ./ (gmmax_inside - gmmin);
% sigmoid transformation to saturate gm over 0.5 
gm = 1 ./ (1 + exp(- 10 * (gm - 0.5)));

end

%%
function [EFx, EFy, EFMap] = TFMWithElasticNet(divisionNumber, G, BDx, BDy, a)

[ef, FitInfo] = lasso(G, reshape([BDx, BDy]', 1, [])', 'CV', 10, 'Alpha', a);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
ef = ef(:, idxLambdaMinMSE)';
EFx = ef(1:2:end);
EFy = ef(2:2:end);
EFx = reshape(EFx, divisionNumber, divisionNumber);
EFy = reshape(EFy, divisionNumber, divisionNumber);
EFMap = sqrt(EFx .^ 2 + EFy .^ 2);

end

%%
function [EFx, EFy, EFMap] = TFMWithRidge(divisionNumber, G, BDx, BDy)
[EFx, EFy, EFMap] = TFMWithElasticNet(divisionNumber, G, BDx, BDy, 1e-4);
end

%%
function [EFx, EFy, EFMap] = TFMWithLASSO(divisionNumber, G, BDx, BDy)
[EFx, EFy, EFMap] = TFMWithElasticNet(divisionNumber, G, BDx, BDy, 1);
end

%%
function [EFx, EFy, EFMap] = TFMWithBayesianEstimation(divisionNumber, G, BDx, BDy, cellID, forceScale, A, B, beadDensityMap)

BD = reshape([BDx, BDy]', 1, [])';
[~, ~, Prior, II] = initial_estimation(cellID, forceScale);
Prior = Prior';

if exist("beadDensityMap", 'var')
    Prior(1:2:end) = Prior(1:2:end) .* reshape(beadDensityMap, [], 1);
    Prior(2:2:end) = Prior(2:2:end) .* reshape(beadDensityMap, [], 1);
end

III = reshape([2 * II - 1; 2 * II], 1, []);

bayestfm = BayesTFM2D(A, B, 10000);
bayestfm = bayestfm.fit(BD, G(:, III), Prior(III));
EF = zeros(1, 2 * divisionNumber ^ 2);
EF(III) = bayestfm.force;

EFx = reshape(EF(1:2:end), divisionNumber, divisionNumber);
EFy = reshape(EF(2:2:end), divisionNumber, divisionNumber);

EFMap = sqrt(EFx .^ 2 + EFy .^ 2);
end

%%

% function kfolds = kfoldsplitter(G, d, b, k)
% p = randperm(b);
% kfolds = cell(1, k);
% for i=1:k
%     s = (b / k) * (i - 1) + 1;
%     e = (b / k) * i;
%     l = false(1, b);
%     l(s:e) = true;
%     testp = p(l);
%     testp = reshape([2 * testp - 1; 2 * testp], 1, [])';
%     trainp = p(~l);
%     trainp = reshape([2 * trainp - 1; 2 * trainp], 1, [])';
%     kfolds{i} = struct;
%     kfolds{i}.Xtrain = G(trainp);
%     kfolds{i}.Xtest = G(testp);
%     kfolds{i}.ytrain = d(trainp);
%     kfolds{i}.ytest = d(testp);
% end
% end
% 
% function bestparam = RidgeCV(kfolds, paramlist)
% n = length(paramlist);
% fold = length(kfolds);
% cvMSE = zeros(fold, n);
% for i=1:n
%     for k=1:fold
%         ef = ridge(kfolds{k}.ytrain, kfolds{k}.Xtrain, paramlist(i), 0);
%         ef = ef(2:end);% + ef(1);
%         disperr = kfolds{k}.ytest - kfolds{k}.Xtest * ef;
%         cvMSE(k, i) = mean(disperr .^ 2, 'all') * 2;
%     end
% end
% figure;
% plot(1:n, mean(cvMSE));
% [~, I] = min(mean(cvMSE));
% bestparam = paramlist(I);
% end
