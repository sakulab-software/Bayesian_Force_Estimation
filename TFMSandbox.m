%% 

visible = true;

divisionNumber = 30;
BeadNumber = 100;
w = 15;
h = w;
forceScale = 5e-1;
cellID = 2;
forceID = 10;
randomState = 50; % For reproductivity

[X, Y, SFx, SFy, IBLx, IBLy, ~, BDx, BDy, dBD, dBDx, dBDy, tfmc] = ... 
    generateSyntheticData(divisionNumber, BeadNumber, w, h, forceScale, cellID, forceID, 1e-2, randomState);

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

U = reshape([dBDx, dBDy]', 1, []);
SFMap = sqrt(SFx .^ 2 + SFy .^ 2);

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

%% Ridge

visible = true;

if visible
    tic;
end

[REFx, REFy, REFMap] = TFMWithRidge(divisionNumber, tfmc.G, dBDx, dBDy);
[REFx, REFy, ~] = backGroundCutOff(linearForce(REFx, REFy), tfmc.G', U);

if visible
    disp(mean((SFx - REFx) .^ 2 + (SFy - REFy) .^ 2, 'all') * 2);
    toc;
    
    figure;
    subplot(1, 2, 1);
    quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, REFx, REFy, 'DisplayName', 'estimated');
    title("Ridge");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
    subplot(1, 2, 2);
    ax = gca;
    plotROC(ax, divisionNumber, cellID, forceID, forceScale, tfmc.G, dBDx, dBDy, REFx, REFy)
end

%% Lasso

visible = true;

if visible
    tic;
end

[LEFx, LEFy, LEFMap] = TFMWithLASSO(divisionNumber, tfmc.G, dBDx, dBDy);
[LEFx, LEFy, ~] = backGroundCutOff(linearForce(LEFx, LEFy), tfmc.G', U);

if visible
    disp(mean((SFx - LEFx) .^ 2 + (SFy - LEFy) .^ 2, 'all') * 2);
    toc;
    
    figure;
    subplot(1, 2, 1);
    quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, LEFx, LEFy, 'DisplayName', 'estimated');
    title("Lasso");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
    subplot(1, 2, 2);
    ax = gca;
    plotROC(ax, divisionNumber, cellID, forceID, forceScale, tfmc.G, dBDx, dBDy, LEFx, LEFy)
end

%% Elastic Net

visible = true;

if visible
    tic;
end

[EEFx, EEFy, EEFMap] = TFMWithElasticNet(divisionNumber, tfmc.G, dBDx, dBDy, 1e-3);
[EEFx, EEFy, ~] = backGroundCutOff(linearForce(EEFx, EEFy), tfmc.G', U);

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

%% Bayes

visible = true;

if visible
    tic;
end

[BEFx, BEFy, BEFMap] = TFMWithBayesianEstimation(divisionNumber, tfmc.G, dBDx2, dBDy2, cellID, forceScale, 1, 1e0, gm);
[BEFx, BEFy, ~] = backGroundCutOff(linearForce(BEFx, BEFy), tfmc.G', U);

if visible
    disp(mean((SFx - BEFx) .^ 2 + (SFy - BEFy) .^ 2, 'all') * 2);
    toc;

    [~, ~, Pr, ~] = initial_estimation(cellID, forceScale);
    Pr = Pr';
    Prx = reshape(Pr(1:2:end), divisionNumber, divisionNumber);
    Pry = reshape(Pr(2:2:end), divisionNumber, divisionNumber);
    
    figure;
    subplot(1, 2, 1);
    quiver(X, Y, SFx, SFy, 'DisplayName', 'force');
    hold on; quiver(X, Y, BEFx, BEFy, 'DisplayName', 'estimated');
%     hold on; quiver(X, Y, Prx, Pry, 'DisplayName', 'prior');
    title("Bayes");
    xlim([0 w]);
    ylim([0 h]);
    daspect([1 1 1]);
    legend();
    subplot(1, 2, 2);
    ax = gca;
    plotROC(ax, divisionNumber, cellID, forceID, forceScale, tfmc.G, dBDx, dBDy, BEFx, BEFy)
end

rng('default');

%% Plot Results of Above All but EN.

maxM = max([SFMap, REFMap, LEFMap, BEFMap], [], 'all');
forceFiguresSBRL(X, Y, SFx, SFy, BEFx, BEFy, REFx, REFy, LEFx, LEFy, maxM);

%%
figure('WindowStyle', 'docked'); hold on;
scatter(IBLx, IBLy, 'DisplayName', 'Initial Bead Locations', 'MarkerFaceColor', 'b', 'MarkerEdgeAlpha', 0);
cell = cell_shape(divisionNumber, cellID);
cell = polyshape({cell(:, 1) * w}, {cell(:, 2) * h});
plot(cell);
title("Bead");
xlim([0 w]);
xticks([]);
ylim([0 h]);
yticks([]);
daspect([1 1 1]);

function F = linearForce(Fx, Fy)
Fx = reshape(Fx, 1, []);
Fy = reshape(Fy, 1, []);
F = reshape([Fx; Fy], 1, []);
end

function plotROC(ax, divisionNumber, cellID, forceID, forceScale, G, dBDx, dBDy, EFx, EFy)

EFx = reshape(EFx, 1, []);
EFy = reshape(EFy, 1, []);
ef = reshape([EFx; EFy], 1, []);
[tpfs, fpfs, ths] = tfmroc(ef, divisionNumber, cellID, forceID, forceScale);
[~, ~, th] = backGroundCutOff(ef, G', reshape([dBDx'; dBDy'], 1, []));
[~, I] = min(abs(ths - th));

ax.FontSize = 16;
hold on;
plot(ax, fpfs, tpfs, 'Marker', 'o');
scatter(ax, fpfs(I), tpfs(I), 'MarkerFaceColor', 'k');
title('ROC');
xlabel('FPF');
xlim([0, 1]);
ylabel('TPF');
ylim([0, 1]);
% legend('Location', 'eastoutside');
pbaspect([1 1 1]);
display(tfmauc(tpfs, fpfs))
end