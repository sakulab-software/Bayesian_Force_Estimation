load('mainBayes210106-210838.mat');
% b c f n l t

%%
n = divisionNumber;
m = n;
X = linspace(0, w, n);
Y = linspace(0, h, m);
[X, Y] = meshgrid(X, Y);

%%

AUCMatrix = zeros(size(estimated_force_cell_matrix));
n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpAUC = zeros(size(AUCMatrix(b, c, f, :, :, :)));
            tmpEst = estimated_force_cell_matrix(b, c, f, :, :, :);
            tmpIBL = initial_bead_location_cell_matrix(b, c, f, :, :, :);
            tmpBD = bead_displacement_cell_matrix(b, c, f, :, :, :);
            
            parfor i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                EF = tmpEst{i};
                EFx = EF(1:2:end);
                EFy = EF(2:2:end);
                EFx = reshape(EFx, n, m);
                EFy = reshape(EFy, n, m);
                IBL = tmpIBL{i};
                IBLx = IBL(1:2:end);
                IBLy = IBL(2:2:end);
                Sim = TFM_computation(X, Y, EFx, EFy, IBLx(:), IBLy(:));
                Sim = Sim.simulate();
                BD = tmpBD{i};
                [EFx, EFy, ~] = backGroundCutOff(EF, Sim.G', BD);
                EF = reshape([reshape(EFx, 1, []); reshape(EFy, 1, [])], [], 1);
                [tpfs, fpfs] = tfmroc(EF, divisionNumber, cellIDList(c), forceIDList(f), forceScale);
                tmpAUC(i) = tfmauc(tpfs, fpfs);
            end
            AUCMatrix(b, c, f, :, :, :) = tmpAUC;
        end
    end
end

%% 
% tmp = auc_matrix(5, 5, 1, 6, :, :);
% tmp = tmp(:);
% plot(tmp);
% mean(tmp)

%%

% ibl = initial_bead_location_cell_matrix{8, 5, 1, 1, 4, 1}(:);
% iblx = ibl(1:2:length(ibl));
% ibly = ibl(2:2:length(ibl));
% bd = bead_displacement_cell_matrix{8, 5, 1, 1, 4, 1}(:);
% bdx = bd(1:2:length(bd));
% bdy = bd(2:2:length(bd));
% scatter(iblx, ibly);
% hold on;
% quiver(iblx, ibly, bdx, bdy);

%%
% ef = estimated_force_cell_matrix{15, 5, 1, 6, 4, 1}(:);
% efx = ef(1:2:length(ef));
% efy = ef(2:2:length(ef));
% quiver(X, Y, efx, efy);

%%

M = AUCMatrix;

SMB= zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));
SUBB = zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));
SLBB = zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for l=1:length(lamList)
        for r=1:length(noiseRadiusList)
            ab = betafit(reshape(round(M(b, :, :, r, l, :), 4), 1, []));
            SMB(b, l, r) = ab(1) / (ab(1) + ab(2));
            SUBB(b, l, r) = betaincinv(0.95, ab(1), ab(2));
            SLBB(b, l, r) = betaincinv(0.05, ab(1), ab(2));
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=3:length(lamList)
    errorbar((beadNumberList + 5 * (l - 1)) / 255, ...
        SMB(:, l, 1), ...
        SLBB(:, l, 1) - SMB(:, l, 1), ...
        SUBB(:, l, 1) - SMB(:, l, 1), ...
        'DisplayName', "log\lambda=" + lamList(l));
end
title('AUC');
xlabel('Bead density [-]');
ylabel('AUC');
ylim([0, 1]);
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=4:length(noiseRadiusList) - 2
    errorbar((beadNumberList + 5 * (r - 7 - 1)) / 255, ...
        SMB(:, 4, r), ...
        SLBB(:, 4, r) - SMB(:, 4, r), ...
        SUBB(:, 4, r) - SMB(:, 4, r), ...
        'DisplayName', "noise radius=" + noiseRadiusList(r));
end
title('AUC');
xlabel('Bead density [-]');
ylabel('AUC');
ylim([0, 1]);
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
% clear;
load('mat/main_Ridge200929-164743.mat');
% b c f n t

%%

AUCMatrix = zeros(size(estimated_force_cell_matrix));
n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpAUC = zeros(size(AUCMatrix(b, c, f, :, :)));
            tmpEst = estimated_force_cell_matrix(b, c, f, :, :);
            tmpIBL = initial_bead_location_cell_matrix(b, c, f, :, :);
            tmpBD = bead_displacement_cell_matrix(b, c, f, :, :);
            
            parfor i=1:numel(estimated_force_cell_matrix(b, c, f, :, :))
                EF = tmpEst{i};
                EFx = EF(1:2:end);
                EFy = EF(2:2:end);
                EFx = reshape(EFx, n, m);
                EFy = reshape(EFy, n, m);
                IBL = tmpIBL{i};
                IBLx = IBL(1:2:end);
                IBLy = IBL(2:2:end);
                Sim = TFM_computation(X, Y, EFx, EFy, IBLx(:), IBLy(:));
                Sim = Sim.simulate();
                BD = tmpBD{i};
                [EFx, EFy, ~] = backGroundCutOff(EF, Sim.G', BD);
                EF = reshape([reshape(EFx, 1, []); reshape(EFy, 1, [])], [], 1);
                [tpfs, fpfs] = tfmroc(EF, divisionNumber, cellIDList(c), forceIDList(f), forceScale);
                tmpAUC(i) = tfmauc(tpfs, fpfs);
            end
            AUCMatrix(b, c, f, :, :) = tmpAUC;
        end
    end
end

%%

M = AUCMatrix;

SMR = zeros(length(beadNumberList), length(noiseRadiusList));
SUBR = zeros(length(beadNumberList), length(noiseRadiusList));
SLBR = zeros(length(beadNumberList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for r=1:length(noiseRadiusList)
        if std(reshape(M(b, :, :, r, :), 1, [])) ~= 0
            ab = betafit(reshape(round(M(b, :, :, r, :), 4), 1, []));
            SMR(b, r) = ab(1) / (ab(1) + ab(2));
            SUBR(b, r) = betaincinv(0.95, ab(1), ab(2));
            SLBR(b, r) = betaincinv(0.05, ab(1), ab(2));
        else
            SMR(b, r) = mean(M(b, :, :, r, :), 'all');
            SUBR(b, r) = 0;
            SLBR(b, r) = 0;
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
l = 1;
errorbar((beadNumberList + 5 * (l - 1)) / 255, ...
    SMR(:, 1), ...
    SLBR(:, 1) - SMR(:, 1), ...
    SUBR(:, 1) - SMR(:, 1));
title('AUC');
xlabel('Bead density [-]');
ylabel('AUC');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);


%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=5:length(noiseRadiusList)
    errorbar((beadNumberList + 5 * (r - 7 - 1)) / 255, ...
        SMR(:, r), ...
        SLBR(:, r) - SMR(:, r), ...
        SUBR(:, r) - SMR(:, r), ...
        'DisplayName', "noise radius=" + noiseRadiusList(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
clear;
load('mat/main_Ridge200930-163904.mat'); %Lasso
% b c f n t

%%
auc_matrix = zeros(size(estimated_force_cell_matrix));

n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpauc = zeros(size(auc_matrix(b, c, f, :, :, :)));
            tmpest = estimated_force_cell_matrix(b, c, f, :, :, :);
            for i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                ef = tmpest{i};
                [tpfs, fpfs] = tfmroc(ef, divisionNumber, cellIDList(c), forceIDList(f), forceScale);
                tmpauc(i) = tfmauc(tpfs, fpfs);
            end
            auc_matrix(b, c, f, :, :, :) = tmpauc;
        end
    end
end

%%

M = AUCMatrix;

SML = zeros(length(beadNumberList), length(noiseRadiusList));
SUBL = zeros(length(beadNumberList), length(noiseRadiusList));
SLBL = zeros(length(beadNumberList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for r=1:length(noiseRadiusList)
        if std(reshape(M(b, :, :, r, :), 1, [])) ~= 0
            ab = betafit(reshape(round(M(b, :, :, r, :), 4), 1, []));
            SML(b, r) = ab(1) / (ab(1) + ab(2));
            SUBL(b, r) = betaincinv(0.95, ab(1), ab(2));
            SLBL(b, r) = betaincinv(0.05, ab(1), ab(2));
        else
            SML(b, r) = mean(M(b, :, :, r, :), 'all');
            SUBL(b, r) = 0;
            SLBL(b, r) = 0;
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
errorbar(beadNumberList / 255, ...
    SML(:, 1), ...
    SLBL(:, 1) - SML(:, 1), ...
    SUBL(:, 1) - SML(:, 1), ...
    'DisplayName', "without noise");
title('AUC');
xlabel('Bead density [-]');
ylabel('AUC');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);


%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=1:length(noiseRadiusList)
    errorbar((beadNumberList + 5 * (r - 7 - 1)) / 255, ...
        SML(:, r), ...
        SLBL(:, r) - SML(:, r), ...
        SUBL(:, r) - SML(:, r), ...
        'DisplayName', "noise radius" + noiseRadiusList(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%

load('mat/RidgeAUC200718-210820.mat');
SMR = SM;
SUBR = SUB;
SLBR = SLB;
clear SM SUB SLB;
load('mat/LassoAUC200718-192141.mat');
SML = SM;
SUBL = SUB;
SLBL = SLB;
clear SM SUB SLB;
load('BayesAUC.mat');
SMB = SM;
SUBB = SUB;
SLBB = SLB;
clear SM SUB SLB;

%%
figure; hold on;
ax = gca;
ax.FontSize = 16;
errorbar((beadNumberList) / 255, ...
    SMB(:, 4, 1), ...
    SLBB(:, 4, 1) - SMB(:, 4, 1), ...
    SUBB(:, 4, 1) - SMB(:, 4, 1), ...
    'DisplayName', 'Bayes');
errorbar((beadNumberList + 5) / 255, ...
    SMR(:, 4, 1), ...
    SLBR(:, 4, 1) - SMR(:, 4, 1), ...
    SUBR(:, 4, 1) - SMR(:, 4, 1), ...
    'DisplayName', 'Ridge');
errorbar((beadNumberList + 10) / 255, ...
    SML(:, 3, 1), ...
    SLBL(:, 3, 1) - SML(:, 3, 1), ...
    SUBL(:, 3, 1) - SML(:, 3, 1), ...
    'DisplayName', 'Lasso');
title('AUC');
xlabel('Bead density [\mum^{-2}]');
ylabel('AUC');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
figure; hold on;
ax = gca;
ax.FontSize = 16;
errorbar((beadNumberList) / 255, ...
    SMB(:, 4, 7), ...
    SLBB(:, 4, 7) - SMB(:, 4, 7), ...
    SUBB(:, 4, 7) - SMB(:, 4, 7), ...
    'DisplayName', 'Bayes');
errorbar((beadNumberList + 5) / 255, ...
    SMR(:, 4, 7), ...
    SLBR(:, 4, 7) - SMR(:, 4, 7), ...
    SUBR(:, 4, 7) - SMR(:, 4, 7), ...
    'DisplayName', 'Ridge');
errorbar((beadNumberList + 10) / 255, ...
    SML(:, 3, 7), ...
    SLBL(:, 3, 7) - SML(:, 3, 7), ...
    SUBL(:, 3, 7) - SML(:, 3, 7), ...
    'DisplayName', 'Lasso');
title('AUC');
xlabel('Bead density [\mum^{-2}]');
ylabel('AUC');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%

load BayesAUC200929-220023.mat
load RidgeAUC200929-164743.mat
load LassoAUC200930-163904.mat

%%
figure; hold on;
ax = gca;
ax.FontSize = 16;
ind = 1:7;
noiseLevel = 8;
e = errorbar((beadNumberList(ind)) / 255, ...
    SMB(ind, 4, noiseLevel), ...
    SLBB(ind, 4, noiseLevel) - SMB(ind, 4, noiseLevel), ...
    SUBB(ind, 4, noiseLevel) - SMB(ind, 4, noiseLevel), ...
    'DisplayName', 'Bayes');
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;

e = errorbar((beadNumberList(ind) + 1) / 255, ...
    SMR(ind, noiseLevel), ...
    SLBR(ind, noiseLevel) - SMR(ind, noiseLevel), ...
    SUBR(ind, noiseLevel) - SMR(ind, noiseLevel), ...
    'DisplayName', 'Ridge');
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;

e = errorbar((beadNumberList(ind) + 2) / 255, ...
    SML(ind, noiseLevel), ...
    SLBL(ind, noiseLevel) - SML(ind, noiseLevel), ...
    SUBL(ind, noiseLevel) - SML(ind, noiseLevel), ...
    'DisplayName', 'Lasso');
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;

xlabel('Bead density [\mum^{-2}]');
xlim([0.07 0.33]);
ylabel('AUC');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);
