load('mainBayes210106-210838.mat');
% b c f n l t

%%
n = divisionNumber;
m = n;
X = linspace(0, w, n);
Y = linspace(0, h, m);
[X, Y] = meshgrid(X, Y);

%%
angMatrix = zeros(size(estimated_force_cell_matrix));
n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpAng = zeros(size(angMatrix(b, c, f, :, :, :)));
            tmpEst = estimated_force_cell_matrix(b, c, f, :, :, :);
            tmpIBL = initial_bead_location_cell_matrix(b, c, f, :, :, :);
            tmpBD = bead_displacement_cell_matrix(b, c, f, :, :, :);
            
            for i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                [SF, I] = cell_force_field(cellIDList(c), forceIDList(f), forceScale);
                SFx = SF(2 * I - 1);
                SFy = SF(2 * I);
                SFAng = angle(SFx + 1i * SFy);
                SFAng = SFAng(:);
                SFAng = exp(1i * SFAng);
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
                EFx = EFx(I);
                EFy = EFy(I);
                EFAng = angle(EFx + 1i * EFy);
                EFAng = EFAng(:);
                EFAng = exp(1i * EFAng);
                
                tmpAng(i) = angle(sum(exp(1i * acos(real(SFAng) .* real(EFAng) + imag(SFAng) .* imag(EFAng)))));
            end
            angMatrix(b, c, f, :, :, :) = tmpAng;
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

M = angMatrix;

SMB = zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));
SUBB = zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));
SLBB = zeros(length(beadNumberList), length(lamList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for l=1:length(lamList)
        for r=1:length(noiseRadiusList)
            SMB(b, l, r) = mean(M(b, :, :, r, l, :), 'all');
            SUBB(b, l, r) = std(M(b, :, :, r, l, :), [], 'all');
            SLBB(b, l, r) = SUBB(b, l, r);
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
% ylim([0, 1]);
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
% ylim([0, 1]);
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
% clear
load('mainRidge210107-151954.mat');
% b c f n t

%%
angMatrix = zeros(size(estimated_force_cell_matrix));

n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpAng = zeros(size(angMatrix(b, c, f, :, :)));
            tmpEst = estimated_force_cell_matrix(b, c, f, :, :);
            tmpIBL = initial_bead_location_cell_matrix(b, c, f, :, :);
            tmpBD = bead_displacement_cell_matrix(b, c, f, :, :);
            
            parfor i=1:numel(estimated_force_cell_matrix(b, c, f, :, :))
                [SF, I] = cell_force_field(cellIDList(c), forceIDList(f), forceScale);
                SFx = SF(2 * I - 1);
                SFy = SF(2 * I);
                SFAng = angle(SFx + 1i * SFy);
                SFAng = SFAng(:);
                SFAng = exp(1i * SFAng);
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
                EFx = EFx(I);
                EFy = EFy(I);
                EFAng = angle(EFx + 1i * EFy);
                EFAng = EFAng(:);
                EFAng = exp(1i * EFAng);
                
                tmpAng(i) = angle(sum(exp(1i * acos(real(SFAng) .* real(EFAng) + imag(SFAng) .* imag(EFAng)))));
            end
            angMatrix(b, c, f, :, :) = tmpAng;
        end
    end
end

%%

M = angMatrix;

SMR = zeros(length(beadNumberList), length(noiseRadiusList));
SUBR = zeros(length(beadNumberList), length(noiseRadiusList));
SLBR = zeros(length(beadNumberList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for r=1:length(noiseRadiusList)
        SMR(b, r) = mean(M(b, :, :, r, :), 'all');
        SUBR(b, r) = std(M(b, :, :, r, :), [], 'all');
        SLBR(b, r) = SUBR(b, r);
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

M = angMatrix;

SML = zeros(length(beadNumberList), length(noiseRadiusList));
SUBL = zeros(length(beadNumberList), length(noiseRadiusList));
SLBL = zeros(length(beadNumberList), length(noiseRadiusList));

for b=1:length(beadNumberList)
    for r=1:length(noiseRadiusList)
        SML(b, r) = mean(M(b, :, :, r, :), 'all');
        SUBL(b, r) = std(M(b, :, :, r, :), [], 'all');
        SLBL(b, r) = SUBL(b, r);
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
SMB = SMB * 180 / pi;
SMR = SMR * 180 / pi;
SML = SML * 180 / pi;
SUBB = SUBB * 180 / pi;
SUBR = SUBR * 180 /pi;
SUBL = SUBL * 180 / pi;
SLBB = SLBB * 180 / pi;
SLBR = SLBR * 180 / pi;
SLBL = SLBL * 180 / pi;
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
