%%

pRB = 'mainBayes210413-001502.mat';
pRL = 'mainLasso210411-003148.mat';
pRR = 'mainRidge210411-143931.mat';

%%
RB = load(pRB);
RL = load(pRL);
RR = load(pRR);

%%
divisionNumber = 30;

b = 10;
c = 1;
f = 1;
n = 1;
l = 1;
t = 1;

X = linspace(0, 15, divisionNumber);
Y = linspace(0, 15, divisionNumber);
[X, Y] = meshgrid(X, Y);

[tf, ~] = cell_force_field(RB.cellIDList(c), RB.forceIDList(f), RB.forceScale);
SFx = reshape(tf(1:2:end), divisionNumber, divisionNumber);
SFy = reshape(tf(2:2:end), divisionNumber, divisionNumber);
SFMap = sqrt(SFx .^ 2 + SFy .^ 2);

EF = RB.estimated_force_cell_matrix{b, c, f, n, l, t};
BEFx = reshape(EF(1:2:end), divisionNumber, divisionNumber);
BEFy = reshape(EF(2:2:end), divisionNumber, divisionNumber);
BEFMap = sqrt(BEFx .^ 2 + BEFy .^ 2);

EF = RL.estimated_force_cell_matrix{b, c, f, n, t};
LEFx = reshape(EF(1:2:end), divisionNumber, divisionNumber);
LEFy = reshape(EF(2:2:end), divisionNumber, divisionNumber);
LEFMap = sqrt(LEFx .^ 2 + LEFy .^ 2);

EF = RR.estimated_force_cell_matrix{b, c, f, n, t};
REFx = reshape(EF(1:2:end), divisionNumber, divisionNumber);
REFy = reshape(EF(2:2:end), divisionNumber, divisionNumber);
REFMap = sqrt(REFx .^ 2 + REFy .^ 2);

maxM = max(cat(3, SFMap, REFMap, LEFMap, BEFMap), [], 'all');
forceFiguresSBRL(X, Y, SFx, SFy, BEFx, BEFy, REFx, REFy, LEFx, LEFy, maxM);

SFMap = SFMap(:);
BEFMap = BEFMap(:);
LEFMap = LEFMap(:);
REFMap = REFMap(:);

edges = linspace(0, maxM, 50);
upperBound = 0.05;

subplot(4, 1, 1);
gca.YScale = 'log';
c = histcounts(SFMap, edges) / length(SFMap);
bar(edges(1:end-1), c);
title('Synthetic force');
xlim([0 maxM]);
ylim([0 upperBound]);

subplot(4, 1, 2);
gca.YScale = 'log';
c = histcounts(BEFMap, edges) / length(BEFMap);
bar(edges(1:end-1), c);
title('Bayes');
xlim([0 maxM]);
ylim([0 upperBound]);

subplot(4, 1, 3);
gca.YScale = 'log';
c = histcounts(REFMap, edges) / length(REFMap);
bar(edges(1:end-1), c);
title('Ridge');
xlim([0 maxM]);
ylim([0 upperBound]);

subplot(4, 1, 4);
gca.YScale = 'log';
c = histcounts(LEFMap, edges) / length(LEFMap);
bar(edges(1:end-1), c);
title('Lasso');
xlim([0 maxM]);
ylim([0 upperBound]);
