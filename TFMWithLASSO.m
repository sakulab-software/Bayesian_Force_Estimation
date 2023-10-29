function [EFx, EFy, EFMap] = TFMWithLASSO(divisionNumber, G, BDx, BDy)

y = reshape([BDx, BDy]', 1, [])';
h1 = hyperparameters('fitrlinear', G, y);
h2.ShowPlots = false;
h2.Verbose = 0;
[ef, ~] = fitrlinear(G, reshape([BDx, BDy]', 1, [])', ...
    'Regularization', 'lasso', 'Learner', 'leastsquares', 'FitBias', false, ...
    'OptimizeHyperparameters', h1(1), 'HyperparameterOptimizationOptions', h2, ...
    'Solver', 'sparsa');
ef = ef.Beta;
% ef = ef(2:end);
EFx = ef(1:2:end);
EFy = ef(2:2:end);
EFx = reshape(EFx, divisionNumber, divisionNumber);
EFy = reshape(EFy, divisionNumber, divisionNumber);
EFMap = sqrt(EFx .^ 2 + EFy .^ 2);
end
