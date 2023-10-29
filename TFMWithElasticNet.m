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
