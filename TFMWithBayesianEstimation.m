function [EFx, EFy, EFMap] = TFMWithBayesianEstimation(divisionNumber, G, BDx, BDy, cellID, forceScale, A, B, beadDensityMap, thr, X, Y, IBLx, IBLy)

BD = reshape([BDx, BDy]', 1, [])';
% [~, ~, Prior, II] = initial_estimation(cellID, forceScale);
[Prior, II] = priorFromCellID(cellID, forceScale);
Prior = Prior';

if exist("beadDensityMap", 'var')
    Prior(1:2:end) = Prior(1:2:end);% .* reshape(beadDensityMap, [], 1);
    Prior(2:2:end) = Prior(2:2:end);% .* reshape(beadDensityMap, [], 1);
end

III = reshape([2 * II - 1; 2 * II], 1, []);

bayestfm = BayesTFM2D(A, B, thr, 10000);
bayestfm = bayestfm.fit(BD, G(:, III), Prior(III), X(II), Y(II), IBLx, IBLy);
EF = zeros(1, 2 * divisionNumber ^ 2);
EF(III) = bayestfm.force;

EFx = reshape(EF(1:2:end), divisionNumber, divisionNumber);
EFy = reshape(EF(2:2:end), divisionNumber, divisionNumber);

EFMap = sqrt(EFx .^ 2 + EFy .^ 2);
end
