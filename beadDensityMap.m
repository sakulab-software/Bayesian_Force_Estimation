function gm = beadDensityMap(X, Y, BeadNumber, IBLx, IBLy, divisionNumber, cellID)

gm = local_bead_density(BeadNumber, IBLx, IBLy, 1, divisionNumber, X, Y);
[~, II, ~, ~] = cell_shape(divisionNumber, cellID);
gmmax_inside = max(gm(II), [], 'all');
gmmin = min(reshape(gm, 1, []));
gm = (gm - gmmin) ./ (gmmax_inside - gmmin);
% sigmoid transformation to saturate gm over 0.5 
gm = 1 ./ (1 + exp(- 10 * (gm - 0.5)));

end
