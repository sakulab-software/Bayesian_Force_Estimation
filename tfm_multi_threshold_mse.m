function [thresholds, mses] = tfm_multi_threshold_mse(ef, division_number, cell_number, force_number, force_scale)

[tf, I] = cell_force_field(cell_number, force_number, force_scale);
% II = ones(1, division_number ^ 2);
% II(I) = 0;
% II = find(II);
[~, II, EI, ~] = cell_shape(division_number, cell_number);
II = union(II, EI);

tfx = tf(1:2:end);
tfy = tf(2:2:end);
tfnorm = sqrt(tfx(II) .^ 2 + tfy(II) .^ 2);
tfnorm = reshape(tfnorm, 1, []);
maxtfnorm = max(tfnorm);

efx = ef(1:2:end);
efx = reshape(efx, 1, []);
efy = ef(2:2:end);
efy = reshape(efy, 1, []);
efnorm = sqrt(efx(II) .^ 2 + efy(II) .^ 2);
efnorm = reshape(efnorm, 1, []);
% 
% [efnorm, III] = sort(efnorm);
% II = II(III);
% III = setdiff(II, I);

thresholds = efnorm(efnorm <= maxtfnorm);
thresholds = [thresholds, maxtfnorm];
thresholds = unique(thresholds);
mses = zeros(1, length(thresholds));

for i=1:length(thresholds)
    tmpefx = efx;
    tmpefx(efnorm <= thresholds(i)) = 0;
    tmpefy = efy;
    tmpefy(efnorm <= thresholds(i)) = 0;
    tmptfx = tfx;
    tmptfx(tfnorm <= thresholds(i)) = 0;
    tmptfy = tfy;
    tmptfy(tfnorm <= thresholds(i)) = 0;
    
    I = find((tmpefx .^ 2 + tmpefy .^ 2) ~= 0 | (tmptfx .^ 2 + tmptfy .^ 2) ~= 0);
    mses(i) = mean((tmpefx(I) - tmptfx(I)) .^ 2 + (tmpefy(I) - tmptfy(I)) .^ 2);
%     n = find(tmpefx .^ 2 + tmpefy .^ 2);
%     n = union(n, find(tmptfx .^ 2 + tmptfy .^ 2));
%     n = unique(n);
%     n = numel(n);
%     mses(i) = sum((tmpefx - tmptfx) .^ 2 + (tmpefy - tmptfy) .^ 2) / 2 / n;
end

end