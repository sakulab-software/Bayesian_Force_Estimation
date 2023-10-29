function mse = tfmroc(ef, division_number, cell_number, force_number, force_scale)

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

efx = ef(1:2:end);
efy = ef(2:2:end);
efnorm = sqrt(efx(II) .^ 2 + efy(II) .^ 2);
efnorm = reshape(efnorm, 1, []);

mse = mean((efx - tfx) .^ 2 + (efy - tfy) .^ 2) / mean(tfx .^ 2 + tfy .^ 2);
end

% function [mtpfs, mfpfs] = tfmroc(ef, division_number, cell_number, force_number, force_scale)
% 
% [tf, I] = cell_force_field(cell_number, force_number, force_scale);
% % II = ones(1, division_number ^ 2);
% % II(I) = 0;
% % II = find(II);
% [~, II, ~, ~] = cell_shape(division_number, cell_number);
% 
% tfx = tf(1:2:length(tf));
% tfy = tf(2:2:length(tf));
% tfn = sqrt(tfx(I) .^ 2 + tfy(I) .^ 2);
% tfn = reshape(tfn, 1, []);
% mtfn = mean(tfn);
% tfa = angle(complex(tfx(I), tfy(I)));
% btfa = angle(complex(tfx(II), tfy(II)));
% 
% efx = ef(1:2:length(ef));
% efy = ef(2:2:length(ef));
% efn = sqrt(efx(I) .^ 2 + efy(I) .^ 2);
% efn = reshape(efn, 1, []);
% befn = sqrt(efx(II) .^ 2 + efy(II) .^ 2);
% efa = angle(complex(efx(I), efy(I)));
% efa = reshape(efa, 1, []);
% befa = angle(complex(efx(II), efy(II)));
% befa = reshape(befa, 1, []);
% 
% 
% 
% thresholds = 0:0.005:2;
% mtpfs = zeros(1, length(thresholds));
% mfpfs = zeros(1, length(thresholds));
% 
% for i=1:length(thresholds)
%     mtpfs(i) = length(find(efn >= tfn * thresholds(i))) / length(I);
%     mfpfs(i) = length(find(befn >= mtfn * thresholds(i))) / length(II);
% end
% end