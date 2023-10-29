function [Prior, I] = priorFromCellID(cellID, forceScale)
path = fullfile(cd, 'cell_shape', ['cell' num2str(cellID) '.mat']);
data = load(path);
Px = data.Px;
Py = data.Py;
Px = Px(end:-1:1, :);
Py = -Py(end:-1:1, :);
Prior = reshape([reshape(Px, 1, []); ...
                 reshape(Py, 1, [])], 1, []);
Prior = Prior * forceScale;
I = find(Px .^ 2 + Py .^ 2 > 0);
end