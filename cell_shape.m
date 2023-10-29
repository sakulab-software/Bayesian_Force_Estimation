function [CSPs, II, EI, TFin] = cell_shape(division_number, cell_number)
% [CSPs, II, EI, TFin] = cell_shape(divisionNumber, cellID)
% CSPs = Vertices
% II = Indices of grid points inside of the cell
% EI = on edge
% TFin

m = load(fullfile(cd, 'cell_shape', ['force_field_' int2str(cell_number) '.mat']));
CSPs = m.cs.Vertices;
II = find(m.II == 1);
EI = find(m.EI == 1);
TFin = false([1 division_number ^ 2]);
TFin(II) = true;
% 
% GP = grid_points(division_number);
% 
% MAT = open(fullfile(cd, "cell_shape", ['cell' int2str(cell_number) '.mat']));
% cell = MAT.cell;
% CSPs = cell;
% cell = polyshape({cell(:, 1)}, {cell(:, 2)});
% 
% TFin = isinterior(cell, GP(1:2:length(GP)), GP(2:2:length(GP)));
% 
% I = zeros([1 division_number ^ 2]);
% for i=1:division_number ^ 2
%     if TFin(i)
%         I(i) = 1;
%     end
% end
% IMG = reshape(I, division_number, []);
% II = find(IMG);
% II = reshape(II, 1, []);
% EI = edge(IMG);
% EI = find(EI);
% EI = reshape(EI, 1, []);

% IMG2 = zeros([division_number division_number]);
% for i=1:division_number
%     for j=1:division_number
%         if ~isempty(find(II == division_number * (i - 1) + j, 1))
%             IMG2(j, i) = 1;
%         end
%     end
% end
% figure, imshow(IMG2);
end