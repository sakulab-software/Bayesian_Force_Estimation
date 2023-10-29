function [F, I] = cell_force_field(cell_shape, number, power)

path = fullfile(cd, 'cell_force_field', num2str(cell_shape), [num2str(number) '.mat']);
m = load(path);
F = m.f * power;
I = m.i;

end