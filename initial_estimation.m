function [F, M, IE, I] = initial_estimation(cell_number, scale)

path = fullfile(cd, 'cell_shape', ['force_field_' num2str(cell_number) '.mat']);
data = load(path);
M = exp(- 1 / (2 * 4) * (data.d) .^ 2);
M = M ./ max(M, [], 'all');
M = scale .* M;
F = data.f;

I = find(data.EI == 1);

I = union(I, find(data.II == 1));

IE = ones([1, length(F)]) * 1e-12;

for i=I
    IE(2 * i - 1) = M(i) * F(2 * i - 1);
    IE(2 * i) = M(i) * F(2 * i);
end

end