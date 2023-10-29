function [M, A] = mag_dir_map(division_number, V, regularization, magmax, magmin)

[M, A] = mag_dir_vector(V, 2);

A = angle(A(1:2:length(A)) + A(2:2:length(A)) * 1i);
A(isnan(A)) = 0;
A = reshape(A, division_number, []);

if ~exist('magmax', 'var')
    magmax = max(M);
end

if ~exist('magmin', 'var')
    magmin = min(M);
end

M = reshape(M, division_number, []);

if exist('regularization', 'var') && regularization
    M = (M + magmin) ./ (magmax - magmin);
    A = (A + pi) ./ (2 * pi);
end

end