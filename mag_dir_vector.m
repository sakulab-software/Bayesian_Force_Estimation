function [M, D] = mag_dir_vector(V, dim)

M = zeros([1 length(V) / dim]);
D = zeros(size(V));

for i=1:length(M)
    sum = 0;
    for j=1:dim
        sum = sum + V(dim * (i - 1) + j) ^ 2;
    end
    M(i) = sqrt(sum);
    
    for j=1:dim
        D(dim * (i - 1) + j) = V(dim * (i - 1) + j) / M(i);
    end
end

end