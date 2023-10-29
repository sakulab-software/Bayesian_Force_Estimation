
cellID = 4;
path = num2str(cellID);
mkdir(path);

for j=1:100
    [f, iI] = priorFromCellID(cellID, 1);
    i = iI(randperm(length(iI), 20));
    ii = reshape([2 * i - 1; 2 * i], 1, []);
    tmp = f(ii);
    f(:) = 0;
    r = rand(1, 20) * 1.9 + 0.1;
    tmp(1:2:end) = tmp(1:2:end) .* r;
    tmp(2:2:end) = tmp(2:2:end) .* r;
    theta = rand(1, 20) * 60 - 30;
    theta = zeros(1, 20);
    for k=1:20
        tmp(2*k-1:2*k) = rotVec(tmp(2*k-1:2*k), deg2rad(theta(k)));
    end
    f(ii) = tmp;
    save(fullfile(path, num2str(j) + ".mat"), "f", "i");
end

function vec = rotVec(vec, theta)
    s = sin(theta);
    c = cos(theta);
    vec(1) = c * vec(1) - s * vec(2);
    vec(2) = s * vec(1) + c * vec(2);
end