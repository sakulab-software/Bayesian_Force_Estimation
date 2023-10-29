function V = varianceInnerCell(ef, division_number, cell_number, force_number, force_scale)

if iscell(ef)

    [~, I] = cell_force_field(cell_number, force_number, force_scale);
    % II = ones(1, division_number ^ 2);
    % II(I) = 0;
    % II = find(II);
    [~, II, EI, ~] = cell_shape(division_number, cell_number);
    II = union(II, EI);
    
    efx = cellfun(@(e) e(1:2:end), ef, UniformOutput=false);
    efy = cellfun(@(e) e(2:2:end), ef, UniformOutput=false);
    efnorm = cellfun(@(x, y) sqrt(x .^ 2 + y .^ 2), efx, efy, UniformOutput=false);

    V = zeros(1, division_number ^ 2);
    for i=1:numel(II)
        II_efnorm = cellfun(@(e) e(II(i)), efnorm);
        V(II(i)) = std(II_efnorm, 0, "all") / mean(II_efnorm, "all");
        if all(II_efnorm == 0)
            V(II(i)) = 0;
        end
    end
else
    V = nan;
end

end