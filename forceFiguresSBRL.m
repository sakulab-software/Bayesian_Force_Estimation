function forceFiguresSBRL(X, Y, SFx, SFy, BFx, BFy, RFx, RFy, LFx, LFy, maxM)

FxList = {SFx, BFx, RFx, LFx};
FyList = {SFy, BFy, RFy, LFy};
SqFList = cellfun(@(x, y) sqrt(x .^ 2 + y .^ 2), FxList, FyList, UniformOutput=false);

if ~exist("maxM", 'var') || numel(maxM) ~= 4
    maxM = cellfun(@(x) max(x, [], "all"), SqFList);
end

titleList = [ ...
    "Synthetic Force", ... 
    "Bayes", ...
    "Ridge", ...
    "Lasso"];

figure;
tiledlayout(2, 2);
for i=1:4
    nexttile;
    forceMagnitudeMapPlot(X, Y, SqFList{i}, maxM(i));
    title(titleList(i));
end

figure;
tiledlayout(2, 2);
for i=1:4
    nexttile;
    if ~exist('ind', 'var')
        ind = forceDirectionPlot(X, Y, FxList{i}, FyList{i}, SqFList{i}, 0);
    end
    forceDirectionPlot(X, Y, FxList{i}, FyList{i}, SqFList{i}, 0);
    title(titleList(i));
end

end


