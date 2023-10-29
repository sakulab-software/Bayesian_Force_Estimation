
function [th, thBDx, thBDy] = beadDispThreshold(BD, BDx, BDy, t, visible)

if ~exist("t", 'var')
    t = 2;
end

if ~exist("visible", 'var')
    visible = false;
end

[counts, edges] = histcounts(BD, 10);
x = mean([edges(1:end-1); edges(2:end)]);

fun = @(coef, xdata) coef(1) * exp(coef(2) * xdata);
coef0 = [1, -1];
coef = lsqcurvefit(fun, coef0, x, counts);

th = - log(t) / coef(2);
ind = BD < th;

thBDx = BDx;
thBDy = BDy;
thBDx(ind) = 1e-12;
thBDy(ind) = 1e-12;

if visible
    figure; hold on;
    stem(x, counts);
    line(ones(1, 2) * th, [0, max(counts)]);
    plot(x, fun(coef, x));
end

end