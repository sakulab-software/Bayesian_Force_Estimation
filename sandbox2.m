%%
clear

cs = cell_shape(30, 1);

X = linspace(0, 1, 90);
Y = linspace(0, 1, 90);
[X, Y] = meshgrid(X, Y);
xq = X(:);
yq = Y(:);

[in, on] = inpolygon(xq, yq, cs(:, 1), cs(:, 2));

I = false(size(X));
I(in) = true;
I(on) = true;
e = edge(I);

ix = xq(in);
iy = yq(in);
ex = xq(e(:));
ey = yq(e(:));


%%

figure; hold on;
plot(cs(:, 1), cs(:, 2));
scatter(xq(in), yq(in));
scatter(xq(on), yq(on));
scatter(ex, ey);

xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]);

%%

figure; hold on;
[a, b] = size(I);
plot([0, b], [0, a]);
imagesc(e);
xlim([0 b]);
ylim([0 a]);
pbaspect([1 1 1]);


%%
clear

param = [10, 100, 0.1, 1, 0.001];

tspan = [0 10];
y0 = [1; 1; 5e2];
[t, y] = ode45(@(t, y) modelAB(y, param), tspan, y0);

sim_dt = t(2:end) - t(1:end-1);
sim_dAdt = y(:, 1);
sim_dAdt = sim_dAdt(2:end) - sim_dAdt(1:end-1);
sim_dAdt = sim_dAdt ./ sim_dt;
sim_dBdt = y(:, 2);
sim_dBdt = sim_dBdt(2:end) - sim_dBdt(1:end-1);
sim_dBdt = sim_dBdt ./ sim_dt;

dAdt = zeros(size(sim_dAdt));
dBdt = zeros(size(sim_dBdt));

for i=1:length(dAdt)
    dAdt(i) = modelA(y(i, :), param(1), param(2));
    dBdt(i) = modelB(y(i, :), param(3), param(4), param(5));
end
%%

figure; hold on;
plot(t, y(:, 1), 'DisplayName', 'A');
plot(t, y(:, 2), 'DisplayName', 'B');
xlim([t(1) t(end)]);
% ylim([0 1]);
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%

figure; hold on;
plot([-100 100], [-100 100], 'DisplayName', 'ref');
% scatter(dAdt, sim_dAdt, 'DisplayName', 'A');
plot(dBdt, sim_dBdt, 'DisplayName', 'B');
xlabel('theoretical');
ylabel('empirical');
legend;
pbaspect([1 1 1]);

%%

figure; hold on;
plot(t(1:end-1), sim_dAdt, 'DisplayName', 'empirical dAdt');
plot(t(1:end-1), sim_dBdt, 'DisplayName', 'empirical dBdt');
xlim([t(1) t(end)]);
% ylim([0 1]);
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
clear
param = [10, 100, 0.1, 1, 0.001];
B = 0:0.0001:1.011;
I = [0:10 20 50 100];
dBdt = zeros(length(I), length(B));

for i=1:length(I)
    for j=1:length(B)
        dBdt(i, j) = modelB([0, B(j), I(i)], param(3), param(4), param(5));
    end
end

figure; hold on;
for i=1:length(I)
    plot(B, dBdt(i, :));
end

%%

modelB(y(1, :), param(3), param(4), param(5))
modelB(y(2, :), param(3), param(4), param(5))
modelB(y(3, :), param(3), param(4), param(5))

function dydt = modelAB(y, param)
grA = param(1);
drA = param(2);
grB = param(3);
drB = param(4);
K = param(5);

dydt = zeros(2, 1);
dydt(1) = modelA(y, grA, drA);
dydt(2) = modelB(y, grB, drB, K);
dydt(3) = 0;
end

function dAdt = modelA(y, gr, dr)
A = y(1);
B = y(2);
input = y(3);

dAdt = gr * input * (1 - A) - dr * A * B;
end

function dBdt = modelB(y, gr, dr, K)
A = y(1);
B = y(2);
input = y(3);

if B < 1 % —‘z‚ÍB < 1 + K - (³‚Ì”÷¬’l)
    dBdt = gr * input * (1 - B) / (K + 1 - B) - dr * B;
else
    slope = (gr * input / K + dr);
    dBdt = - slope * (B - 1) - dr;
    if dBdt < -1
        dBdt = -1;
    end
end

end
