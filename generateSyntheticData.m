function [X, Y, Fx, Fy, IBLx, IBLy, BD, BDx, BDy, dBD, dBDx, dBDy, Sim] = ... 
    generateSyntheticData(divisionNumber, BeadNumber, w, h, forceScale, cellID, forceID, noiseRadius, randomState)

n = divisionNumber;
m = divisionNumber;
N = n * m;

X = linspace(0, w, n);
Y = linspace(0, h, m);
[X, Y] = meshgrid(X, Y);

if forceID == 1
    m = load(fullfile('cell_force_field', num2str(cellID), "cell" + num2str(cellID) + "_new1.mat"));
    X = m.X;
    Y = m.Y;
    Fx = m.U;
    Fy = m.V;
    mF = max(sqrt((Fx .^ 2) + (Fy .^ 2)), [], 'all');
    Fx = Fx / mF * forceScale;
    Fy = Fy / mF * forceScale;
else
    [tf, ~] = cell_force_field(cellID, forceID, forceScale);
    Fx = tf(1:2:length(tf))';
    Fy = tf(2:2:length(tf))';
    Fx = reshape(Fx, n, m);
    Fy = reshape(Fy, n, m);
end

rng(randomState);
IBLx = rand(BeadNumber, 1) * w;
IBLy = rand(BeadNumber, 1) * h;

Sim = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy);
Sim = Sim.simulate();
[BDx, BDy] = Sim.observe_without_noise();
BD = sqrt(BDx .^ 2 + BDy .^ 2);
[dBDx, dBDy] = Sim.observe(noiseRadius, randomState + 1);
dBD = sqrt(dBDx .^ 2 + dBDy .^ 2);

end