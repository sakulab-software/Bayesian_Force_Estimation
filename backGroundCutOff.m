function [Fx, Fy, t] = backGroundCutOff(F, G, U)
n = size(G, 1) / 2;
rsn = round(sqrt(n));
t = fminbnd(@(t) e(cutOff(F, t), G, U), 0, 1);
F = cutOff(F, t);
Fx = reshape(F(1:2:end), rsn, []);
Fy = reshape(F(2:2:end), rsn, []);
end

function E = e(F, G, U)
E = sum((U - F * G) .^ 2);
end

function Z = cutOff(F, t)
Fx = F(1:2:end);
Fy = F(2:2:end);
M = sqrt(Fx .^ 2 + Fy .^ 2);
m = max(M, [], 'all');
ind = M <= t * m;
Fx(ind) = 0;
Fy(ind) = 0;
Z = zeros(size(F));
Z(1:2:end) = Fx;
Z(2:2:end) = Fy;
end