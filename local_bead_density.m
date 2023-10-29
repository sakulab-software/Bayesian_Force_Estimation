% B = 500;
% rng(51);
% IBLx = rand(B, 1);
% IBLy = rand(B, 1);
% bead_density_radius = 0.005;
% g = 0:0.01:1;
% 
% [X, Y] = meshgrid(g, g);
% X = reshape(X, [], 1);
% Y = reshape(Y, [], 1);
% 
% gms = zeros(numel(X), B);
% for b=1:B
%     gms(:, b) = mvnpdf([X, Y], [IBLx(b), IBLy(b)], eye(2) * bead_density_radius);
% end
% 
% gm = sum(gms, 2);
% gm = reshape(gm, length(g), length(g));
% X = reshape(X, length(g), length(g));
% Y = reshape(Y, length(g), length(g));
% 
% figure; hold on;
% p = pcolor(X, Y, gm);
% p.LineStyle = 'none';
% colorbar;
% scatter(IBLx, IBLy);
% xlim([min(g), max(g)]);
% ylim([min(g), max(g)]);
% pbaspect([1 1 1]);

%%
function gm = local_bead_density(B, IBLx, IBLy, bead_density_radius, division_number, X, Y)
X = reshape(X, [], 1);
Y = reshape(Y, [], 1);
gms = zeros(numel(X), B);
for b=1:B
    gms(:, b) = mvnpdf([X, Y], [IBLx(b), IBLy(b)], eye(2) * bead_density_radius);
end

gm = sum(gms, 2);
gm = reshape(gm, division_number, division_number);
end
% %%
% 
% [minAIC,numComponents] = min(AIC);
% gmPDF = @(x,y)reshape(pdf(gm{numComponents},[x(:) y(:)]),size(x));
% figure; hold on;
% scatter(IBLx, IBLy);
% h = fcontour(gmPDF, [0 1 0 1]);
% title('Simulated Data and Contour lines of pdf');
% 
% [minBIC,numComponents] = min(BIC);
% gmPDF = @(x,y)reshape(pdf(gm{numComponents},[x(:) y(:)]),size(x));
% figure; hold on;
% scatter(IBLx, IBLy);
% h = fcontour(gmPDF, [0 1 0 1]);
% title('Simulated Data and Contour lines of pdf');
% 
% %%
% 
% [minABIC,numComponents] = min(AIC + BIC);
% gmPDF = @(x,y)reshape(pdf(gm{numComponents},[x(:) y(:)]),size(x));
% figure; hold on;
% scatter(IBLx, IBLy);
% h = fcontour(gmPDF, [0 1 0 1]);
% title('Simulated Data and Contour lines of pdf');
% 
% %%
% 
% [X, Y] = meshgrid(0:0.01:1, 0:0.01:1);
% C = pdf(gm{numComponents}, [reshape(X, [], 1), reshape(Y, [], 1)]);
% C = reshape(C, size(X));
% surf(X, Y, C);
