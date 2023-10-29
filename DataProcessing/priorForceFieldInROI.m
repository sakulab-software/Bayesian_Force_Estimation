function [X, Y, Px, Py] = priorForceFieldInROI(divisionNumber, CellShapeImage, r, namedProp)
% priorForceFieldInROI: computes prior force field in ROI defined on
% CellShapeImage.
%
%   [X, Y, Px, Py] = priorForceFieldROI(divisionNumber, CellShapeImage, r)
% 
% Inputs:
% CellSpaceImage: 2D Logical Matrix (N x N)
% true -> inside or on the edge
% false -> outside
% r: Rectangle, given by drawRectangle()
%
% Outputs:
% X, Y: meshgrid, dividing r (ROI) into rectangles put on every
% divisionNumber-th of width and height.
% Px, Py: force field defined on all points indicated by X and Y.
%

arguments
    divisionNumber double
    CellShapeImage (:, :) logical
    r images.roi.Rectangle
    namedProp.Visualize = 'off'
    namedProp.Verbose = 'off'
end

mustBeMember(namedProp.Visualize, ["on", "off"]);
mustBeMember(namedProp.Verbose, ["on", "off"]);

[X, Y, Px, Py] = priorField(divisionNumber, r, CellShapeImage, namedProp);

end

%%

function [X, Y, Px, Py] = priorField(divisionNumber, r, CellShapeImage, namedProp)
[X, Y, ind, ~, ~] = meshXY(divisionNumber, r.Position, CellShapeImage);
g = position2grid(divisionNumber, r.Position);
data0 = signedDistance(X, Y, g, ind);
[~, data] = reinit(g, data0, 10, namedProp);
[~, data] = mcf(g, data, 1, 10, namedProp);
[Px, Py] = gradient(data);
M = sqrt(Px .^ 2 + Py .^ 2);
PFM = priorForceMagnitude(data0, length(find(CellShapeImage)));
Px = Px .* (-PFM) ./ M;
Px(data0 > 0) = 0;
Py = Py .* (-PFM) ./ M;
Py(data0 > 0) = 0;

if strcmp(namedProp.Visualize, "on")
    visualizeGrad(X, Y, Px, Py, CellShapeImage);
end

    function PFM = priorForceMagnitude(signedDistance, r)
        PFM = exp(- signedDistance .^ 2 / r);
    end
end

function h = visualizeGrad(X, Y, Gx, Gy, CellShapeImage)
h = figure; hold on;
% pcolor(g.xs{1}', g.xs{2}', Hxy);
imagesc(CellShapeImage);
quiver(X, Y, Gx, Gy);
end

function [X, Y, ind, w, h] = meshXY(divisionNumber, position, CellShapeImage)

[x, y, w, h] = position2xywh(position);

if w ~= h
    error("width and height must be same value pair.");
end

X = linspace(x, x + w, divisionNumber);
Y = linspace(y, y + h, divisionNumber);
[X, Y] = meshgrid(X, Y);
ind = false(size(X));
for i=1:numel(X)
    ind(i) = CellShapeImage(round(Y(i)), round(X(i)));
end

end

function g = position2grid(divisionNumber, position)

[x, y, w, h] = position2xywh(position);

g.dim = 2;
g.min = [y; x];
g.max = [y + h; x + w];
g.dx = [1; 1] .* (g.max - g.min) / (divisionNumber - 1);
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

end

function data = signedDistance(X, Y, g, ind)

data = zeros(size(X));
data(ind) = -1;
data(~ind) = 1;

set(0, 'DefaultFigureVisible', 'off');
[ ~, hdl ] = contour(g.xs{1}, g.xs{2}, data, [0, 0], 'b');
contourMatrix = hdl.ContourMatrix;
delete(gcf);
set(0, 'DefaultFigureVisible', 'on');

for i=1:numel(X)
    M = min(sqrt((contourMatrix(1, :) - X(i)) .^ 2 + ...
                 (contourMatrix(2, :) - Y(i)) .^ 2));
    data(i) = data(i) * M;
end

end

function [x, y, w, h] = position2xywh(position)

position = round(position);
x = position(1);
y = position(2);
w = position(3);
h = position(4);

end

function [t, y] = reinit(g, data0, duration, namedProp)

schemeFunc = @termReinit;
schemeData.grid = g;
schemeData.derivFunc = @upwindFirstENO3;
schemeData.initial = data0;
schemeData.subcell_fix_order = 1;

integratorOptions = odeCFLset('factorCFL', 0.9, 'stats', namedProp.Verbose);
integratorFunc = @odeCFL3;

tSpan = [0, duration];
[t, y] = odeCFL3(schemeFunc, tSpan, data0(:), integratorOptions, schemeData);
y = reshape(y, size(g.xs{1}));

end

function [t, y] = mcf(g, data0, b, duration, namedProp)

schemeFunc = @termCurvature;
schemeData.grid = g;
schemeData.curvatureFunc = @curvatureSecond;
schemeData.b = b;

integratorOptions = odeCFLset('factorCFL', 0.9, 'stats', namedProp.Verbose);
integratorFunc = @odeCFL3;

tSpan = [0, duration];
[t, y] = odeCFL3(schemeFunc, tSpan, data0(:), integratorOptions, schemeData);
y = reshape(y, size(g.xs{1}));
end

