function MTMGridComp(thr)

mat = load('mainBayes2210731-050728.mat');
% b c f n l thr t

% n = divisionNumber;
% m = n;
% X = linspace(0, w, n);
% Y = linspace(0, h, m);
% [X, Y] = meshgrid(X, Y);


MTMMatrixBayes = zeros(size(mat.estimated_force_cell_matrix(:, :, :, :, :, thr, :)));
% AUCMatrixBayes = zeros(size(estimated_force_cell_matrix));
divisionNumber = mat.divisionNumber;
beadNumberList = mat.beadNumberList;
cellIDList = mat.cellIDList;
forceIDList = mat.forceIDList;
forceScale = mat.forceScale;
n = divisionNumber;
m = n;

for c=1:length(cellIDList)
    for f=1:length(forceIDList)
        for b=1:length(beadNumberList)
            tmpMTM = zeros(size(MTMMatrixBayes(b, c, f, :, :, :)));
            tmpEst = mat.estimated_force_cell_matrix(b, c, f, :, :, thr, :);
%             tmpIBL = mat.initial_bead_location_cell_matrix(b, c, f, :, :, thr, :);
%             tmpBD = mat.bead_displacement_cell_matrix(b, c, f, :, :, thr, :);
            
%             parfor i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
            ci = cellIDList(c);
            fi = forceIDList(f);
            parfor i=1:numel(mat.estimated_force_cell_matrix(b, c, f, :, :, thr, :))
                EF = tmpEst{i};
                EFx = EF(1:2:end);
                EFy = EF(2:2:end);
                EFx = reshape(EFx, n, m);
                EFy = reshape(EFy, n, m);
%                 IBL = tmpIBL{i};
%                 IBLx = IBL(1:2:end);
%                 IBLy = IBL(2:2:end);
%                 Sim = TFM_computation(X, Y, EFx, EFy, IBLx(:), IBLy(:));
%                 Sim = Sim.simulate();
%                 BD = tmpBD{i};
%                 [EFx, EFy, ~] = backGroundCutOff(EF, Sim.G', BD);
                EF = reshape([reshape(EFx, 1, []); reshape(EFy, 1, [])], [], 1);
                [thrs, mses] = tfm_multi_threshold_mse(EF, divisionNumber, ci, fi, forceScale);
                tmpMTM(i) = tfmauc(thrs, mses);
            end
            MTMMatrixBayes(b, c, f, :, :, :) = tmpMTM;
        end
    end
end

fn = "MTMGridComp-" + thr + "-" + datestr(datetime(), 'yymmdd-HHMMSS') + ".mat";
save(fn, 'MTMMatrixBayes', '-v7.3');

end