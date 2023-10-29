
beadNumberList = [20:10:80, 100:100:800];
trialNumber = 5;
cellIDList = 1:5;
forceIDList = 1:5;
noiseRadiusList = logspace(-8, 0, 9);
divisionNumber = 30;
forceScale = 1e-1;
randomStateList = 1:length(beadNumberList) * trialNumber;
randomStateList = reshape(randomStateList, trialNumber, length(beadNumberList));

dA = 1e-1;
AList = dA:dA:1 - dA;

w = 15;
h = w;

shp1 = [
    length(beadNumberList), ...
    length(cellIDList), ...
    length(forceIDList), ...
    length(noiseRadiusList), ...
    length(AList), ...
    trialNumber
    ];

shp2 = shp1(2:end);
shp3 = shp2(2:end);
shp4 = shp3(2:end);
shp5 = shp4(2:end);

initial_bead_location_cell_matrix = cell(shp1);
bead_displacement_cell_matrix = cell(shp1);
estimated_force_cell_matrix = cell(shp1);
estimated_force_error_cell_matrix = cell(shp1);

for b=1:length(beadNumberList)
    idata1 = cell(shp2);
    bdata1 = cell(shp2);
    efdata1 = cell(shp2);
    eedata1 = cell(shp2);
    
    for c=1:length(cellIDList)
        idata2 = cell(shp3);
        bdata2 = cell(shp3);
        efdata2 = cell(shp3);
        eedata2 = cell(shp3);
        
        for f=1:length(forceIDList)
            idata3 = cell(shp4);
            bdata3 = cell(shp4);
            efdata3 = cell(shp4);
            eedata3 = cell(shp4);
            
            for r=1:length(noiseRadiusList)
                idata4 = cell(shp5);
                bdata4 = cell(shp5);
                efdata4 = cell(shp5);
                eedata4 = cell(shp5);
                
                for t=1:trialNumber
                    for a=1:length(AList)
                        [idata4{a, t}, bdata4{a, t}, ...
                            efdata4{a, t}, eedata4{a, t}] =  ...
                                TFMRoutineWithElasticNet(...
                                    divisionNumber, beadNumberList(b), w, h, forceScale, ...
                                    cellIDList(c), forceIDList(f), AList(a), noiseRadiusList(r), randomStateList(t, b));
                    end
                end
                idata3(r, :) = idata4;
                bdata3(r, :) = bdata4;
                efdata3(r, :) = efdata4;
                eedata3(r, :) = eedata4;
            end
            idata2(f, :, :) = idata3;
            bdata2(f, :, :) = bdata3;
            efdata2(f, :, :) = efdata3;
            eedata2(f, :, :) = eedata3;
        end
        idata1(c, :, :, :) = idata2;
        bdata1(c, :, :, :) = bdata2;
        efdata1(c, :, :, :) = efdata2;
        eedata1(c, :, :, :) = eedata2;
    end
    initial_bead_location_cell_matrix(b, :, :, :, :) = idata1;
    bead_displacement_cell_matrix(b, : ,:, :, :) = bdata1;
    estimated_force_cell_matrix(b, :, :, :, :) = efdata1;
    estimated_force_error_cell_matrix(b, :, :, :, :) = eedata1;
end

fn = "mainElasticNet" + datestr(datetime(), 'yymmdd-HHMMSS') + ".mat";
save(fn);

% 
% parfor b=1:length(bead_number_list)
%     idata1 = cell(shp2);
%     bdata1 = cell(shp2);
%     efdata1 = cell(shp2);
%     eedata1 = cell(shp2);
%     
%     for c=1:length(cell_id_list)
%         idata2 = cell(shp3);
%         bdata2 = cell(shp3);
%         efdata2 = cell(shp3);
%         eedata2 = cell(shp3);
%         
%         for f=1:length(force_id_list)
%             idata3 = cell(shp4);
%             bdata3 = cell(shp4);
%             efdata3 = cell(shp4);
%             eedata3 = cell(shp4);
%             
%             for r=1:length(noise_radius_list)
%                 idata4 = cell(shp5);
%                 bdata4 = cell(shp5);
%                 efdata4 = cell(shp5);
%                 eedata4 = cell(shp5);
%                 
%                 for logk=1:length(log_k_list)
%                     
%                     for t=1:trial_number
%                         
%                         n = division_number;
%                         m = division_number;
%                         N = n * m;
% 
%                         random_seed = random_seed_list(t, b);
%                         B = bead_number_list(b);
% 
%                         X = linspace(0, 1, n);
%                         Y = linspace(0, 1, m);
%                         [X, Y] = meshgrid(X, Y);
% 
%                         [tf, ~] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
%                         Fx = tf(1:2:length(tf))';
%                         Fy = tf(2:2:length(tf))';
%                         Fx = reshape(Fx, n, m);
%                         Fy = reshape(Fy, n, m);
% 
%                         rng(random_seed);
%                         IBLx = rand(B, 1);
%                         IBLy = rand(B, 1);
%                         IBL = reshape([IBLx, IBLy]', 1, [])';
%                         
%                         idata4{logk, t} = IBL;
% 
%                         tfmc = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy);
%                         tfmc = tfmc.simulate();
%                         [dBDx, dBDy] = tfmc.observe(noise_radius_list(r), random_seed);
%                         
%                         bdata4{logk, t} = reshape([dBDx, dBDy]', 1, [])';
%                         
%                         [dBDx2, dBDy2] = cut_noise_from_BD(dBDx, dBDy);
%                         BD = reshape([dBDx2, dBDy2]', 1, [])';
%                         
%                         ef = lasso(tfmc.G, BD, 'Lambda', 10 ^ log_k_list(logk));
%                         
%                         efdata4{logk, t} = ef;
%                         eedata4{logk, t} = mean((tf - ef) .^ 2, 'all') * 2;
%                     end
%                 end
% 
%                 idata3(r, :, :) = idata4;
%                 bdata3(r, :, :) = bdata4;
%                 efdata3(r, :, :) = efdata4;
%                 eedata3(r, :, :) = eedata4;
%             end
%             idata2(f, :, :, :) = idata3;
%             bdata2(f, :, :, :) = bdata3;
%             efdata2(f, :, :, :) = efdata3;
%             eedata2(f, :, :, :) = eedata3;
%         end
%         idata1(c, :, :, :, :) = idata2;
%         bdata1(c, :, :, :, :) = bdata2;
%         efdata1(c, :, :, :, :) = efdata2;
%         eedata1(c, :, :, :, :) = eedata2;
%     end
%     initial_bead_location_cell_matrix(b, :, :, :, :, :) = idata1;
%     bead_displacement_cell_matrix(b, : ,:, :, :, :) = bdata1;
%     estimated_force_cell_matrix(b, :, :, :, :, :) = efdata1;
%     estimated_force_error_cell_matrix(b, :, :, :, :, :) = eedata1;
% end
% 
% fn = "main_Lasso" + datestr(datetime(), 'yymmdd-HHMMSS') + ".mat";
% save(fn);