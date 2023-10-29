
beadNumberList = [20:10:80, 100:100:800];
beadNumberList = beadNumberList([3, 8, 9:12]);
% beadNumberList = [40, 80, 100, 300, 500];
beadNumberList = [45, 90, 180, 270, 360, 450];
trial_number = 2;
trial_number = 10;
lamList = logspace(-6, 1, 8);
lamList = 0.1;
thrList = 1e-1 * [0.01 0.05 0.1 0.5 1 5 10];
thrList = 0.1;
cellIDList = 1:5;
forceIDList = 1:5;
% noiseRadiusList = logspace(-8, 0, 9);
noiseRadiusList = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1];
noiseRadiusList = 0.01;
divisionNumber = 30;
forceScale = 1e-1;
forceScale = 5e-1;
randomStateList = 1:length(beadNumberList) * trial_number;
randomStateList = reshape(randomStateList, trial_number, length(beadNumberList));

w = 15;
h = w;

shp1 = [
    length(beadNumberList), ...
    length(cellIDList), ...
    length(forceIDList), ...
    length(noiseRadiusList), ...
    length(lamList), ...
    length(thrList), ...
    trial_number
    ];

shp2 = shp1(2:end);
shp3 = shp2(2:end);
shp4 = shp3(2:end);
shp5 = shp4(2:end);
shp6 = shp5(2:end);

initial_bead_location_cell_matrix = cell(shp1);
bead_displacement_cell_matrix = cell(shp1);
estimated_force_cell_matrix = cell(shp1);
estimated_force_error_cell_matrix = cell(shp1);

parfor b=1:length(beadNumberList)
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
                
                for lam=1:length(lamList)
                    idata5 = cell(shp6);
                    bdata5 = cell(shp6);
                    efdata5 = cell(shp6);
                    eedata5 = cell(shp6);
                    
                    for thr=1:length(thrList)
                        for t=1:trial_number
                            [idata5{thr, t}, bdata5{thr, t}, ...
                                efdata5{thr, t}, eedata5{thr, t}] =  ...
                                TFMRoutineWithBayesEstimation(...
                                    divisionNumber, beadNumberList(b), w, h, ...
                                    forceScale, cellIDList(c), forceIDList(f), ...
                                    1, lamList(lam), thrList(thr), noiseRadiusList(r), randomStateList(t, b));
                        end
                    end
                    
                    idata4(lam, :, :) = idata5;
                    bdata4(lam, :, :) = bdata5;
                    efdata4(lam, :, :) = efdata5;
                    eedata4(lam, :, :) = eedata5;
                end

                idata3(r, :, :, :) = idata4;
                bdata3(r, :, :, :) = bdata4;
                efdata3(r, :, :, :) = efdata4;
                eedata3(r, :, :, :) = eedata4;
            end
            idata2(f, :, :, :, :) = idata3;
            bdata2(f, :, :, :, :) = bdata3;
            efdata2(f, :, :, :, :) = efdata3;
            eedata2(f, :, :, :, :) = eedata3;
        end
        idata1(c, :, :, :, :, :) = idata2;
        bdata1(c, :, :, :, :, :) = bdata2;
        efdata1(c, :, :, :, :, :) = efdata2;
        eedata1(c, :, :, :, :, :) = eedata2;
    end
    initial_bead_location_cell_matrix(b, :, :, :, :, :, :) = idata1;
    bead_displacement_cell_matrix(b, : ,:, :, :, :, :) = bdata1;
    estimated_force_cell_matrix(b, :, :, :, :, :, :) = efdata1;
    estimated_force_error_cell_matrix(b, :, :, :, :, :, :) = eedata1;
end

fn = "mainBayes2-" + datestr(datetime(), 'yymmdd-HHMMSS') + ".mat";
save(fn, '-v7.3');
