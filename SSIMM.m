load('main_Bayes201103-200019.mat');
% b c f n l t

ssimm_matrix = zeros(size(estimated_force_cell_matrix));

n = division_number;
m = n;

for c=1:length(cell_id_list)
    for f=1:length(force_id_list)
        [tf, ~] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
        Fx = tf(1:2:length(tf))';
        Fy = tf(2:2:length(tf))';
        Fx = reshape(Fx, n, m);
        Fy = reshape(Fy, n, m);
        tf_img = sqrt(reshape(Fx .^ 2 + Fy .^ 2, n, m));
        
        for b=1:length(bead_number_list)
            tmpssimm = zeros(size(ssimm_matrix(b, c, f, :, :, :)));
            for i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                tmpest = estimated_force_cell_matrix(b, c, f, :, :, :);
                ef = tmpest{i};
                efx = ef(1:2:length(ef));
                efy = ef(2:2:length(ef));
                ef_img = sqrt(reshape(efx .^ 2 + efy .^ 2, n, m));
                tmpssimm(i) = ssim(ef_img, tf_img);
            end
            ssimm_matrix(b, c, f, :, :, :) = tmpssimm;
        end
    end
end

%%

M = ssimm_matrix;

SM = zeros(length(bead_number_list), length(log_lam_list), length(noise_radius_list));
SUB = zeros(length(bead_number_list), length(log_lam_list), length(noise_radius_list));
SLB = zeros(length(bead_number_list), length(log_lam_list), length(noise_radius_list));

for b=1:length(bead_number_list)
    for l=1:length(log_lam_list)
        for r=1:length(noise_radius_list)
            ab = betafit(reshape(M(b, :, :, r, l, :), 1, []));
            SM(b, l, r) = ab(1) / (ab(1) + ab(2));
            SUB(b, l, r) = betaincinv(0.95, ab(1), ab(2));
            SLB(b, l, r) = betaincinv(0.05, ab(1), ab(2));
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=3:length(log_lam_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, ...
        SM(:, l, 1), ...
        SLB(:, l, 1) - SM(:, l, 1), ...
        SUB(:, l, 1) - SM(:, l, 1), ...
        'DisplayName', "log\lambda=" + log_lam_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);


%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=7:length(noise_radius_list)
    errorbar((bead_number_list + 5 * (r - 7 - 1)) / 255, ...
        SM(:, 4, r), ...
        SLB(:, 4, r) - SM(:, 4, r), ...
        SUB(:, 4, r) - SM(:, 4, r), ...
        'DisplayName', "log\lambda=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
load('mat/main_Ridge200718-210820.mat');
% b c f n l t

ssimm_matrix = zeros(size(estimated_force_cell_matrix));

n = division_number;
m = n;

for c=1:length(cell_id_list)
    for f=1:length(force_id_list)
        [tf, ~] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
        Fx = tf(1:2:length(tf))';
        Fy = tf(2:2:length(tf))';
        Fx = reshape(Fx, n, m);
        Fy = reshape(Fy, n, m);
        tf_img = sqrt(reshape(Fx .^ 2 + Fy .^ 2, n, m));
        
        for b=1:length(bead_number_list)
            tmpssimm = zeros(size(ssimm_matrix(b, c, f, :, :, :)));
            for i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                tmpest = estimated_force_cell_matrix(b, c, f, :, :, :);
                ef = tmpest{i};
                efx = ef(1:2:length(ef));
                efy = ef(2:2:length(ef));
                ef_img = sqrt(reshape(efx .^ 2 + efy .^ 2, n, m));
                tmpssimm(i) = ssim(ef_img, tf_img);
            end
            ssimm_matrix(b, c, f, :, :, :) = tmpssimm;
        end
    end
end

%%

M = ssimm_matrix;

SM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
SUB = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
SLB = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));

for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            ab = betafit(reshape(M(b, :, :, r, l, :), 1, []));
            SM(b, l, r) = ab(1) / (ab(1) + ab(2));
            SUB(b, l, r) = betaincinv(0.95, ab(1), ab(2));
            SLB(b, l, r) = betaincinv(0.05, ab(1), ab(2));
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=3:length(log_k_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, ...
        SM(:, l, 1), ...
        SLB(:, l, 1) - SM(:, l, 1), ...
        SUB(:, l, 1) - SM(:, l, 1), ...
        'DisplayName', "logk=" + log_k_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);


%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=7:length(noise_radius_list)
    errorbar((bead_number_list + 5 * (r - 7 - 1)) / 255, ...
        SM(:, 4, r), ...
        SLB(:, 4, r) - SM(:, 4, r), ...
        SUB(:, 4, r) - SM(:, 4, r), ...
        'DisplayName', "log\lambda=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
load('mat/main_Lasso200718-192141.mat');
% b c f n l t

ssimm_matrix = zeros(size(estimated_force_cell_matrix));

n = division_number;
m = n;

for c=1:length(cell_id_list)
    for f=1:length(force_id_list)
        [tf, ~] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
        Fx = tf(1:2:length(tf))';
        Fy = tf(2:2:length(tf))';
        Fx = reshape(Fx, n, m);
        Fy = reshape(Fy, n, m);
        tf_img = sqrt(reshape(Fx .^ 2 + Fy .^ 2, n, m));
        
        for b=1:length(bead_number_list)
            tmpssimm = zeros(size(ssimm_matrix(b, c, f, :, :, :)));
            for i=1:numel(estimated_force_cell_matrix(b, c, f, :, :, :))
                tmpest = estimated_force_cell_matrix(b, c, f, :, :, :);
                ef = tmpest{i};
                efx = ef(1:2:length(ef));
                efy = ef(2:2:length(ef));
                ef_img = sqrt(reshape(efx .^ 2 + efy .^ 2, n, m));
                tmpssimm(i) = ssim(ef_img, tf_img);
            end
            ssimm_matrix(b, c, f, :, :, :) = tmpssimm;
        end
    end
end

%%

M = ssimm_matrix;

SM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
SUB = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
SLB = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));

for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            ab = betafit(reshape(M(b, :, :, r, l, :), 1, []));
            SM(b, l, r) = ab(1) / (ab(1) + ab(2));
            SUB(b, l, r) = betaincinv(0.95, ab(1), ab(2));
            SLB(b, l, r) = betaincinv(0.05, ab(1), ab(2));
        end
    end
end

%% without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=3:length(log_k_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, ...
        SM(:, l, 1), ...
        SLB(:, l, 1) - SM(:, l, 1), ...
        SUB(:, l, 1) - SM(:, l, 1), ...
        'DisplayName', "logk=" + log_k_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);


%% with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for r=7:length(noise_radius_list)
    errorbar((bead_number_list + 5 * (r - 7 - 1)) / 255, ...
        SM(:, 4, r), ...
        SLB(:, 4, r) - SM(:, 4, r), ...
        SUB(:, 4, r) - SM(:, 4, r), ...
        'DisplayName', "log\lambda=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
clear;
load('BayesSSIMM_200719-184907.mat');
SMB = SM;
SLBB = SLB;
SUBB = SUB;
clear SM SLB SUB;
load('RidgeSSIMM200718-210820.mat');
SMR= SM;
SLBR = SLB;
SUBR = SUB;
clear SM SLB SUB;
load('LassoSSIMM200718-192141.mat');
SML = SM;
SLBL = SLB;
SUBL = SUB;
clear SM SLB SUB;

%%
figure; hold on;
ax = gca;
ax.FontSize = 15;
r = 7;
bead_number_list = [20:10:80, 100:100:800];
e = errorbar(bead_number_list / 255, ...
        SMB(:, 4, r), ...
        SLBB(:, 4, r) - SMB(:, 4, r), ...
        SUBB(:, 4, r) - SMB(:, 4, r), ...
        'DisplayName', "Bayes");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% e = errorbar((bead_number_list + 3) / 255, sqrt(EMR(:, 4, 1)), sqrt(ESR(:, 4, 1)), 'DisplayName', "Ridge");
e = errorbar(bead_number_list / 255, ...
        SMR(:, 4, r), ...
        SLBR(:, 4, r) - SMR(:, 4, r), ...
        SUBR(:, 4, r) - SMR(:, 4, r), ...
        'DisplayName', "Ridge");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% e.Marker = 'o';
% e.MarkerFaceColor = 'white';
% e.LineWidth = 1;
% e = errorbar((bead_number_list + 6) / 255, sqrt(EML(:, 4, 1)), sqrt(ESL(:, 4, 1)), 'DisplayName', "Lasso");
e = errorbar(bead_number_list / 255, ...
        SML(:, 4, r), ...
        SLBL(:, 4, r) - SML(:, 4, r), ...
        SUBL(:, 4, r) - SML(:, 4, r), ...
        'DisplayName', "Ridge");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% e.Marker = 'o';
% e.MarkerFaceColor = 'white';
% e.LineWidth = 1;
% title('Root Mean Squared Error');
xlabel('Bead density [-]');
xlim([0 3.3]);
ylim([0.5 1]);
ylabel('SSIM');
legend('Location', 'southeast');
pbaspect([1 1 1]);
% text(0.2, 4.8e-4, "Error Bar: SD (n = 30)", 'FontSize', 13);
%%
disp(datestr(datetime(), 'yymmdd-HHMMSS'));

