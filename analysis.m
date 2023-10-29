load('mat/out_new_code_bayes.mat');

mse_matrix = zeros(size(force_err_criterias_matrix));
for i=1:numel(mse_matrix)
    mse_matrix(i) = force_err_criterias_matrix{i}.MSE;
end

mean_mse = zeros(size(bead_number_list));
std_mse = zeros(size(bead_number_list));

for i=1:length(bead_number_list)
    mean_mse(i) = mean(mse_matrix(:, :, :, :, i), 'all');
    std_mse(i) = std(reshape(mse_matrix(:, :, :, :, i), 1, []));
end

errorbar(bead_number_list, mean_mse, std_mse);

data_violin = zeros(...
    numel(force_err_criterias_matrix(:, :, :, :, 1)), ...
    length(bead_number_list));

for i=1:length(bead_number_list)
    data_violin(:, i) = reshape(mse_matrix(:, :, :, :, i), [], 1);
end

violin(data_violin);

mean_mse_matrix = zeros(length(log_k_list), length(bead_number_list));
std_mse_matrix = zeros(length(log_k_list), length(bead_number_list));

for i=1:length(log_k_list)
    for j=1:length(bead_number_list)
        mean_mse_matrix(i, j) = mean(mse_matrix(:, i, :, :, j), 'all');
        std_mse_matrix(i, j) = std(reshape(mse_matrix(:, i, :, :, j), 1, []));
    end
end

figure, heatmap(mean_mse_matrix, 'Title', 'AVERAGE');
figure, heatmap(std_mse_matrix, 'Title', 'STD');

%%
clear;
load('main_Bayes200929-220023.mat');
% b c f r l t

log_lambda_list = log_lam_list;
M = cell2mat(estimated_force_error_cell_matrix);

EM = zeros(length(bead_number_list), length(log_lambda_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(log_lambda_list), length(noise_radius_list));
for b=1:length(bead_number_list)
    for l=1:length(log_lambda_list)
        for r=1:length(noise_radius_list)
            EM(b, l, r) = mean(M(b, :, :, r, l, :), 'all');
            ES(b, l, r) = std(reshape(M(b, :, :, r, l, :), 1, []));
        end
    end
end
% 
% M = estimated_force_cell_matrix;
% 
% EM = zeros(length(bead_number_list), length(log_lambda_list), length(noise_radius_list));
% ES = zeros(length(bead_number_list), length(log_lambda_list), length(noise_radius_list));
% 
% for b=1:length(bead_number_list)
%     for l=1:length(log_lambda_list)
%         for r=1:length(noise_radius_list)
%             tmp = zeros(length(cell_id_list), length(force_id_list));
%             for c=1:length(cell_id_list)
%                 for f=1:length(force_id_list)
%                     [F, I] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
%                     tmp2 = estimated_force_cell_matrix{b, c, f, r, l, 1};
%                     F = sqrt(F(1:2:length(F)) .^ 2 + F(2:2:length(F)) .^ 2);
%                     tmp2 = sqrt(tmp2(1:2:length(tmp2)) .^ 2 + tmp2(2:2:length(tmp2)) .^ 2);
%                     tmp2 = abs((tmp2(I) - F(I)) / F(I));
%                     tmp(c, f) = mean(tmp2, 'all');
%                 end
%             end
%             EM(b, l, r) = mean(tmp, 'all');
%             ES(b, l, r) = std(tmp, 0, 'all');
%         end
%     end
% end

%% analysis in the conditions without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=1:length(log_lambda_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, sqrt(EM(:, l, 1)), sqrt(ES(:, l, 1)), 'DisplayName', "log\lambda=" + log_lambda_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%% analysis in the conditions with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
% ax.YScale = 'log';
for r=1:length(noise_radius_list) - 2
    errorbar((bead_number_list + 5 * (r - 1)) / 255, sqrt(EM(:, 4, r)), sqrt(ES(:, 4, r)), 'DisplayName', "noise_radius=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%% 
clear;
load("out_EN.mat");
% t, lk, a, f, c, b

M = zeros(size(force_err_criterias_matrix));
for i=1:numel(force_err_criterias_matrix)
    M(i) = force_err_criterias_matrix{i}.MSE;
end

EM = zeros(length(bead_number_list), length(log_k_list), length(a_list));
ES = zeros(length(bead_number_list), length(log_k_list), length(a_list));

for b=1:length(bead_number_list)
    for lk=1:length(log_k_list)
        for a=1:length(a_list)
            EM(b, lk, a) = mean(M(:, lk, a, :, :, b), 'all');
            ES(b, lk, a) = std(reshape(M(:, lk, a, :, :, b), 1, []));
        end
    end
end

%%

for b=1:length(bead_number_list)
    figure;
    heatmap(a_list, log_k_list, reshape(sqrt(EM(b, :, :)), length(log_k_list), length(a_list)));
    colormap viridis(100);
end

%%
clear;
load('mat/main_Ridge200718-210820.mat');

% b c f r l t

M = cell2mat(estimated_force_error_cell_matrix);

EM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            EM(b, l, r) = mean(M(b, :, :, r, l, :), 'all');
            ES(b, l, r) = std(reshape(M(b, :, :, r, l, :), 1, []));
        end
    end
end

M = estimated_force_cell_matrix;

EM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));

for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            tmp = zeros(length(cell_id_list), length(force_id_list));
            for c=1:length(cell_id_list)
                for f=1:length(force_id_list)
                    [F, I] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
                    tmp2 = estimated_force_cell_matrix{b, c, f, r, l, 1};
                    F = sqrt(F(1:2:length(F)) .^ 2 + F(2:2:length(F)) .^ 2);
                    tmp2 = sqrt(tmp2(1:2:length(tmp2)) .^ 2 + tmp2(2:2:length(tmp2)) .^ 2);
                    tmp2 = abs((tmp2(I) - F(I)) / F(I));
                    tmp(c, f) = mean(tmp2, 'all');
                end
            end
            EM(b, l, r) = mean(tmp, 'all');
            ES(b, l, r) = std(tmp, 0, 'all');
        end
    end
end

%% analysis in the conditions without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=1:length(log_k_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, sqrt(EM(:, l, 1)), sqrt(ES(:, l, 1)), 'DisplayName', "log\lambda=" + log_k_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%% analysis in the conditions with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
ax.YScale = 'log';
for r=1:length(noise_radius_list) - 2
    errorbar((bead_number_list + 5 * (r - 1)) / 255, sqrt(EM(:, 4, r)), sqrt(ES(:, 4, r)), 'DisplayName', "noise_radius=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
clear;
load('mat/main_Lasso200718-192141.mat');

% b c f r l t

M = cell2mat(estimated_force_error_cell_matrix);

EM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            EM(b, l, r) = mean(M(b, :, :, r, l, :), 'all');
            ES(b, l, r) = std(reshape(M(b, :, :, r, l, :), 1, []));
        end
    end
end

M = estimated_force_cell_matrix;

EM = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(log_k_list), length(noise_radius_list));

for b=1:length(bead_number_list)
    for l=1:length(log_k_list)
        for r=1:length(noise_radius_list)
            tmp = zeros(length(cell_id_list), length(force_id_list));
            for c=1:length(cell_id_list)
                for f=1:length(force_id_list)
                    [F, I] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
                    tmp2 = estimated_force_cell_matrix{b, c, f, r, l, 1};
                    F = sqrt(F(1:2:length(F)) .^ 2 + F(2:2:length(F)) .^ 2);
                    tmp2 = sqrt(tmp2(1:2:length(tmp2)) .^ 2 + tmp2(2:2:length(tmp2)) .^ 2);
                    tmp2 = abs((tmp2(I) - F(I)) / F(I));
                    tmp(c, f) = mean(tmp2, 'all');
                end
            end
            EM(b, l, r) = mean(tmp, 'all');
            ES(b, l, r) = std(tmp, 0, 'all');
        end
    end
end

%% analysis in the conditions without noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
for l=1:length(log_k_list)
    errorbar((bead_number_list + 5 * (l - 1)) / 255, sqrt(EM(:, l, 1)), sqrt(ES(:, l, 1)), 'DisplayName', "log\lambda=" + log_k_list(l));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%% analysis in the conditions with noise

figure; hold on;
ax = gca;
ax.FontSize = 16;
ax.YScale = 'log';
for r=1:length(noise_radius_list) - 2
    errorbar((bead_number_list + 5 * (r - 1)) / 255, sqrt(EM(:, 4, r)), sqrt(ES(:, 4, r)), 'DisplayName', "noise_radius=" + noise_radius_list(r));
end
title('Root Mean Squared Error');
xlabel('Bead density [-]');
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);

%%
clear;
load('mat/main_Bayes200718-184731.mat');
load E;

%%
figure; hold on;
ax = gca;
ax.FontSize = 15;
bead_number_list = [20:10:80, 100:100:800];
e = errorbar(bead_number_list / 255, sqrt(EMB(:, 4, 1)), sqrt(ESB(:, 4, 1)), 'DisplayName', "Bayes");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 3) / 255, sqrt(EMR(:, 4, 1)), sqrt(ESR(:, 4, 1)), 'DisplayName', "Ridge");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 6) / 255, sqrt(EML(:, 4, 1)), sqrt(ESL(:, 4, 1)), 'DisplayName', "Lasso");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% title('Root Mean Squared Error');
xlabel('Bead density [-]');
xlim([0 3.3]);
% ylim([0 40e-4]);
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);
% text(1.7, 3.5e-3, "Error Bar: SD (n = 30)", 'FontSize', 13);

%%
clear;
load('main_Ridge200929-164743.mat');
% b c f n t

M = cell2mat(estimated_force_error_cell_matrix);

EM = zeros(length(bead_number_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(noise_radius_list));
for b=1:length(bead_number_list)
    for r=1:length(noise_radius_list)
        EM(b, r) = mean(M(b, :, :, r, :), 'all');
        ES(b, r) = std(reshape(M(b, :, :, r, :), 1, []));
    end
end
% 
% M = estimated_force_cell_matrix;
% 
% EM = zeros(length(bead_number_list), length(noise_radius_list));
% ES = zeros(length(bead_number_list), length(noise_radius_list));
% 
% for b=1:length(bead_number_list)
%     for r=1:length(noise_radius_list)
%         tmp = zeros(length(cell_id_list), length(force_id_list));
%         for c=1:length(cell_id_list)
%             for f=1:length(force_id_list)
%                 [F, I] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
%                 tmp2 = estimated_force_cell_matrix{b, c, f, r, 1};
%                 F = sqrt(F(1:2:length(F)) .^ 2 + F(2:2:length(F)) .^ 2);
%                 tmp2 = sqrt(tmp2(1:2:length(tmp2)) .^ 2 + tmp2(2:2:length(tmp2)) .^ 2);
%                 tmp2 = abs((tmp2(I) - F(I)) / F(I));
%                 tmp(c, f) = mean(tmp2, 'all');
%             end
%         end
%         EM(b, r) = mean(tmp, 'all');
%         ES(b, r) = std(tmp, 0, 'all');
%     end
% end

%%
clear;
load('main_Ridge200930-163904.mat'); %Lasso
% b c f n t

M = cell2mat(estimated_force_error_cell_matrix);

EM = zeros(length(bead_number_list), length(noise_radius_list));
ES = zeros(length(bead_number_list), length(noise_radius_list));
for b=1:length(bead_number_list)
    for r=1:length(noise_radius_list)
        EM(b, r) = mean(M(b, :, :, r, :), 'all');
        ES(b, r) = std(reshape(M(b, :, :, r, :), 1, []));
    end
end

% M = estimated_force_cell_matrix;
% 
% EM = zeros(length(bead_number_list), length(noise_radius_list));
% ES = zeros(length(bead_number_list), length(noise_radius_list));
% 
% for b=1:length(bead_number_list)
%     for r=1:length(noise_radius_list)
%         tmp = zeros(length(cell_id_list), length(force_id_list));
%         for c=1:length(cell_id_list)
%             for f=1:length(force_id_list)
%                 [F, I] = cell_force_field(cell_id_list(c), force_id_list(f), force_scale);
%                 tmp2 = estimated_force_cell_matrix{b, c, f, r, 1};
%                 F = sqrt(F(1:2:length(F)) .^ 2 + F(2:2:length(F)) .^ 2);
%                 tmp2 = sqrt(tmp2(1:2:length(tmp2)) .^ 2 + tmp2(2:2:length(tmp2)) .^ 2);
%                 tmp2 = abs((tmp2(I) - F(I)) / F(I));
%                 tmp(c, f) = mean(tmp2, 'all');
%             end
%         end
%         EM(b, r) = mean(tmp, 'all');
%         ES(b, r) = std(tmp, 0, 'all');
%     end
% end

%%
load BayesER200929-220023.mat
load RidgeER200929-164743.mat
load LassoER200930-163904.mat

%%
figure; hold on;
ax = gca;
ax.FontSize = 15;
bead_number_list = [20:10:80, 100:100:800];
e = errorbar(bead_number_list / 255, EMB(:, 4, 1), ESB(:, 4, 1), 'DisplayName', "Bayes");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 3) / 255, EMR(:, 4, 1), ESR(:, 4, 1), 'DisplayName', "Ridge");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 6) / 255, EML(:, 4, 1), ESL(:, 4, 1), 'DisplayName', "Lasso");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% title('Root Mean Squared Error');
xlabel('Bead density [-]');
xlim([0 3.3]);
% ylim([0 40e-4]);
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);
% text(1.7, 3.5e-3, "Error Bar: SD (n = 30)", 'FontSize', 13);

%%
load BayesMSE200929-220023.mat
load RidgeMSE200929-164743.mat
load LassoMSE200930-163904.mat

%%
figure; hold on;
ax = gca;
ax.FontSize = 15;
bead_number_list = [20:10:80, 100:100:800];
e = errorbar(bead_number_list / 255, sqrt(EMB(:, 4, 1)), sqrt(ESB(:, 4, 1)), 'DisplayName', "Bayes");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 3) / 255, sqrt(EMR(:, 4, 1)), sqrt(ESR(:, 4, 1)), 'DisplayName', "Ridge");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
e = errorbar((bead_number_list + 6) / 255, sqrt(EML(:, 4, 1)), sqrt(ESL(:, 4, 1)), 'DisplayName', "Lasso");
e.Marker = 'o';
e.MarkerFaceColor = 'white';
e.LineWidth = 1;
% title('Root Mean Squared Error');
xlabel('Bead density [-]');
xlim([0 3.3]);
% ylim([0 40e-4]);
ylabel('RMSE [pN]');
legend('Location', 'eastoutside');
pbaspect([1 1 1]);
% text(1.7, 3.5e-3, "Error Bar: SD (n = 30)", 'FontSize', 13);