ef = reshape(ef, 1, []);
[tpfs, fpfs, ths] = tfmroc(ef, division_number, cell_number, force_number, force_scale);
[~, ~, th] = backGroundCutOff(ef, tfmc.G', reshape([dBDx'; dBDy'], 1, []));
[~, I] = min(abs(ths - th));

figure; hold on;
ax = gca;
ax.FontSize = 16;
plot(fpfs, tpfs, 'Marker', 'o');
scatter(fpfs(I), tpfs(I), 'MarkerFaceColor', 'k');
title('ROC');
xlabel('FPF');
xlim([0, 1]);
ylabel('TPF');
ylim([0, 1]);
% legend('Location', 'eastoutside');
pbaspect([1 1 1]);
display(tfmauc(tpfs, fpfs))