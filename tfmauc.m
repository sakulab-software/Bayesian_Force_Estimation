function auc = tfmauc(tpfs, fpfs)
auc = trapz(fpfs, tpfs);
end