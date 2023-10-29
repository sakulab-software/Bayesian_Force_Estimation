function score = dtm(EFMap, TFMap)
ind = TFMap > min(TFMap);
EFMap = EFMap(ind);
EFMap = EFMap(:);
TFMap = TFMap(ind);
TFMap = TFMap(:);
score = ((EFMap - TFMap) ./ TFMap) .^ 2;
score = mean(score);
end