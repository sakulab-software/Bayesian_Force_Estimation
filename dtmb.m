function score = dtmb(EFMap, TFMap)
ind = TFMap == min(TFMap);
EFBackMap = EFMap(ind);
EFBackMap = EFBackMap(:);
TFBackMap = TFMap(ind);
TFBackMap = TFBackMap(:);
TFMap = TFMap(~ind);
TFMap = TFMap(:);
score = mean(EFBackMap - TFBackMap) / mean(TFMap);
end