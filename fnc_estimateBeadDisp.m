function [eBDx, eBDy] = fnc_estimateBeadDisp(BX, BY, BDx, BDy)

feed = [BX(:), BY(:)];

modelX = fnc_trainSVM(feed, BDx(:));
modelY = fnc_trainSVM(feed, BDy(:));

eBDx = modelX.predictFcn(feed);
eBDy = modelY.predictFcn(feed);

end