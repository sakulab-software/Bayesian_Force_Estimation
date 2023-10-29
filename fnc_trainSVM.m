function [trainedModel] = fnc_trainSVM(trainingData, responseData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% 学習済み回帰モデルとその RMSE を返します。このコードは回帰学習器アプリで学習させたモデル
% を再作成します。生成されるコードを使用して、同じモデルでの新規データを使用した学習の自動化
% や、プログラミングによってモデルに学習させる方法の調査を行います。
%
%  入力:
%      trainingData: アプリにインポートされた行列と同じ列数とデータ型をもつ行列。
%
%      responseData: アプリにインポートされたベクトルと同じデータ型をもつベクト
%       ル。responseData の長さと trainingData の行数は等しくなければなりません。
%
%  出力:
%      trainedModel: 学習済みの回帰モデルを含む struct。この struct には、学習済みの
%       モデルに関する情報をもつさまざまなフィールドが含まれています。
%
%      trainedModel.predictFcn: 新規データに関する予測を行う関数。
%
%      validationRMSE: RMSE を表す double。アプリでは [モデル] ペインにモデルごとの
%       RMSE が表示されます。
%
% このコードを使用して新規データでモデルに学習させます。モデルに再学習させるには、元のデータ
% または新規データを入力引数 trainingData および responseData として指定し、コマンド
% ラインから関数を呼び出します。
%
% たとえば、元のデータセット T と応答 Y で学習させた回帰モデルに再学習させるには、次を入力
% します:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% 返された 'trainedModel' を使用して新規データ T2 の予測を行うには、次を使用します
%   yfit = trainedModel.predictFcn(T2)
%
% T2 は、学習に使用した予測列のみを含む行列でなければなりません。詳細については、次のように
% 入力してください:
%   trainedModel.HowToPredict

% MATLAB からの自動生成日: 2022/11/27 22:31:53


% 予測子と応答の抽出
% このコードは、データを処理して、モデルに学習させるのに適した
% 形状にします。
% 入力をテーブルへ変換
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2'});

predictorNames = {'column_1', 'column_2'};
predictors = inputTable(:, predictorNames);
response = responseData;
isCategoricalPredictor = [false, false];

% 回帰モデルの学習
% このコードは、すべてのモデル オプションを指定してモデルに学習させます。
responseScale = iqr(response);
if ~isfinite(responseScale) || responseScale == 0.0
    responseScale = 1.0;
end
boxConstraint = responseScale/1.349;
epsilon = responseScale/13.49;
regressionSVM = fitrsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 1.4, ...
    'BoxConstraint', boxConstraint, ...
    'Epsilon', epsilon, ...
    'Standardize', true);

% 関数 predict で結果の構造体を作成
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(regressionSVM, x);
trainedModel.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% 結果の構造体にさらにフィールドを追加
trainedModel.RegressionSVM = regressionSVM;
trainedModel.About = 'この構造体は、回帰学習器 R2022a からエクスポートされた学習済みのモデルです。';
trainedModel.HowToPredict = sprintf('新しい予測子列行列 X についての予測を行うには、次を使用します: \n yfit = c.predictFcn(X) \n''c'' をこの構造体の変数の名前 (''trainedModel'' など) に置き換えます。 \n \nこのモデルは 2 個の予測子を使用して学習を行ったため、X は厳密に 2 列を含んでいなければなりません。 \nX は、順序と形式が学習データと厳密に同じ予測子列のみを含んでいなければなりません。 \n応答列や、アプリにインポートしなかった列を含めないでください。 \n \n詳細については、<a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a> を参照してください。');

% % 予測子と応答の抽出
% % このコードは、データを処理して、モデルに学習させるのに適した
% % 形状にします。
% % 入力をテーブルへ変換
% inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2'});
% 
% predictorNames = {'column_1', 'column_2'};
% predictors = inputTable(:, predictorNames);
% response = responseData;
% isCategoricalPredictor = [false, false];
% 
% % 交差検証の実行
% KFolds = 5;
% cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% % 予測を適切なサイズに初期化
% validationPredictions = response;
% for fold = 1:KFolds
%     trainingPredictors = predictors(cvp.training(fold), :);
%     trainingResponse = response(cvp.training(fold), :);
%     foldIsCategoricalPredictor = isCategoricalPredictor;
% 
%     % 回帰モデルの学習
%     % このコードは、すべてのモデル オプションを指定してモデルに学習させます。
%     responseScale = iqr(trainingResponse);
%     if ~isfinite(responseScale) || responseScale == 0.0
%         responseScale = 1.0;
%     end
%     boxConstraint = responseScale/1.349;
%     epsilon = responseScale/13.49;
%     regressionSVM = fitrsvm(...
%         trainingPredictors, ...
%         trainingResponse, ...
%         'KernelFunction', 'gaussian', ...
%         'PolynomialOrder', [], ...
%         'KernelScale', 1.4, ...
%         'BoxConstraint', boxConstraint, ...
%         'Epsilon', epsilon, ...
%         'Standardize', true);
% 
%     % 関数 predict で結果の構造体を作成
%     svmPredictFcn = @(x) predict(regressionSVM, x);
%     validationPredictFcn = @(x) svmPredictFcn(x);
% 
%     % 結果の構造体にさらにフィールドを追加
% 
%     % 検証予測の計算
%     validationPredictors = predictors(cvp.test(fold), :);
%     foldPredictions = validationPredictFcn(validationPredictors);
% 
%     % 予測を元の順序で保存
%     validationPredictions(cvp.test(fold), :) = foldPredictions;
% end
% 
% % RMSE 検証の計算
% isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
% validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
