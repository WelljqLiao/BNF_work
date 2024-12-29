clc,clear all
%% ------- Data ------- %%
% -*- coding: GBK -*-
% To create the final SNF and FNF maps, we used an ensemble approach,  
% whereby we averaged the global predictions from the 100 best random-forest models based on bootstrap procedure.

datapath = {'1_Data\SNF_extract.csv','1_Data\FNFdata_extract.csv'};
VarName = {'SNF','FNF'};
Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool();
end

% SNF-datapath{1}, FNF-datapath{2}
BNF_ob = readtable(datapath{1});
BNF_ob = BNF_ob(:,1:39);
BNF = table2array(BNF_ob);
BNF = fillmissing(BNF,"movmean",235);
X = BNF(:,5:39);
Y = log(BNF(:,4));

%% Initial model trian
for i = 1:100
    % 80% train data 20% test data
    cv = cvpartition(size(X, 1), 'HoldOut', 0.2);
    idxTrain = training(cv);
    Xtrain = X(idxTrain,:);
    ytrain = Y(idxTrain,:);
    Xtest = X(~idxTrain,:);
    ytest = Y(~idxTrain,:);
    % model train
    numTrees = 50;
    RFLeaf= 5; 
    RFModel = TreeBagger(numTrees, Xtrain, ytrain, 'Method', 'regression', ...
        'OOBPredictorImportance','on',...
        'MinLeafSize',5);

    % variable importance
    importance = RFModel.OOBPermutedPredictorDeltaError;
    [~,idx] = sort(importance,'descend');
    numFeatures = 10;
    selectedFeatures(i,:) = idx(1:numFeatures);
    VarIm(i,:) = RFModel.OOBPermutedPredictorDeltaError;
end
% top 10 drivers
[counts,binEdges] = histcounts(selectedFeatures);
[sortedCounts, idx] = sort(counts, 'descend');
Features = idx(1:10);

%% Further model trian (top 10 drivers)
bestModels = cell(100, 1);
bestR2s = zeros(100, 1);
parfor j = 1:100
    bestR2 = -inf;
    bestModel = []; 
    for i = 1:100
        % 80% train data 20% test data
        cv = cvpartition(size(X, 1), 'HoldOut', 0.2);
        idxTrain = training(cv);
        Xtrain = X(idxTrain,:);
        ytrain = Y(idxTrain,:);
        Xtest = X(~idxTrain,:);
        ytest = Y(~idxTrain,:);
        % model train
        numTrees = 50;
        RFModelSelected = TreeBagger(numTrees, Xtrain(:, Features), ytrain, ...
            'Method', 'regression','OOBPredictorImportance','on',...
        'MinLeafSize',5);
        % model test
        ypred2 = predict(RFModelSelected, Xtest(:, Features));
        mse = mean((ytest - ypred2).^2);
        rsquared = 1 - mse/var(ytest);
        if rsquared > bestR2
            bestR2 = rsquared;
            bestModel = RFModelSelected;
        end
    end
    bestModels{j} = bestModel;
    bestR2s(j) = bestR2;
end
[maxR2, maxIndex] = max(bestR2s);
bestModel = bestModels{maxIndex};
% Model performance
Y_trained = bestModel.Y;
X_trained = bestModel.X;
Y_predicted = predict(bestModel,X_trained);
biasModel = fitlm(Y_predicted, Y_trained);
y_pred = predict(biasModel, Y_predicted);
figure()
scatter(Y_trained,y_pred, 'o', ...
    'MarkerEdgeColor', 'k', 'SizeData', 50,'LineWidth', 1);
axis square
box on;
hold on;
xlim = get(gca, 'XLim');
line([xlim(1), xlim(2)], [xlim(1), xlim(2)], 'Color', 'k', 'LineStyle', ':','LineWidth',1);
hold off;
cor = corr(Y_trained,y_pred);
title_str = [VarName{1},' Model Performance'];
title(title_str,'FontSize',14);
ylabel('Predicted Value (ln, kg ha^-^1 yr^-^1)','FontSize',12,'FontName','Times');
xlabel('Observed Value (ln, kg ha^-^1 yr^-^1)','FontSize',12,'FontName','Times');
legend(num2str(maxR2), '1:1 Line', 'Location', 'northwest');
set(gca,'FontName','Times');
grid on;

%% Model prediction
load '2_ML model'\X_predict(35)_inter.mat
Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);
varNames = Variable_Name(Features);
predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end
Y_predict_all = zeros(size(predict_X, 1), 100);

% Model ensemble result (100 runs)
tic;
parfor j = 1:100
    Y_predict = predict(bestModels{j}, predict_X);
    Y_predict(Landcover_2020 <1 | Landcover_2020 >14) = nan;

    Y_predict_all(:, j) = Y_predict;
end
toc;
% save SNF100model Y_predict_all
% save FNF100model Y_predict_all

Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;

Y_predict_std = std(Y_predict_all, 0, 2);
Y_predict_std = reshape(Y_predict_std,[360,720]);

anss = prctile(Y_predict_avg ,[5,95],'all')
meanY = mean(Y_predict_avg,"all","omitnan");

Y_predict_avg = reshape(Y_predict_avg,[360,720]);
Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
histogram(Y_predict_avg)

Y_predict_avg = exp(Y_predict_avg);
histogram(Y_predict_avg)

BNF_predict = Y_predict_avg;
load '2_ML model'\Area_WGS_1984_720_360.mat  % unit m2
Area = Area_WGS_1984/10000; % unit ha

area_BNF = BNF_predict.*Area;  
total_BNF = sum(area_BNF,'all','omitnan'); 
total_BNF = total_BNF*1000*1e-12; 
disp(['Global amount = ',num2str(total_BNF)]);
