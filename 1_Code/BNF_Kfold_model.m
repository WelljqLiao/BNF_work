% 2025/5/20 by jiaqiang Liao
clc,clear
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool();
end

% select variable VIF<10
varNames = {'MAP_season','pH','SWC','TP','NPP','BNPP','Ndep','LNC', ...
    'MAT','MAP','AET','Nfixer','TN','SOC'};

%% ------------------- 1.Data Input ------------------- %%
cd('..\0_Data\')
VarName = {'Symbiotic','Free-living'};
ID = 1; % choose SNF{1}, FNF{2}
BNF_ob = readtable('BNFdata_use.xlsx','sheet',VarName{ID});
X = table2array(BNF_ob(:, varNames));
Y0 = table2array(BNF_ob(:, 'BNF'));
Y = log(Y0+1);

latitudes = table2array(BNF_ob(:, 'lat'));
longitudes = table2array(BNF_ob(:, 'lon'));

num_vars = size(X, 2);
VIF = zeros(num_vars, 1);
for i = 1:num_vars
    y = X(:, i); 
    X_other = X(:, [1:i-1, i+1:end]); 
    mdl = fitlm(X_other, y);
    Rsq = mdl.Rsquared.Ordinary;
    VIF(i) = 1 / (1 - Rsq);
end
threshold = 10;
collinear_vars = find(VIF < threshold);
disp(' VIF < 10：');
disp(varNames(collinear_vars));
Features = collinear_vars;

%% ------------------- 2.Model Train ------------------- %%
% ------------------- 2.1 Train based on selected predictors ------------------- %%
K = 5; 
bestR2s = zeros(100, 1);
bestModels = cell(100, 1);
parfor j = 1:100
    cv = cvpartition(size(X, 1), 'KFold', K);
    r2_scores = zeros(K, 1);
    bestR2 = -inf;
    bestModel = [];
    for k = 1:K
        idxTrain = training(cv, k);
        idxTest = test(cv, k);
        Xtrain = X(idxTrain, :);
        ytrain = Y(idxTrain, :);
        Xtest = X(idxTest, :);
        ytest = Y(idxTest, :);
        RFModel = TreeBagger(50, Xtrain(:, Features), ytrain, ...
            'Method', 'regression', 'OOBPredictorImportance', 'on', 'MinLeafSize', 5);
        ypred = predict(RFModel, Xtest(:, Features));
        mse = mean((ytest - ypred).^2);
        r2_scores(k) = 1 - mse / var(ytest);
        if r2_scores(k) > bestR2
            bestR2 = r2_scores(k);
            bestModel = RFModel;
        end
    end
    bestModels{j} = bestModel;
    bestR2s(j) = mean(r2_scores); 
end
[maxR2, maxIndex] = max(bestR2s);
bestModel = bestModels{maxIndex};

cd('..\2_Interim\')
% save SNFallmodel_0520 bestR2s bestModel bestModels maxR2
% save FNFallmodel_0520 bestR2s bestModel bestModels maxR2

%% ------------------- 2.3 Model performance ------------------- %%
load SNFallmodel_0520.mat
% load FNFallmodel_0520.mat
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
axis([-0.9 6.5 -0.9 6.5]) % SNF
% axis([-0.2 4.5 -0.2 4.5]) % FNF
xlim = get(gca, 'XLim');
line([xlim(1), xlim(2)], [xlim(1), xlim(2)], 'Color', 'k', 'LineStyle', ':','LineWidth',1);
hold off;
cor = corr(Y_trained,y_pred);
title_str = [VarName{ID},' Model Performance'];
title(title_str,'FontSize',14);
ylabel('Predicted Value (ln, kg ha^-^1 yr^-^1)','FontSize',12,'FontName','Times');
xlabel('Observed Value (ln, kg ha^-^1 yr^-^1)','FontSize',12,'FontName','Times');
legend(num2str(maxR2), '1:1 Line', 'Location', 'northwest');
set(gca,'FontName','Times');
grid on;

%% ------------------- 3 Model Predict ------------------- %%
load X_predict(36)_inter.mat
load Landcover_2020.mat
load mycolor.mat

varNames = {'MAP_season','pH','SWC','TP','NPP','BNPP','Ndep','LNC_MODIS', ...
    'MAT','MAP','AET','Nfixer','TN','SOC'};
predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end

Y_predict_all = zeros(size(predict_X, 1), 10);
parfor j = 1:100
    Y_predict = predict(bestModels{j}, predict_X);
    Y_predict_all(:, j) = Y_predict;
end
Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
Y_predict_std = reshape(Y_predict_std,[360,720]);
anss = prctile(Y_predict_avg ,[5,95],'all')
meanY = mean(Y_predict_avg,"all","omitnan");
Y_predict_avg = reshape(Y_predict_avg,[360,720]);
Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
BNF_predict = exp(Y_predict_avg)-1;

% save SNF_predict0520 BNF_predict Y_predict_all Y_predict_std
% save FNF_predict0520 BNF_predict Y_predict_all Y_predict_std

%  Global Size
load Area_WGS_1984_720_360.mat  % unit m2
Area = Area_WGS_1984/10000; % unit ha
area_BNF = BNF_predict.*Area;
total_BNF = sum(area_BNF,'all','omitnan');
total_BNF = total_BNF*1000*1e-12;
disp(['Global amount = ',num2str(total_BNF)]);

%% ------------------- 4 AOA spatial extrapolation degree test ------------------- %%
% Calculate the maximum and minimum values of each variable
obs_var = X(:, Features);
max_values = zeros(1, length(Features));
min_values = zeros(1, length(Features));

for i = 1:length(Features)
    max_values(i) = max(obs_var(:, i));
    min_values(i) = min(obs_var(:, i));
end

predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end
predict_X(Landcover_2020 <1 | Landcover_2020 >14,:) = nan;

for i = 1:length(Features)
    results(:, i) = predict_X(:, i) >= min_values(i) & predict_X(:, i) <= max_values(i);
end

% Convert logical values to 0 and 1
results(results == false) = 0;
results(results == true) = 1;

% Average the variable result and calculate the proportion
interpolation = sum(results,2,'omitnan')./height(Features);
interpolation(Landcover_2020 <1 | Landcover_2020 >14,:) = nan;
inter_map = reshape(interpolation,[360,720]);

% figure plot
mapdata = flipud(inter_map);
figure()
lat = [-89.5:0.5:90];
lon = [-179.5:0.5:180];
[lon, lat] = meshgrid(lon,lat);
gca = axesm('MapProjection','robinson','MapLatLimit',[-90 90],'Frame','on','Grid','on', ...
    'FontName','Times','FontSize',24,'FEdgeColor','white');
setm(gca, 'FLineWidth', 0.5);
setm(gca, 'GLineWidth', 0.2);
setm(gca,'GLineStyle','-.');
setm(gca, 'GColor', '#545454');
i1 = surfm(lat, lon, mapdata);
load coastlines
plotm(coastlat,coastlon,'Color','k')
h1 = colorbar('FontName', 'Times', 'FontSize', 12,'Location','southoutside',...
    'Position',[0.25 0.065 0.53 0.05]);
colormap("winter")
axis off;
tightmap;

nanID = ~isnan(interpolation);
inter_box = interpolation(nanID);
figure(),histogram(inter_box,'Normalization','probability')
