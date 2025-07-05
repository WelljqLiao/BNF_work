% 2025/5/20 by jiaqiang Liao
clc, clear all
% colormap
mycolorpoint=[[253 231 36];...
    [91 200 98];
    [32 143 140];
    [28 82 139];
    [68 1 84]];

mycolorposition=[1 32 64 96 128];
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:128,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:128,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:128,'linear','extrap');
mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
mycolor=round(mycolor*10^4)/10^4;

%% global size and uncertainty
cd('..\2_Interim\')
load Area_WGS_1984_720_360.mat  
Area = Area_WGS_1984/10000; % unit ha
load Landcover_2020.mat

% SNF-globalsize
load SNF_predict0520.mat
area_SNF = BNF_predict.*Area;  
total_BNF = sum(area_SNF,'all','omitnan'); 
total_BNF = total_BNF*1000*1e-12;
disp(['Global SNF amount = ',num2str(total_BNF)]);

SNF_predict = flipud(BNF_predict);
R = georefcells([-90 90], [-180 180], size(SNF_predict));
% geotiffwrite('SNF_predict_0520.tif', SNF_predict, R);

% cv
Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
SNF_CV = Y_predict_std./Y_predict_avg;
SNF_CV = reshape(SNF_CV,[360,720]);
SNF_CV(Landcover_2020 <1 | Landcover_2020 >14) = nan;
mean_cv = mean(SNF_CV,"all","omitnan");

SNF_CV = flipud(SNF_CV);
R = georefcells([-90 90], [-180 180], size(SNF_CV));
% geotiffwrite('SNF_uncer_0520.tif', SNF_CV, R);

% Model ensemble result (100 runs) 
for i = 1:100
    bnf_i = Y_predict_all(:,i);
    Y_predict_avg = reshape(bnf_i,[360,720]);
    Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
    BNF_i = exp(Y_predict_avg)-1;
    mean_bnf(i) = mean(BNF_i,"all","omitnan");
    area_bnf = BNF_i.*Area;
    total_bnf = sum(area_bnf,'all','omitnan');
    resultbnf(i) = total_bnf*1000*1e-12;
end
part_bnf = prctile(mean_bnf,[5,95],'all');
std_rate_bnf = std(mean_bnf);
part_size_bnf = prctile(resultbnf,[5,95],'all');
std_size_bnf = std(resultbnf);
meanbnf = mean(BNF_predict,"all","omitnan");
area_BNF = BNF_predict.*Area;
total_BNF = sum(area_BNF,'all','omitnan');
total_BNF = total_BNF*1000*1e-12;

T = table( ...
    [part_bnf(1); part_size_bnf(1)], ...
    [part_bnf(2); part_size_bnf(2)], ...
    [std_rate_bnf; std_size_bnf], ...
    [meanbnf; total_BNF], ...
    'VariableNames', {'5th_Percentile', '95th_Percentile', 'Standard_Deviation', 'Mean_or_Total'}, ...
    'RowNames', {'Rate', 'Quantity'});
disp('SNF Estimation Results:');
disp(T);

%% FNF
load FNF_predict0520.mat
area_FNF = BNF_predict.*Area;  
total_BNF = sum(area_FNF,'all','omitnan'); 
total_BNF = total_BNF*1000*1e-12; 
disp(['Global FNF amount = ',num2str(total_BNF)]);

BNF_predict = flipud(BNF_predict);
R = georefcells([-90 90], [-180 180], size(BNF_predict));
% geotiffwrite('FNF_predict_0520.tif', BNF_predict, R);

% FNF-CV
Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
FNF_CV = Y_predict_std./Y_predict_avg;
FNF_CV = reshape(FNF_CV,[360,720]);
FNF_CV(Landcover_2020 <1 | Landcover_2020 >14) = nan;
mean_cv = mean(FNF_CV,"all","omitnan");

FNF_CV = flipud(FNF_CV);
R = georefcells([-90 90], [-180 180], size(FNF_CV));
% geotiffwrite('FNF_uncer_0520.tif', FNF_CV, R);

% Model ensemble result (100 runs) 
for i = 1:100
    bnf_i = Y_predict_all(:,i);
    Y_predict_avg = reshape(bnf_i,[360,720]);
    Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
    BNF_i = exp(Y_predict_avg)-1;
    mean_bnf(i) = mean(BNF_i,"all","omitnan");
    area_bnf = BNF_i.*Area;
    total_bnf = sum(area_bnf,'all','omitnan');
    resultbnf(i) = total_bnf*1000*1e-12;
end
part_bnf = prctile(mean_bnf,[5,95],'all');
std_rate_bnf = std(mean_bnf);
part_size_bnf = prctile(resultbnf,[5,95],'all');
std_size_bnf = std(resultbnf);
meanbnf = mean(BNF_predict,"all","omitnan");
area_BNF = BNF_predict.*Area;
total_BNF = sum(area_BNF,'all','omitnan');
total_BNF = total_BNF*1000*1e-12;

T = table( ...
    [part_bnf(1); part_size_bnf(1)], ...
    [part_bnf(2); part_size_bnf(2)], ...
    [std_rate_bnf; std_size_bnf], ...
    [meanbnf; total_BNF], ...
    'VariableNames', {'5th_Percentile', '95th_Percentile', 'Standard_Deviation', 'Mean_or_Total'}, ...
    'RowNames', {'Rate', 'Quantity'});
disp('FNF Estimation Results:');
disp(T);

%% ecosystem statics
[Landcover_2020 R] = readgeoraster('Landcover_WGS84.tif');
Land = Landcover_2020;
Land = imresize(Land,[360,720],'nearest');

fracdata = {'area_SNF','area_FNF'};
for n = 1:17
    cover_idx = find(Land == n);
    for i = 1:2
        Eco_data = eval(fracdata{i});
        cover_map = Eco_data(cover_idx);
        cover_sum = sum(cover_map,"all",'omitnan');
        Ecomean(i,n) = cover_sum*1000*1e-12; 
    end
end
disp('Global sum')
sum(Ecomean,2,'omitnan')

colNames = {'ENF','EBF','DNF','DBF','MF','CS','OS',...
    'WSava','Sava','Grass','Perma','Crop','Urban','Crop&Vet','Snow','Barren','Unclass'};
EcoT = array2table(Ecomean, 'VariableNames', colNames,'RowNames', fracdata);
disp(EcoT)

%% pie plot
% gsum = sum(Ecomean,1,'omitnan'); % bnf
gsum = Ecomean(1,:); % snf
% gsum = Ecomean(2,:); % fnf

gsum(:,6) = gsum(:,6) + gsum(:,7); %shrub
gsum(:,7) = [];
gsum(:,10:16) = []; % non_natural

Label = {'ENF','EBF','DNF','DBF','MIX','SHB','WAS','SAV','GRS'};
p = pie(gsum);
colormap(flipud(mycolor))
hLegend = legend(Label, 'Position', [0.87 0.2 0.1 0.3]);
hLegend.ItemTokenSize = [5 5];
legend('boxoff');
th = findobj(gca, 'Type', 'text');
set(th, 'FontName', 'Times', 'FontSize', 13)
set(hLegend, 'FontName',  'Times', 'FontSize', 11)
set(gcf,'Color',[1 1 1])