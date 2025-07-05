% 2025/5/20 by jiaqiang Liao
clear all,clc

%% dataread
cd('..\2_Interim\')
SNF = imread("SNF_predict_0520.tif");
FNF = imread("FNF_predict_0520.tif");
subplot(1,2,1),imagesc(SNF),title("SNF")
subplot(1,2,2),imagesc(FNF),title("FNF")

% global change factor data from CMIP6
load Delta_globalChange.mat 
CO2 = imresize(delta_CO2,[360 720],'nearest');
Ta = imresize(delta_Ta,[360 720],"nearest");
Pre = imresize(Precent_Pre,[360 720],"nearest");
Ndep = imresize(Precent_Ndep,[360 720],"nearest");

Ta = [Ta(:,361:720),Ta(:,1:360)];
Pre = [Pre(:,361:720),Pre(:,1:360)];
Ndep = [Ndep(:,361:720),Ndep(:,1:360)];

subplot(2,2,1),imagesc(CO2),title('delta CO2'),colorbar
subplot(2,2,2),imagesc(Ta),title('delta Ta'),colorbar
subplot(2,2,3),imagesc(Pre),title('delta Pre'),colorbar
subplot(2,2,4),imagesc(Ndep),title('delta Ndep'),colorbar

%% Calculate future BNF response based on control experiment
% Weighted response coefficient (k) of BNF global change control experiment
% eCO2(FNF,0.18%; SNF,0.01%)
% Ta(FNF,86.20%; SNF,19.93%)
% Pre(FNF,0.97%; SNF,0.93%)
% Ndep(FNF,-2.29%; SNF,-0.38%)

% eCO2
FNF_eCO2 = (FNF.*CO2.*0.18)./100;
SNF_eCO2 = (SNF.*CO2.*0.01)./100;
FNF_eCO2 = flipud(FNF_eCO2);
SNF_eCO2 = flipud(SNF_eCO2);
R = georefcells([-90 90], [-180 180], size(FNF_eCO2));
% geotiffwrite('FNF_eCO2.tif', FNF_eCO2, R);
% geotiffwrite('SNF_eCO2.tif', SNF_eCO2, R);

% Warming
FNF_eT = (FNF.*Ta.*86.20)./100;
SNF_eT = (SNF.*Ta.*19.93)./100;
FNF_eT = flipud(FNF_eT);
SNF_eT = flipud(SNF_eT);
% geotiffwrite('FNF_eT.tif', FNF_eT, R);
% geotiffwrite('SNF_eT.tif', SNF_eT, R);

% Pre
FNF_Pre = (FNF.*Pre.*0.97)./100;
SNF_Pre = (SNF.*Pre.*0.93)./100;
FNF_Pre = flipud(FNF_Pre);
SNF_Pre = flipud(SNF_Pre);
% geotiffwrite('FNF_Pre.tif', FNF_Pre, R);
% geotiffwrite('SNF_Pre.tif', SNF_Pre, R);

% Ndep
FNF_Ndep = (FNF.*Ndep.*(-2.29))./100;
SNF_Ndep = (SNF.*Ndep.*(-0.38))./100;
FNF_Ndep = flipud(FNF_Ndep);
SNF_Ndep = flipud(SNF_Ndep);
% geotiffwrite('FNF_Ndep.tif', FNF_Ndep, R);
% geotiffwrite('SNF_Ndep.tif', SNF_Ndep, R);


%% calculate and plot
clear all, clc
cd('..\2_Interim\')

%% Fig.4B result
load Delta_globalChange.mat 
CO2 = imresize(delta_CO2,[360 720],'nearest');
Ta = imresize(delta_Ta,[360 720],"nearest");
Pre = imresize(Precent_Pre,[360 720],"nearest");
Ndep = imresize(Precent_Ndep,[360 720],"nearest");
Ta = [Ta(:,361:720),Ta(:,1:360)];
Pre = [Pre(:,361:720),Pre(:,1:360)];
Ndep = [Ndep(:,361:720),Ndep(:,1:360)];

Land = imread('Landcover_WGS84.tif');
Land = imresize(Land,[360 720],'nearest');

RR_vector = [0.01, 0.18, 19.93, 86.20, 0.93, 0.97, -0.38, -2.29];
RR_name = {CO2,CO2,Ta,Ta,Pre,Pre,Ndep,Ndep};
Var_name = {'eCOSNF','eCOFNF','eTSNF','eTFNF','preSNF','preFNF','NdepSNF','NdepFNF'};

for i = 1:8
    data = RR_name{i}.*RR_vector(i);
data(Land <1 | Land >= 12) = nan;

dataVector = data(:);
dataVector = dataVector(~isnan(dataVector));
meanValue = mean(dataVector);
SD = std(dataVector);
result(1,i) = meanValue;
result(2,i) = SD;
end
T = table(result(1,:)',result(2,:)','VariableNames',{'Mean (%)','SD'},'RowNames',Var_name);
disp(T);

%% %% Fig.4a result
SNFeCO2 = imread('SNF_eCO2.tif');
FNFeCO2 = imread('FNF_eCO2.tif');
SNFeT = imread('SNF_eT.tif');
FNFeT = imread('FNF_eT.tif');
SNFPre = imread('SNF_Pre.tif');
FNFPre = imread('FNF_Pre.tif');
SNFNdep = imread('SNF_Ndep.tif');
FNFNdep = imread('FNF_Ndep.tif');
load Area_WGS_1984_720_360.mat  %Area,m2
Area = Area_WGS_1984/10000;

Change_BNF = {SNFeCO2,FNFeCO2,SNFeT,FNFeT,SNFPre,FNFPre,SNFNdep,FNFNdep};
for i = 1:8
data = Change_BNF{i};
data(data < -10000) = nan;
data = imresize(data,[360 720],'nearest');
area_BNF = data.*Area;  % kg/yr-1
total_BNF = sum(area_BNF,'all','omitnan'); 
total_BNF = total_BNF*1000*1e-12;% Convert to Tg/yr-1
change_size(i) = total_BNF;
end
T2 = table(change_size(1,:)','VariableNames',{'Size,Tg N'},'RowNames',Var_name);
disp(T2);

%% Calculate SD for future response (Monte Carlo method)
SNF = imread("SNF_predict_0520.tif");
FNF = imread("FNF_predict_0520.tif");
Var_name = {'eCOSNF','eCOFNF','eTSNF','eTFNF','preSNF','preFNF','NdepSNF','NdepFNF'};
RR_vector = [0.01, 0.18, 19.93, 83.56, 0.93, 0.97, -0.38, -2.29];
RR_sd = [0.09, 0.40, 20.94, 101.93, 1.80, 3.70, 1.00, 8.03];
RR_se = [0.01, 0.07, 3.76, 9.34, 0.22, 0.46, 0.15, 0.65];
BNF_map = {SNF,FNF,SNF,FNF,SNF,FNF,SNF,FNF};
GCF_map = {CO2,CO2,Ta,Ta,Pre,Pre,Ndep,Ndep};

for i = 1:8
    BNF_map{i} = imresize(BNF_map{i},[360 720],'nearest');
end
for i = 1:8
    GCF_map{i} = imresize(GCF_map{i},[360 720],'nearest');
end

rng(2024)
num_simulations = 1000;  
change_size_samples = zeros(num_simulations, 8);
for j = 1:num_simulations
    for i = 1:8
        mean_RR = RR_vector(i);
        sd_RR = RR_se(i);
        RR = normrnd(mean_RR,sd_RR);
        sampled_data = (BNF_map{i}.*GCF_map{i}.*RR)./100;
        sampled_data(Land <1 | Land >= 12) = nan;
        sampled_data(sampled_data < -10000) = nan;

        area_BNF = sampled_data .* Area;
        total_BNF = sum(area_BNF, 'all', 'omitnan');
        total_BNF = total_BNF * 1000 * 1e-12;
        change_size_samples(j, i) = total_BNF;
    end
end

TotalchangeSize_samples = sum(change_size_samples, 2); 
mean_total_change = mean(TotalchangeSize_samples);
std_total_change = std(TotalchangeSize_samples);
conf_interval = prctile(TotalchangeSize_samples, [5, 95]);
fprintf('Mean Total Change: %.4f Tg\n', mean_total_change);
fprintf('Standard Deviation: %.4f Tg\n', std_total_change);
fprintf('95%% Confidence Interval: [%.4f, %.4f] Tg\n', conf_interval(1), conf_interval(2));

