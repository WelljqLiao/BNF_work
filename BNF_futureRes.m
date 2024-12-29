% 2024/11/22 by jiaqiang Liao
clear all,clc
cd '2_ML model'\
SNF = imread("SNF_predict.tif");
FNF = imread("FNF_predict.tif");

load Delta_globalChange.mat 
CO2 = imresize(delta_CO2,[1800 3600],'nearest');
Ta = imresize(delta_Ta,[1800 3600],"nearest");
Pre = imresize(Precent_Pre,[1800 3600],"nearest");
Ndep = imresize(Precent_Ndep,[1800 3600],"nearest");

Ta = [Ta(:,1801:3600),Ta(:,1:1800)];
Pre = [Pre(:,1801:3600),Pre(:,1:1800)];
Ndep = [Ndep(:,1801:3600),Ndep(:,1:1800)];

subplot(2,4,1),imagesc(CO2),title('delta CO2'),colorbar
subplot(2,4,2),imagesc(Ta),title('delta Ta'),colorbar
subplot(2,4,3),imagesc(Pre),title('delta Pre'),colorbar
subplot(2,4,4),imagesc(Ndep),title('delta Ndep'),colorbar
subplot(2,4,5),imagesc(SNF),title("SNF"),colorbar
subplot(2,4,6),imagesc(FNF),title("FNF"),colorbar

%% Fig4.a-h
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
%
clear all,clc

%% Fig.4i
load Delta_globalChange.mat 
CO2 = imresize(delta_CO2,[1800 3600],'nearest');
Ta = imresize(delta_Ta,[1800 3600],"nearest");
Pre = imresize(Precent_Pre,[1800 3600],"nearest");
Ndep = imresize(Precent_Ndep,[1800 3600],"nearest");
Ta = [Ta(:,1801:3600),Ta(:,1:1800)];
Pre = [Pre(:,1801:3600),Pre(:,1:1800)];
Ndep = [Ndep(:,1801:3600),Ndep(:,1:1800)];
% eCO2(FNF,0.18%; SNF,0.01%)
% Ta(FNF,86.20%; SNF,19.93%)
% Pre(FNF,0.97%; SNF,0.93%)
% Ndep(FNF,-2.29%; SNF,-0.38%)

Land = imread('Landcover_WGS84.tif');
Land = imresize(Land,[1800 3600],'nearest');

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

%% %% Fig.4j
SNFeCO2 = imread('SNFeCO2-0.1.tif');
FNFeCO2 = imread('FNFeCO2-0.1.tif');
SNFeT = imread('SNFeT-0.1.tif');
FNFeT = imread('FNFeT-0.1.tif');
SNFPre = imread('SNFPre-0.1.tif');
FNFPre = imread('FNFPre-0.1.tif');
SNFNdep = imread('SNFNdep-0.1.tif');
FNFNdep = imread('FNFNdep-0.1.tif');

load Area_WGS_1984_720_360.mat  
Area = Area_WGS_1984/10000; 

Change_BNF = {SNFeCO2,FNFeCO2,SNFeT,FNFeT,SNFPre,FNFPre,SNFNdep,FNFNdep};
for i = 1:8
data = Change_BNF{i};
data(data < -10000) = nan;
data = imresize(data,[360 720],'nearest');
area_BNF = data.*Area;  
total_BNF = sum(area_BNF,'all','omitnan'); 
total_BNF = total_BNF*1000*1e-12; 
change_size(i) = total_BNF;
end
T2 = table(change_size(1,:)','VariableNames',{'Size,Tg N'},'RowNames',Var_name);
disp(T2);
