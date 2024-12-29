clc,clear all
% ------- 1 load environmental variables  ------- %%
% Create a geographic cells reference object
G = georefcells([-90 90],[-180 180],[360 720],'ColumnsStartFrom','north');

cd("2_ML model\")
load Climate_data.mat
load Soil_data.mat
load Vegetation_data.mat
load Ecosystem_data.mat
Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',... 
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);

%% ------- Biological N fixation data  ------- %%
filename = '\1_Data\SNF.csv';   % SNF or FNF
BNF_siteTable = readtable(filename);
BNF_site = table2array(BNF_siteTable);
BNF_RowsNum = height(BNF_site); 

for i = 1:BNF_RowsNum     
latitude = BNF_site(i,2); 
longitude = BNF_site(i,3); 
[row, col] = geographicToDiscrete(G, latitude, longitude);
    for j = 1:variablesNum 
        eval(['value = ',char(Variable_Name(j)),';']);
        if row >= 1 && row <= size(value, 1) && col >= 1 && col <= size(value, 2)
            BNF_grid(i,j) = value(row, col);
        else
            BNF_grid(i,j) = nan;
        end
    end
end

BNF_Variable = array2table(BNF_grid, 'VariableNames', Variable_Name);
BNFdataTable = [BNF_siteTable,BNF_Variable];
savePath = '\1_Data\'; 
filename = 'SNF_extract.csv';
fullPath = fullfile(savePath, filename);
writetable(BNFdataTable, fullPath);
