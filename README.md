# 🌱 Increasing Biological Nitrogen Fixation Under Climate Change

**🧑‍🔬 Author**: Jiaqiang Liao  
**📅 Date**: July 4, 2025  

---

## 🧭 Overview

This repository supports the study **“Increasing Biological Nitrogen Fixation Under Climate Change.”**  
It explores:

- Biological nitrogen fixation (BNF) rates
- Responses to global change factors
- Future predictions under climate change scenarios

📁 Includes: datasets, analysis scripts, model outputs, and final results.

---

## 📁 Repository Structure

### `0_Data/`
- `BNFdata_use.xlsx`: Observed nitrogen fixation rates for various nitrogen-fixer types, with references.
- `BNFGlobalchange_experiment.xlsx`: BNF responses (weight k) to global change factors from field control experiment.

### `1_Code/`
- `BNF_Variable_importance.R`: VIF selection and predictor importance (R)
- `BNF_Kfold_model.m`: Train/predict BNF using a random forest model (MATLAB)
- `BNF_statics.m`: Calculate global BNF size and ecological metrics (MATLAB)
- `BNF_futureRes.m`: Future BNF responses under climate scenarios (MATLAB)
- `BNF_Currentmap.R`: Generate current BNF prediction maps (R)
- `BNF_Futuremap.R`: Generate future BNF prediction maps (R)

### `2_Interim/`
- Intermediate variables and model outputs
- Please download these five large files output from [https://drive.google.com/drive/folders/19cKfZqTG4hHND93NPL4ASK3JP3EDmgYL?usp=drive_link] and then put them back into ./2-interim.
`FNF_predict0520.mat`,`FNFallmodel_0520.mat`,`SNF_predict0520.mat`,`SNFallmodel_0520.mat`,`X_predict(36)_inter.mat`

### `3_Results/`
- Final figures (Figures 1 to 4)

---

## 🛠️ Software Requirements

| Software | Version  |
|----------|----------|
| R        | 4.3.3    |
| MATLAB   | R2021b   |
