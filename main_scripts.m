%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This file is a master script that contains various functions for how to perform subsequent analysis in the paper %%%%%%%%
%%%%%% This script is broken down into various parts described below. Some steps require previous steps/analysis so     %%%%%%%%
%%%%%% please navigate to the appropriate steps to perform the required analysis. Modify variables as appropriate for   %%%%%%%%
%%%%%% your required study/analysis											     %%%%%%%
%%%%%% Step 1: Build base models using ftINIT as part of the RAVEN package						     %%%%%%%
%%%%%% Step 2: Determine what reactions are present in each model as a subset of the original base model		     %%%%%%%
%%%%%% Step 3: Determine what metabolic functions are performed in each model					     %%%%%%%
%%%%%% Step 4: Determine what genes are essential									     %%%%%%%
%%%%%% Step 5: Build enzyme constrained models and simulate them under different conditions			             %%%%%%%
%%%%%% Step 6: Run blood simulation curves for various values							     %%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pre-requisitie packages (please follow instructions on the repositories for installation in MATLAB:
% 1. RAVEN Tooblox (https://github.com/SysBioChalmers/RAVEN)
% 2. GECKO-3 (https://github.com/SysBioChalmers/GECKO)
% 3. Human1 (https://github.com/SysBioChalmers/Human-GEM) 

%%% Step 0: Set the current working directory as you will need to come back to this folder. NB: This must be the folder where all 
%%% the scripts are stored
cwd = pwd;

%%% Step 1: Build base models using ftINIT as part of the RAVEN package. NB: these models are not enzyme constrained and these will need to be done below
modelpath = 'Models/TCGA_Pan'; %Set the model path directory to a folder of interest; models will be saved in this folder and can be loaded again for use later
mkdir(modelpath);
transcriptomic_file = 'TPM_Data/TCGA_Pan/TCGA_Pan_TPM.csv'; %Modify this to the appropriate csv with TPM data

%prepData is required by ftINIT to be able to build models quickly
load('Reference_Models/prepData_new.mat') 
build_models(transcriptomic_file, prepData_new, modelpath)
%%%

%%% Step 2: Determine what reactions are present in each model as a subset of the original base model; this was used for analysis in Figure 2
%First we need to load models from the appropriate directory; You can skip these two steps if you are coming from step 1
load('Reference_Models/prepData_new.mat') 
modelpath = 'Models/TCGA_Pan'; %Models will be retrieved from this directory

%load the models
[all_models, modelids] = loadmodels(modelpath);
all_models = struct2cell(all_models);
cd(cwd);

%Run the command to find what reactions are present and absent in each model
rxn_composition_path = 'Results/rxn_composition/TCGA_Pan'; %output path for this analysis
mkdir(rxn_composition_path)
rxn_composition(prepData_new, all_models, rxn_composition_path)
%%%

%%% Step 3: Determine what metabolic functions are performed in each model; this was used for analysis in Figure 3
%First we need to load models from the appropriate directory; You can skip this step if you have loaded the models
modelpath = 'Models/TCGA_Pan'; %Models will be retrieved from this directory

%load the models
[all_models, modelids] = loadmodels(modelpath);
all_models = struct2cell(all_models);
cd(cwd);

%for metabolic function analysis boundary metabolites are required
for i = 1:length(all_models)
    all_models{i} = addBoundaryMets(all_models{i});
end

%Run the analysis with the function compareMultipleModels; NB this function is part of the RAVEN tooblox that is required to be installed before
model_func_file = ['Results/model_struc_comparison/TCGA_Pan' '/' 'Metabolic_Functions.csv']; %define output file for this
taskFileName = 'Reference_Models/Human1_metabolicTasks/metabolicTasks_Full.txt'; %contains a list of all the metabolic tasks and required outputs
res_func = compareMultipleModels(all_models, false, false, [], true, taskFileName); %run the function to compare models and see what tasks they can perform NB: this can take some time
writematrix(res_func.funcComp.matrix, model_func_file, 'Delimiter', ',') %write a binary mxn matrix in which m corresponds to 257 (the number of metabolic tasks) and n corresponds to the length of the models
%%%

%%% Step 4: Determine what genes are essential; this was used for analysis in Figure 3
%WARNING: This analysis can take a long time. It is recommended to do this on a smaller number of models
%First we need to load models from the appropriate directory; You can skip this step if you have loaded the models
modelpath = 'Models/TCGA_Pan'; %Models will be retrieved from this directory
load('Reference_Models/prepData_new.mat') 

%load the models
[all_models, modelids] = loadmodels(modelpath);
all_models = struct2cell(all_models);
cd(cwd);

%for gene essentiality analysis boundary metabolites are required
for i = 1:length(all_models)
    all_models{i} = addBoundaryMets(all_models{i});
end

%Get essential tasks from file and parse the tasklist before running the analysis
essential_tasks = 'Reference_Models/Human1_metabolicTasks/metabolicTasks_Essential.txt';
essential_tasks_parsed = parseTaskList(essential_tasks);

%Output folder for gene essentiality analysis
essential_gene_path = 'Results/gene_essentiality/TCGA_Pan';
mkdir(essential_gene_path)

%run the function NB: this will take a long time
gene_essentiality(essential_gene_path, essential_tasks_parsed, all_models, modelids, prepData_new)
%%%

%%% Step 5: Build enzyme constrained models and simulate them under different conditions e.g. nutrient rich and nutrient poor; this was used for analysis in Figure 4, 6
%First we need to load models from the appropriate directory; You can skip this step if you have loaded the models
modelpath = 'Models/TCGA_Pan'; %Models will be retrieved from this directory

%load the models
[all_models, modelids] = loadmodels(modelpath);
all_models = struct2cell(all_models);
cd(cwd);

%define the output directory
FBA_ecPath = 'Results/FBA_ecModels/TCGA_Pan';
mkdir(FBA_ecPath)
[all_ecModels] = build_ecModel_blood_constrained(all_models, modelids, FBA_ecPath); %run the function
%%%

%%% Step 6: Run blood simulation curves for various values ; this was used for analysis in Supplementary Figure 5
%First we need to load models from the appropriate directory; You can skip this step if you have loaded the models
modelpath = 'Models/TCGA_Pan'; %Models will be retrieved from this directory

%load the models
[all_models, modelids] = loadmodels(modelpath);
all_models = struct2cell(all_models);
cd(cwd);

%define the output directory
growth_simulation_path = 'Results/growth_simulation/TCGA_Pan';
mkdir(growth_simulation_path)
growth_simulation_blood(all_models, modelids, growth_simulation_path); %run the function
%%%

