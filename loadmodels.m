%%% This is a function to load models built from a directory

%Arguments: directory -> directory for where models need to be loaded
%Ouputs: all_models -> struct object that will have all the models
%        modelids -> cell array thas has all the model id names in order of
%                    how the models are loaded

function [all_models, modelids] = loadmodels(directory)
        
    cd(directory) %change to the relevant directory
    files = dir('*'); %retrieve file names
    all_models = struct(); %initialise structure that will have all the models
    modelids = {}; %initialise modelid names for later
    for k = 3:length(files) %index starts from 3 due to how matlab retrieves file names
        current_model = getfield(files, {k}, 'name'); 
        model_name = erase(current_model, '_model.mat');
        T = load(current_model);
        T = T.built_model;
        all_models.(model_name) = T; %save model to struct
        modelids = [modelids; T.id]; %save model id to cell array
    end
end