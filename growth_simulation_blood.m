%%% This is a function to determine growth for various values of a to
%%% determine the optimal constants to use for later

%Arguments: models -> mx1 cell array of models
%           modelsids -> mx1 cell array of modelids
%           output_path -> where data is saved
%Ouputs: model_biomass_vals.csv -> list of biomass values

function growth_simulation_blood(models, modelids, output_path)
    
    %load the necessary script for GECKO enzyme constraint
    adapterLocation = fullfile(findGECKOroot,'tutorials','light_ecModel','HumanGEMAdapter.m');
    adapter = ModelAdapterManager.setDefault(adapterLocation);
    ecModel_base = loadEcModel; %load base enzyme constrained model
    bloodData = prepBloodData(); %prep blood constraints
    cwd = pwd;

    for k = 1:length(models)
        disp(k)
        %load the model
        current_model = models{k};
        current_model_id = modelids{k};
        model_output_dir = [cwd '/' output_path '/' current_model_id];
        mkdir(model_output_dir)
        
        %generate the enzyme constrained model
        ec_Model = getSubsetEcModel(ecModel_base, current_model);
        ec_Model_irrev = convertToIrreversible(ec_Model);
        
        %define the range for a as well as other parameters for simulation
        a = (0:2e-6:2e-4);
        nPoints = length(a);
        params = struct();
        params.relGap = 0.4;
        params.FeasibilityTol = 1e-9;
        params.OptimalityTol = 1e-9;
        modelgrowth = {};
        %loop over various values of a
        for i = 1:nPoints
            %for each value of a set the appropriate constraints for each
            %nutrient
            ux = zeros(length(bloodData.totDxC),2);
            ux(:,2) = bloodData.totDxC*a(i);
            %constraint the medium appropriately
            modelGrowth = constrainMedium(ec_Model_irrev, bloodData.totMets, ux, false, true);
            modelGrowth = setParam(modelGrowth, 'ub', 'prot_pool_exchange_r', 20.8356);
            %generate a biomass value for growth given the constraints
            biomass_initial_sol = solveLP(modelGrowth, 0, params);
            if biomass_initial_sol.stat == 1
                modelgrowth{i} = biomass_initial_sol.f;
            else
                modelgrowth{i} = NaN;
            end
        end
        %write data to an output file
        fileID = fopen([model_output_dir '/' current_model_id '_' 'biomass_vals.txt'], 'w');
        for i = 1:length(modelgrowth)
            fprintf(fileID, '%s\n', mat2str(modelgrowth{i})); % Convert value to string and write
        end
        fclose(fileID);
    end
end