%%% This is a function to build enzyme constrained models and then use
%%% blood medium to properly constrain them. It produces three output files
%%% corresponding to differing levels of nutrient constrain (see paper for
%%% more details)

%Arguments: models -> models to be used 
%           modelids -> model ids
%           output_path -> where all the data should be stored
%Ouputs: ecModels_constrained -> final enzyme constrained models
%        model_1_flux_solution_Prot_restriction_ecModel.txt -> flux
%        solutions for nutrient restricted (a = 1e-5)
%        model_2_flux_solution_Prot_restriction_ecModel.txt -> flux
%        solutions for intermediate nutrient profile (a = 1e-4)
%        model_3_flux_solution_Prot_restriction_ecModel.txt -> flux
%        solutions for nutrient rich (a = 2e-4)
%        model_ecmodel.mat -> saved model
%        model_rxn_list.txt -> List of reaction IDs
function [ecModels_constrained] = build_ecModel_blood_constrained(models, modelids, output_path)
    
    %load the necessary script for GECKO enzyme constraint
    adapterLocation = fullfile(findGECKOroot,'tutorials','light_ecModel','HumanGEMAdapter.m');
    adapter = ModelAdapterManager.setDefault(adapterLocation);
    ecModel_base = loadEcModel; %load base enzyme constrained model which will make building models faster
    ecModels_constrained = {}; %initialise models
    bloodData = prepBloodData(); %prep blood constraints
    cwd = pwd; 
    for k = 1:length(models)
        disp(k)
        current_model = models{k};
        current_model_id = modelids{k};
        model_output_dir = [cwd '/' output_path '/' current_model_id];
        mkdir(model_output_dir)
        %cd(model_output_dir)
        
        %make the enzyme constraint model
        ec_Model = getSubsetEcModel(ecModel_base, current_model);
        ec_Model_irrev = convertToIrreversible(ec_Model);
        write_path = [model_output_dir '/' current_model_id '_ecmodel.mat'];
        writecell(ec_Model_irrev.rxns, [model_output_dir '/' current_model_id '_' 'rxn_list.txt'])
        save(write_path, 'ec_Model_irrev');
        
        %set constants for the radial diffusion model; these have been
        %determined earlier
        a = [1e-5, 1e-4, 2e-4];
        nPoints = length(a);
        params = struct();
        params.relGap = 0.4;
        params.FeasibilityTol = 1e-9;
        params.OptimalityTol = 1e-9;
        for i = 1:nPoints
            ux = zeros(length(bloodData.totDxC),2);
            ux(:,2) = bloodData.totDxC*a(i);

            %insert constraints from blood conc 
            modelGrowth = constrainMedium(ec_Model_irrev, bloodData.totMets, ux, false, true);
            
            %make an initial solution of highest possible biomass
            %accumulation with relaxed protein constraints
            modelGrowth = setParam(modelGrowth, 'ub', 'prot_pool_exchange_r', 20.8356);
            biomass_initial_sol = solveLP(modelGrowth, 1);
            biomass_val = biomass_initial_sol.f;
            
            %now take the previous value and try to minimise the total
            %protein usage to obtain better fluxes
            modelGrowth = setParam(modelGrowth, 'lb', 'MAR13082', 0.99*biomass_val);
            modelGrowth = setParam(modelGrowth, 'obj', 'prot_pool_exchange_r', -1);
            prot_pool_sim = solveLP(modelGrowth, 1);
            %if a viable solution can be found, great. Otherwise, try to
            %relax the constrain a bit further to find an optimal solution
            if prot_pool_sim.stat == 1
                prot_pool_val = -prot_pool_sim.f;
            else
                modelGrowth = setParam(modelGrowth, 'lb', 'MAR13082', 0.98*biomass_val);
                modelGrowth = setParam(modelGrowth, 'obj', 'prot_pool_exchange_r', -1);
                prot_pool_sim = solveLP(modelGrowth, 1);
                %if you still cannot find a solution, relax the protein
                %constrain to the maximum possible value
                if prot_pool_sim.stat == -1
                    prot_pool_val = 20.8356;
                else
                   prot_pool_val = -prot_pool_sim.f;
                end

            end
            %now try to get a viable flux solution while having the
            %previously calculated protein pool. It is slightly relaxed as
            %otherwise mathematically it becomes hard to find a solution
            modelGrowth = setParam(modelGrowth, 'ub', 'prot_pool_exchange_r', 1.03*prot_pool_val);
            modelGrowth = setParam(modelGrowth, 'obj', 'MAR13082', 1);
            biomass_prot_constr = solveLP(modelGrowth, 1, params);
            %similarly relax a bit further if you cannot find a solution
            if biomass_prot_constr.stat == -1
                modelGrowth = setParam(modelGrowth, 'ub', 'prot_pool_exchange_r', 1.1*prot_pool_val);
                modelGrowth = setParam(modelGrowth, 'obj', 'MAR13082', 1);
                biomass_prot_constr = solveLP(modelGrowth, 1, params);
                %if all else fails then just fully relax the protein pool
                %constraint
                if biomass_prot_constr.stat == -1
                    modelGrowth = setParam(modelGrowth, 'ub', 'prot_pool_exchange_r', 20.8356);
                    modelGrowth = setParam(modelGrowth, 'obj', 'MAR13082', 1);
                    biomass_prot_constr = solveLP(modelGrowth, 1, params);
                end
            end
            
            %write the output
            writematrix(biomass_prot_constr.x,[model_output_dir '/' current_model_id '_' num2str(i) '_' 'flux_solution_Prot_restriction_ecModel.txt'])

        end

    ecModels_constrained.(current_model_id)= ec_Model_irrev;
    end



end