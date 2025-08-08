%%% This is a function to determine essential genes for essential metabolic
%%% functions

%Arguments: output_path -> directory for where results are saved NB it
%                          will contain sub_folders corresponding to each model
%           taskStruc -> tasks parsed that will be used to evaluate gene
%           essentiality
%           models -> mx1 cell array of models
%           modelsids -> mx1 cell array of modelids
%           prepData -> data structure with reference model
%Ouputs: x_genes_all.csv -> genes evaluated during gene essentiality
%        gene_essentiality_matrix -> binary matrix of n rows of genes anx 
%                                    57 columns for essential gene tasks

function gene_essentiality(output_path, taskStruc, models, modelids, prepData)
    % add boundary metabolites for reference model
    refModel = addBoundaryMets(prepData.refModel);
    % loop through models
    for i =1:length(models)
        current_model = addBoundaryMets(models{i}); %add boundary metabolites (if they have been added this will safely pass over anyway)
       
        regen_model = regen_tINIT_model(current_model, refModel, modelids{i}); %regenerate the original tINIT instead of the ftINIT model; this is necessary to avoid bugs

        [~, essentialGenes] = checkTasksGenes(regen_model, [], true, false, true, taskStruc); %check essential genes
        
        %display model id for output progress
        disp(modelids{i})
        
        %write output for model
        mkdir([output_path '/' modelids{i}])
        write_path_genes = [output_path '/' modelids{i} '/' 'genes_all.csv'];
        writecell(regen_model.genes, write_path_genes)
        write_path_essential_genes = [output_path '/' modelids{i} '/' 'gene_essentiality_matrix.csv'];
        writematrix(essentialGenes, write_path_essential_genes, 'Delimiter', ',')
    end

end