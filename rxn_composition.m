%%% This is a function in order to be able to find reactions present in
%%% each model. 

%Arguments: models -> cell array with all the models as a (nx1) shaped
%                     array
%           prepData -> input object required that contains the base model
%           output_path -> directory for where models need to be saved
%Outputs: compMat -> mxn binary matrix for m reactions and n models saved
%                    to a file in output_path folder
%         baseModel.rxns -> names of reactions for rows in the matrix saved
%                           to a file in output_path folder

function rxn_composition(prepData, models, output_path)

    baseModel = prepData.refModel; %retrieve the base model from prepData
    compMat = false(length(baseModel.rxns), length(models)); %initialise empty matrix

    for i = 1:size(compMat, 2)
        compMat(:,i) = ismember(baseModel.rxns, models{i}.rxns); %populate matrix
    end
    
    %write outputs
    writematrix(compMat,[pwd '/' output_path '/' '_rxn_comp_current.txt'])
    writecell(baseModel.rxns, [pwd '/' output_path '/' 'rxn_list_basemodel.txt'])
end