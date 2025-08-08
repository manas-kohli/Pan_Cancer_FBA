%%% This is a function in order to be able to build models using the ftINIT
%%% algorithm. 

%Arguments: transcriptomic_file -> the file which has TPM data for each
%                                  sample
%           prepData -> input object required to be prepared for ftINIT
%           output_path -> directory for where models need to be saved
%           threshold_vector -> optional argument that is by default set to
%                               1. It can either be set to a different
%                               value or set as a vector equal to the
%                               number of rows in the transcriptomic_file

function build_models(transcriptomic_file, prepData, output_path, threshold_vector)
    % set threshold_vector to a default value of 1 if it's not an argument
    if nargin < 4
        threshold_vector = 1.0;
    end
    transcriptomic_data = readtable(transcriptomic_file, 'ReadVariableNames',true);
    % extract the tissue and gene names
    data_struct.tissues = transcriptomic_data.Properties.VariableNames(2:end)';  % sample (tissue) names
    data_struct.genes = transcriptomic_data.Gene_ID;  % gene names
    data_struct.levels = table2array(transcriptomic_data(:, 2:end));  % gene TPM values
    data_struct.threshold = threshold_vector;

    for k = 1:length(data_struct.tissues)
        tissue = data_struct.tissues{k};  % must match the tissue name in data_struct.tissues
        disp(['Model: ' num2str(k) ' of ' num2str(length(data_struct.tissues))])
       
        built_model = ftINIT(prepData, tissue, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
        built_model.id = tissue;
        
        write_path = [pwd '/' output_path '/' tissue '_model.mat'];
        save(write_path, 'built_model');
    end

end