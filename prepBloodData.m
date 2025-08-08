function bloodData = prepBloodData(showPlot, exportData, useBloodFlowModel)
% Prepares the blood data by loading it from disk, fitting a linear model to 
% the diffusion coefficients, etc.
%
% Make sure the path is in the project base before continuing
%

if (nargin < 1)
    showPlot = false;
end

if (nargin < 2)
    exportData = false;
end

if (nargin < 3)
    useBloodFlowModel = false;
end

%read metabolites and constraints from the excel file with metabolite
%concentrations in blood
bloodConcTable = readtable('Reference_Models/bloodData.csv');

metabolites = bloodConcTable{:,1};
bloodConc = bloodConcTable{:,2};
diffCoeff = bloodConcTable{:,3};
mw = bloodConcTable{:,4};
%We do not have all metabolites. Estimate those from the others using a
%linear model on MW
filt = ~isnan(diffCoeff);
%don't use diffusion coefficients for oxygen and albumin in the model, they
%have very different MW - the linearity only holds in small intervals
filt2 = filt & ~ismember(metabolites, {'O2','NEFA blood pool in','cholesterol','albumin'});
x = mw(filt2);
y = diffCoeff(filt2);
fit = fitlm(x,y);
fit %show R^2 for the fit, 0.665
%save to file for plot in R

if showPlot
    figure
    plot(fit)
    xlabel('Molecular weight')
    ylabel('Diffusion constant')
end
if exportData
    %export the x and y to file so we can plot it in R
    plotExp = struct;
    plotExp.x = x;
    plotExp.y = y;
    plotExp.mets = metabolites(filt2);
    save('data/diffCoeff.mat', 'plotExp')
end


diffCoeff(~filt) = predict(fit, mw(~filt));
if ~useBloodFlowModel
    CxD = bloodConc.*diffCoeff;
else
    %Here, we don't multiply with the diffusion coefficient.
    %We also increase the oxygen concentration
    CxD = bloodConc;
    CxD(strcmp(metabolites,'O2')) = 9200; %9200 (9.2 mmolar)comes from "The oxygen status of the arterial blood revised: relevant oxygen parameters for monitoring the arterial oxygen availability"
end




% hamsMediaMets ={'glucose'
%             'arginine'
%             'histidine'
%             'lysine'
%             'methionine'
%             'phenylalanine'
%             'tryptophan'
%             'tyrosine'
%             'alanine'
%             'glycine'
%             'serine'
%             'threonine'
%             'aspartate'
%             'glutamate'
%             'asparagine'
%             'glutamine'
%             'isoleucine'
%             'leucine'
%             'proline'
%             'valine'
%             'cysteine'
%             'thiamin'
%             'hypoxanthine'
%             'folate'
%             'biotin'
%             'pantothenate'
%             'choline'
%             'inositol'
%             'nicotinamide'
%             'pyridoxine'
%             'riboflavin'
%             'thymidine'
%             'aquacob(III)alamin'
%             'lipoic acid'
%             'sulfate'
%             'linoleate'
%             'linolenate'
%             'O2'
%             'H2O'
%             'retinoate'
%             'Fe2+'
%             'Pi'
%             'alpha-tocopherol'
%             'gamma-tocopherol'};

%check which metabolites in ham's media that are not present in the blood
%data
%uniqueInHam = setdiff(hamsMediaMets, metabolites);
%uniqueInHam 
%There are a lot of vitamins etc in the Ham's media that causes problems if left
%unconstrained. We use the strategy to set these to a small value, such as
%the flux of the metabolite with the smallest detected concentration in blood 
%times 10 times the proportionality constant.

constrainedHamMets={
    'alpha-tocopherol'
    'aquacob(III)alamin'
    'biotin'
    'folate'
    'gamma-tocopherol'
    'inositol'
    'linoleate'
    'linolenate'
    'lipoic acid'
    'nicotinamide'
    'pantothenate'
    'pyridoxine'
    'retinoate'
    'riboflavin'
    'thiamin'
    'thymidine'
    'hypoxanthine'};

ions = {'Na+'
        'K+'
        'Ca2+'
        'NH4+'
        'chloride'
        'HCO3-'};

bloodData.totDxC = [CxD;repmat(min(CxD)*10, length(constrainedHamMets),1);repmat(1000, length(ions),1)];
bloodData.totMets = [metabolites; constrainedHamMets; ions];
bloodData.totMets = strcat(bloodData.totMets, '[e]');
%table(bloodData.totMets, bloodData.totDxC)
