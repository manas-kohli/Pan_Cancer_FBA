function exchModel = constrainMedium(model, metaboliteNames, fluxLimits, checkFeas, includesComp)
% constrainMedium
%
% Constrain ec model uptake. Run setGrowthMedium first on the model.
%
% Input:
%
%   model          A model derived from Human-GEM.
%
%   metaboliteNames names of the metabolites to constrain           
%
%   fluxLimits      matrix with two columns, min and max; one row per
%                   metabolite
%
% Output:
%
%   exchModel      Model with updated growth medium constraints.
%

if nargin < 5 || isempty(includesComp)
    includesComp = false; %if [e] should be added to the metabolites or not
end

if nargin < 4 || isempty(checkFeas)
    checkFeas = true;
end

[exchRxns, exchIndxs] = getExchangeRxns(model);
% exclude protein pool exchange
[~,prot_exch_indx] = ismember({'prot_pool_exchange', 'f_prot_pool_exchange', 'o_prot_pool_exchange'}, model.rxns);
if (length(prot_exch_indx) == 0)
    error('Expected protein pool exchange reaction named "prot_pool_exchange" was not found.');
else
    exchIndxs = setdiff(exchIndxs, prot_exch_indx, 'stable');
    exchRxns  = setdiff(exchRxns, {'prot_pool_exchange', 'f_prot_pool_exchange', 'o_prot_pool_exchange'}, 'stable');
end

% differentiate between uptake and production reactions
uptkIndxs = exchIndxs(contains(exchRxns, '_b'));
prodIndxs = exchIndxs(find(~contains(exchRxns,'_b')));
exchModel = setParam(model,'ub',prodIndxs,1000);
%close all uptakes
exchModel = setParam(exchModel,'ub',uptkIndxs,0);

% open uptake of media components one by one
unusedMets = [];
for i = 1:length(metaboliteNames)

    % get metabolite index
    if (includesComp)
        metIndx = getIndexes(model, metaboliteNames{i}, 'metcomps');
    else
        metIndx = getIndexes(model, strcat(metaboliteNames{i},'[e]'), 'metcomps');
    end

    % get rxns for metabolite
    metRxns = find(model.S(metIndx,:));

    % get the uptake reaction for the metabolite
    metUptakeRxn = intersect(metRxns,uptkIndxs);
    if isempty(metUptakeRxn)
        unusedMets = [unusedMets; metaboliteNames(i)];
    else
        exchModel.lb(metUptakeRxn) = fluxLimits(i,1);
        exchModel.ub(metUptakeRxn) = fluxLimits(i,2);
    end
end


% report unused metabolites
if ~isempty(unusedMets)
    fprintf('WARNING: The following metabolites are either not in the model or do not have exchange reactions:\n');
    fprintf('\t%s\n',unusedMets{:});
end

%Check if model is feasible
if ( checkFeas )
    sol = solveLP(exchModel);
    if ~isempty(sol.x)
        disp(['Constrained ec model "' exchModel.id '" is feasible'])
    else
        disp(['*** Constrained ec model "' exchModel.id '" is INFEASIBLE ***'])
        exchModel = [];
    end
end