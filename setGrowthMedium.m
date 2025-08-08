function exchModel = setGrowthMedium(model, isEC, medium, measuredMets, fluxes, fluxTol, maxBound, checkFeas)
% setGrowthMedium
%
% Set the growth medium for a model derived from Human-GEM. This function
% works for either standard or enzyme constrained-GEMs. The code was supplied from 
% Jonathan Robinson and modified for this project.
%
% Input:
%
%   model          A model derived from Human-GEM.
%
%   isEC           TRUE if the model is an enzyme-constrained GEM (ecGEM)
%                  generated using the GECKO framework.
%                  (Default = FALSE)
%
%   medium         Specifies the growth conditions.
%                   'Hams'  (Default) Only allow uptake of metabolites
%                           present in Ham's media.
%                   'All'   Allow uptake of all extracellular metabolites.
%
%   measuredMets   (Optional) A list of metabolite names (corresponding
%                  to those in model.metNames) for which exchange fluxes
%                  are provided (in the FLUXES input).
%
%   fluxes         (Optional) Measured metabolite exchange fluxes [mmol/gDw h].
%                  If provided, the bounds of the corresponding exchange
%                  reactions will be fixed to the provided flux value.
%                   NOTE: NEGATIVE flux = CONSUME
%                         POSITIVE flux = PRODUCE
%
%   fluxTol        (Optional) Fraction by which the metabolite exchange
%                  flux bounds are allowed to differ from the measured
%                  value. If set to zero, the bounds of the corresponding
%                  exchange reactions will be set to EXACTLY the measured
%                  flux value (not recommended, as this can lead to solver
%                  errors).
%                  (Default = 0.001; i.e., allows a 0.1% error)
%
%   maxBound       (Optional) Default absolute maximum bound value to set
%                  for open/unbounded exchange reactions.
%                  (Default = model.annotation.defaultUB, if it exists,
%                             otherwise 1000)
%                  NOTE! The default max bound of other reactions will NOT
%                  be changed.
%
%   checkFeas      (Optional) TRUE if a test should be run to verify if the
%                  resulting model is feasible.
%                  (Default = TRUE)
%
% Output:
%
%   exchModel      Model with updated growth medium constraints.
%


if nargin < 2 || isempty(isEC)
    isEC = false;
end
if nargin < 3 || isempty(medium)
    medium = 'Hams';
end
if nargin < 4 || isempty(measuredMets)
    measuredMets = {};
end
if nargin < 5
    fluxes = [];
end
if nargin < 6 || isempty(fluxTol)
    fluxTol = 0.001;
end
if nargin < 7 || isempty(maxBound)
    if isfield(model, 'annotation') && isfield(model.annotation, 'defaultUB')
        maxBound = model.annotation.defaultUB;
    else
        maxBound = 1000;
    end
end
if nargin < 8 || isempty(checkFeas)
    checkFeas = true;
end

% check if boundary metabolites are present, if so then remove them
boundaryIndx = find(strcmpi(model.compNames, 'Boundary'));
if ~isempty(boundaryIndx)
    boundary     = find(model.metComps==boundaryIndx);
    if ~isempty(boundary)
        model = removeMets(model, boundary, false, false, false, true);
    end
end
% remove unconstrained field, if present
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
end

if startsWith(lower(medium), 'ham')
    mediaMets ={'glucose'
                'arginine'
                'histidine'
                'lysine'
                'methionine'
                'phenylalanine'
                'tryptophan'
                'tyrosine'
                'alanine'
                'glycine'
                'serine'
                'threonine'
                'aspartate'
                'glutamate'
                'asparagine'
                'glutamine'
                'isoleucine'
                'leucine'
                'proline'
                'valine'
                'cysteine'
                'thiamin'
                'hypoxanthine'
                'folate'
                'biotin'
                'pantothenate'
                'choline'
                'inositol'
                'nicotinamide'
                'pyridoxine'
                'riboflavin'
                'thymidine'
                'aquacob(III)alamin'
                'lipoic acid'
                'sulfate'
                'linoleate'
                'linolenate'
                'O2'
                'H2O'
                'retinoate'
                'Fe2+'
                'Pi'
                'alpha-tocopherol'
                'gamma-tocopherol'};
elseif strcmpi(medium, 'all')
    % open exchange of all metabolites
    [~, comp_ind] = ismember('e', model.comps);
    mediaMets = unique(model.metNames(model.metComps == comp_ind));
else
    error('MEDIUM option not recognized. Valid options are "Hams" or "All".');
end

% set default flux bounds (LB, UB)
fluxBounds = [-ones(length(mediaMets),1), ones(length(mediaMets),1)] * maxBound;

% check if provided mets are part of media's formulation
if ~isempty(measuredMets)
    % modify fluxBounds with the provided flux measurements
    [iA,iB] = ismember(measuredMets, mediaMets);
    % force flux to be equal to measured value (within tolerance fluxTol)
    fluxBounds(iB(iA),:) = [fluxes(iA) - abs(fluxTol*fluxes(iA)), fluxes(iA) + abs(fluxTol*fluxes(iA))];
    if any(~iA)
        %If measured mets are not in media formulation, then add them
        mediaMets = [mediaMets; measuredMets(~iA)];
        fluxBounds = [fluxBounds; [fluxes(~iA) - abs(fluxTol*fluxes(~iA)), fluxes(~iA) + abs(fluxTol*fluxes(~iA))] ];
    end
end

if ~isEC
    
    modelStr = 'model';
    % set uptake fluxes for media mets
    [exchModel,unusedMets] = setExchangeBounds(model, mediaMets, fluxBounds(:,1), fluxBounds(:,2), true);

else
    
    modelStr = 'ecModel';
    [exchRxns, exchIndxs] = getExchangeRxns(model);
    
    % exclude protein pool exchange
    [~,prot_exch_indx] = ismember('prot_pool_exchange', model.rxns);
    if (prot_exch_indx == 0)
        error('Expected protein pool exchange reaction named "prot_pool_exchange" was not found.')
    else
        exchIndxs = setdiff(exchIndxs, prot_exch_indx, 'stable');
        exchRxns  = setdiff(exchRxns, 'prot_pool_exchange', 'stable');
    end
    
    % differentiate between uptake and production reactions
    uptkIndxs = exchIndxs(contains(exchRxns, '_REV'));
    prodIndxs = exchIndxs(~contains(exchRxns, '_REV'));
    
    % open all production reactions
    exchModel = setParam(model, 'ub', prodIndxs, maxBound);
    
    % close all uptakes
    exchModel = setParam(exchModel, 'ub', uptkIndxs, 0);
    
    % open uptake of media components one by one
    unusedMets = [];
    for i = 1:length(mediaMets)
        
        % get metabolite indx
        metIndx = getIndexes(model, strcat(mediaMets{i},'[e]'), 'metcomps');
        
        % get rxns for metabolite
        metRxns = find(model.S(metIndx,:));
        
        % get the uptake reaction for the metabolite
        metUptakeRxn = intersect(metRxns,uptkIndxs);
        if isempty(metUptakeRxn)
            if ~strcmpi(medium, 'all') || ismember(mediaMets(i), measuredMets)
                unusedMets = [unusedMets; mediaMets(i)];
            end
        else
            % set metabolite uptake bounds
            if fluxBounds(i,1) < 0
                % maximum allowed uptake flux
                exchModel.ub(metUptakeRxn) = abs(fluxBounds(i,1));
                if fluxBounds(i,2) < 0
                    % minimum allowed uptake flux (for measured mets)
                    exchModel.lb(metUptakeRxn) = abs(fluxBounds(i,2));
                end
            else
                exchModel.ub(metUptakeRxn) = 0;
            end
            
            % set metabolite production bounds
            metProdRxn = intersect(metRxns, prodIndxs);
            if fluxBounds(i,2) > 0
                % maximum allowed production flux
                exchModel.ub(metProdRxn) = abs(fluxBounds(i,2));
                if fluxBounds(i,1) > 0
                    % minimum allowed production flux (for measured mets)
                    exchModel.lb(metProdRxn) = abs(fluxBounds(i,1));
                end
            else
                exchModel.ub(metProdRxn) = 0;
            end
        end
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
        disp(['Constrained ' modelStr ' "' exchModel.id '" is feasible'])
    else
        disp(['*** Constrained ' modelStr ' "' exchModel.id '" is INFEASIBLE ***'])
        exchModel = [];
    end
end