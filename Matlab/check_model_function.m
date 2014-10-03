function [requiredMetsStatus,time] = checkModelFunction(model, varargin)

% Inputs
% - model:
% - optional inputs:
%       'requiredMets',metList
%       'biomass',biomassRxn
%       'media',mediaDef

%% Parse input parameters

if numel(varargin)
    if rem(numel(varargin),2) == 0;
        options = {'requiredMets', 'biomass', 'media'};
        for i = 1:2:numel(varargin)
            argname = varargin{i}; argval = varargin{i + 1};
            option = find(strncmpi(argname,options,numel(argname)));
            if ~isempty(option)
                switch option
                    case 1
                        metList = argval;
                    case 2
                        biomassRxn = argval;
                    case 3
                        mediaDef = argval;
                end
            else error('Unknown option %s.', argname)
            end
        end
    else error('Incorrect number of input arguments to function %s.', ...
            mfilename);
    end
end

%%
t0 = clock;

% Identify exchange reactions in the model
indOutRxns = sum(model.S > 0) > 0 & sum(model.S < 0) == 0;
indInRxns = sum(model.S < 0) > 0 & sum(model.S > 0) == 0;
indExRxns = indOutRxns | indInRxns;
exRxns = model.rxns(indExRxns);

if exist('mediaDef')
    model = set_media_ex_bounds(model);
else
    model = set_organic_met_bounds(model);
end

% Add demand reactions for required metabolites
if exist('metList')
    [model, requiredRxns] = specify_required_rxns(model, metList);
else requiredRxns = {};
end

if exist('biomassRxn')
    requiredRxns = [requiredRxns, biomassRxn];
end

inactiveRequired = check_rxn_flux(model, requiredRxns);

requiredMetsStatus = ~numel(inactiveRequired);
time = etime(clock, t0);
