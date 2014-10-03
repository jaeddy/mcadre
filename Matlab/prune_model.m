NC_removed = 0; C_removed = 0;
cTime = zeros(3000, 1); cRes = zeros(3000, 1);
count = 1;

while numel(P) %&& count < 3 % was just for testing, removed BH
    display(['Reaction no. ', num2str(count)])
    r = P(1);
    modelR = removeRxns(PM, r);

    % First check precursor production; if this test fails, no need to
    % check model consistency with FVA (time-saving step)
    [rStatus, time] = check_model_function(modelR, ...
        'biomass', biomassRxn,...
        'media', mediaDef);
    if rStatus

        % Check for inactive reactions after removal of r
        [inactive_G, time, result] = check_model_consistency(PM, r, C, 1);
        cTime(count) = time; cRes(count) = result;
        if ~numel(inactive_G)
            save mCADRE_error
            break;
        end

        % However, don't remove biomass or exchange reactions, even if
        % they're inactive
        inactive_G = setdiff(inactive_G, {'biomass'});

        inactive_G = setdiff(inactive_G, PM.rxns(findExcRxns(PM)));

        inactive_C = intersect(inactive_G, C);
        inactive_NC = setdiff(inactive_G, inactive_C);

        % Remove reactions with zero expression (previously penalized in
        % rankReactions) and corresponding inactive core reactions, only if
        % sufficiently more non-core reactions are removed
        if ismember(r, Z)

            % Check model function with all inactive reactions removed
            modelTmp = removeRxns(PM, inactive_G);
            tmpStatus = check_model_function(modelTmp, ...
                'biomass', biomassRxn,...
                'media', mediaDef);

            if numel(inactive_C) / numel(inactive_NC) <= eta && tmpStatus
                R_P = setdiff(R_P,inactive_G);
                PM = removeRxns(PM,inactive_G);
                P(ismember(P,inactive_G)) = [];
                NC_removed = NC_removed + numel(inactive_NC);
                C_removed = C_removed + numel(inactive_C);
                num_removed = NC_removed + C_removed;
                display('Removed all inactive reactions')

                % result = -1.x indicates that reaction r had zero
                % expression evidence and was removed along with any
                % consequently inactivated reactions; x = 1 or 2 depending
                % on whether removal of r created dead-end metabolites
                result = -1 - result / 10;
            else
                % Note: no reactions (core or otherwise) are actually
                % removed in this step, but it is necessary to update the
                % total number of removed reactions to avoid errors below
                num_removed = NC_removed + C_removed;
                P(1) = [];

                % result = 1.x.y indicates that no reactions were removed
                % because removal of r either led to a ratio of inactivated
                % core vs. non-core reactions above the specified threshold
                % eta (y = 1) or the removal of r and consequently
                % inactivated reactions prevented production of required
                % metabolites (y = 0)
                result = 1 + result / 10 + tmpStatus / 100;
            end

        % If reaction has expression evidence, only attempt to remove
        % inactive non-core reactions
        else
            % Check model function with non-core inactive reactions removed
            modelTmp = removeRxns(PM, inactive_NC);
            tmpStatus = check_model_function(modelTmp, ...
                'biomass', biomassRxn,...
                'media', mediaDef);

            if numel(inactive_C) == 0 && tmpStatus
                R_P = setdiff(R_P, inactive_NC);
                PM = removeRxns(PM, inactive_NC);
                P(ismember(P, inactive_NC)) = [];
                NC_removed = NC_removed + numel(inactive_NC);
                num_removed = NC_removed + C_removed;
                display('Removed non-core inactive reactions')

                % result = -2.x indicates that reaction r had expression
                % evidence and was removed along with (only) non-core
                % inactivated reactions; x = 1 or 2 depending on whether
                % removal of r created dead-end metabolites
                result = -2 - result / 10;
            else
                num_removed = NC_removed + C_removed;
                P(1) = [];

                % result = 2.x.y indicates that no reactions were removed
                % because removal of r either led to inactivated core
                % reactions (y = 1) or the removal of r and consequently
                % inactivated reactions prevented production of required
                % metabolites (y = 0)
                result = 2 + result / 10 + tmpStatus / 100;
            end
        end
    else
        num_removed = NC_removed + C_removed;
        P(1) = [];

        % result = 3 indicates that no reactions were removed because
        % removal of r by itself prevented production of required
        % metabolites
        result = 3;
    end
    cTime(count) = time; cRes(count) = result;
    count = count + 1;
    display(sprintf(['Num. removed: ', num2str(num_removed), ...
        ' (', num2str(C_removed), ' core, ', ...
        num2str(NC_removed), ' non-core); ', ...
        'Num. remaining: ', num2str(numel(P)), '\n']))
end
cTime(count:end) = []; cRes(count:end) = [];
