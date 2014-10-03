function deadend_C = check_core_deadends(model, C)

deadEndMets = [];
for metNo = 1:length(model.mets)
    bothRxns = find(model.S(metNo, :) ~= 0 & model.rev(:)');
    prodRxns = union(find(model.S(metNo,:) > 0), bothRxns);
    consRxns = union(find(model.S(metNo,:) < 0), bothRxns);
    if (isempty(consRxns))
        deadEndMets(end + 1) = metNo;
    elseif (isempty(prodRxns))
        deadEndMets(end + 1) = metNo;
    else
        if (length(prodRxns) == 1 & length(consRxns) == 1 & prodRxns == consRxns)
            deadEndMets(end + 1) = metNo;
        end
    end
end

deadEndMets = detectDeadEnds_fast(model);

deadEnd = sum(full(model.S(deadEndMets, :) ~= 0), 1) > 0;
deadEndRxns = model.rxns(deadEnd);
deadEnd_C = intersect(rxnList, deadEndRxns);
