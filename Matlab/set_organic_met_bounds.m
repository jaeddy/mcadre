function model = set_organic_met_bounds(model)

% Add a warning if the metFormulas field is empty
carbonMets = ~cellfun('isempty', regexp(model.metFormulas, 'C'));
hydrogenMets = ~cellfun('isempty', regexp(model.metFormulas, 'H'));
is_organic = carbonMets & hydrogenMets;
organicMets = model.mets(is_organic);

organicRxns = findRxnsFromMets(model, organicMets);
organicExRxns = intersect(organicRxns, exRxns);
organicExRxns = [organicExRxns; ...
    'EX_Rtotal(e)'; 'EX_Rtotal2(e)'; 'EX_Rtotal3(e)'; 'EX_Tyr_ggn(e)'; ...
    'UP_Tyr_ggn[c]'];
model = changeRxnBounds(model, organicExRxns, 0, 'l');

model = changeRxnBounds(model, 'EX_glc(e)', -5, 'l');
model = changeRxnBounds(model, 'EX_co2(e)', -1000, 'l');
