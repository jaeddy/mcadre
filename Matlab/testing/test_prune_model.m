
load('HR1_CbModel');
load('precursorMets');
changeCobraSolver('glpk');

%% Test update_model

display('## Testing prune_model...');
try
    [PM, cRes] = prune_model(GM, P, C, eta, precursorMets, method);
    display(['PASS...', ...
        'Function prune_model ran without error']);
catch err
    display(['FAIL...', ...
        'Function prune_model was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');