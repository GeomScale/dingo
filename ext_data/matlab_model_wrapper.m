path_to_model = '';
model = load(strcat(path_to_model, '.mat'));

fldnames = fieldnames(model);

dingo_model = struct;
dingo_model.S = model.(fldnames{1}).S;
dingo_model.lb = model.(fldnames{1}).lb;
dingo_model.ub = model.(fldnames{1}).ub;
dingo_model.c = model.(fldnames{1}).c;
dingo_model.index_obj = find(dingo_model.c == 1);
dingo_model.rxns = model.(fldnames{1}).rxns;
dingo_model.mets = model.(fldnames{1}).mets;

save('dingo_model.mat', 'dingo_model')

