function [drxns]= update_reaction_information(old_model, drxns)

original_indices = 1:length(old_model.rxns);
removed_indices = original_indices(~ismember(old_model.rxns,drxns.rxns));
original_rules= old_model.rules;
original_rxnEC=old_model.rxnECNumbers;
original_names=old_model.rxnNames;
original_Notes= old_model.rxnNotes;

R1=(removed_indices)';
R2=sort(R1, "descend");

% to remove reactions associated with removed_indices from rules field 
for i=1:length(R2)
    idx= R2(i);
    original_rules(idx)=[];
    original_rxnEC(idx)=[];
    original_names(idx)=[];
    original_Notes(idx)=[];
end 
drxns.rules=original_rules; 
drxns.rxnECNumbers=original_rxnEC; 
drxns.rxnNames=original_names; 
drxns.rxnNotes=original_Notes; 

