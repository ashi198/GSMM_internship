function [model_after_rules]= create_rules(old_model, new_model)

original_indices = 1:length(old_model.rxns);
removed_indices = original_indices(~ismember(old_model.rxns,new_model.rxns));
original_rules= new_model.rules;
R1=(removed_indices)';
R2=sort(R1, "descend");

% to remove reactions associated with removed_indices from rules field 
for i=1:length(R2)
    idx= R2(i);
    original_rules(idx)=[];
end 
new_model.rules=original_rules; 
model_after_rules=new_model; 