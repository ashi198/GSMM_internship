%% Flux balance analysis workflow 
% in order to perform FBA, uptake reactions corresponding to the metabolite
% of interests are required. The following workflow finds that. 
%% 
% initialize cobra toolbox 
initCobraToolbox(false)
solverOK=changeCobraSolver('gurobi','LP');

%% upload models
model_names= {'1A01_Kbase.sbml'; 'ModelSeed_1A01.sbml';....
    'carveMe_gapfilled.xml';'Raven_gapfilled.xml'; 'modelSeed_gapfilled.xml'}; 
for i=1:length(model_names)
    models{i, 1}= readCbModel(model_names{i}); % load all models in one cell array 
end

models_original=models;
model_raven=importModel('Raven_gapfilled.xml');
models{4}.c=model_raven.c;
models_original{4}.c= model_raven.c;

%% Find the index of the metabolite of interest

% add the names of all substrates
metabolites_of_interest={'Acetate'; 'glucose';'fructose';'glycerol'; 'galactose'; 'glucosamine'; 'alanine'; 'mannose'; 'Glycine'; 'Histidine'; 'Citrate'; 'Arginine';'Glutamine'; 'Lactate'; 'proline'; 'Lactose'; 'maltose'; 'mannitol'; 'succinate'; 'taurine'; 'cellobiose'; 'ribose'; 'xylose';....
    'Arabinose'; 'Sucrose'; 'Urea'; 'Valine'; 'Acetaldehyde'; 'Ethanol'; 'Adenine'; 'Sorbitol'; 'melibiose'; 'threonine'; 'Oxaloacetate'; 'Propionate'; 'Pyruvate'; 'Beta-alanine'; 'cysteine'; 'Galactitol'; 'lysine'; 'Leucine'; 'Methionine';....
    'tyrosine'; 'methanol'; 'O2'};

% findMetsFromList retrieve the indices and names of metabolites that have
% a name similar to the metabolite of interest. "Metabolites" and
% "metabolite_names" will produce a list of all metabolite that contains the
% name of the metabolite z.b the metabolites file for "acetate" will give
% results such as "'Oxaloacetate','Acetate', '4_Maleylacetoacetate', '2_Hydroxyphenylacetate'". 
% Strict_metabolites and strict_metabolite_names
% produce a list of metabolites that matches exactly with the query
% metabolite. 

%metabolite_for_uptake should be updated such that only one entry per
%metabolite is present in every cell. 

[metabolites, metabolite_names, strict_metabolites, strict_metabolite_names]= findMetsFromList(models, metabolites_of_interest);
metabolites_for_uptake=strict_metabolites; %% select all the appropriate metabolites before proceeding further

%% findUptakeRxns will produce a list of the uptake reactions present in all models  
[uptake_rxns]= findUptakeRxns(models); 

%% find out all uptake reactions that are associated with only the metabolites of interest 

matching_uptake_reactions=cell(length(metabolites_for_uptake), length(models)); 

for i=1:length(models)
  % have all metabolites for one models sequentially
  meta_temp=cell(length(metabolites_for_uptake), 1);

   for s=1:length(metabolites_for_uptake)
    meta_temp{s}=metabolites_for_uptake{s, i}; 
   end 

   temp_uptake_reactions= uptake_rxns{i}; %put uptake rxns from only one model 
 
% remove extra characters from the metabolite name
   for l=1:length(meta_temp)
    meta_temp = meta_temp(~cellfun('isempty', meta_temp));
    meta_temp{l} = regexprep(meta_temp{l}, '\[c0\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[e0\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_c\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_e\]$', ''); 
    meta_temp{l} = regexprep(meta_temp{l}, '\[c\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[e\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[p\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_p\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[c\]$', '');
   end

% pick all the uptake reactions associated with desired metabolites
    for j=1:length(meta_temp)
       idx=meta_temp{j};
       matching_reactions= {}; 
     for k=1:length(temp_uptake_reactions)
       if contains(lower(temp_uptake_reactions{k}), lower(idx))
            matching_reactions{end+1}= temp_uptake_reactions{k};
       end
     end
    matching_uptake_reactions{j, i} = matching_reactions; % this is the array with all the matching reactions
   end
end

%% Flux balance analysis
models=models_original;
biomass_yield=cell(length(matching_uptake_reactions), length(models)); 

%% adding biomass reaction for Raven
models{4} = addReaction(models{4},'bio1','metaboliteList', {'C00197[c]', 'C00024[c]', 'C00002[c]', 'C00279[c]', 'C00085[c]', 'C00118[c]', ...
               'C00092[c]', 'C00064[c]', 'C00025[c]', 'C00001[c]', 'C00004[c]', 'C00006[c]', ...
               'C00036[c]', 'C00074[c]', 'C00022[c]', 'C00117[c]', 'C00008[c]', 'C00026[c]', ...
               'C00010[c]', 'C00080[c]', 'C00003[c]', 'C00005[c]', 'C00013[c]'},'stoichCoeffList',[-1.496, -3.7478, -59.81, -0.361, 1, 0.129, 0.205, 0.2557, 4.9414, 59.81, ...
               -3.547, -13.0279, -1.7867, -0.5191, -2.8328, -0.8977, 59.81, 4.1182, 3.7478, 59.81, ...
               3.547, 13.0279, 59.81], 'reversible',false);

%set biomass reactions 
biomass_reactions={'bio1';'bio1';'Growth';'bio1'; 'bio1'};
%% peforming FBA 
%change [] to 0
   for i=1:length(models)
    for j = 1:(length(matching_uptake_reactions))
        %models{i}=models_original{i};
        if matching_uptake_reactions{j, i} == '0'
        biomass_yield{j, i}= 'NA';
        else 
        models{i} = changeRxnBounds(models{i}, matching_uptake_reactions{j, i}, -10, 'l'); % Set lower bound to -10
        models{i}=changeRxnBounds(models{i}, matching_uptake_reactions{45, i}, -1000, 'l'); %set to index of oxygen 
        models{i} = changeObjective(models{i}, biomass_reactions{i}); % Set objective to biomass
        sol = optimizeCbModel(models{i}, 'max'); % Perform FBA
        fprintf('Biomass yield on %s: %.6f\n', matching_uptake_reactions{j, i}, sol.f);
        biomass_yield{j, i}=sol.f;
        end
    end
  end 
 %% map a heat plot 
names=model_names;
temp_1 = zeros(size(biomass_yield));  % Initialize temp_1 with zeros

for i = 1:length(models)
    for j = 1:length(biomass_yield)
        if isempty(biomass_yield{j, i}) || strcmp(biomass_yield{j, i}, 'NA')  % Check for empty or 'NA' values
            temp_1(j, i) = 0;
        else
            temp_1(j, i) = biomass_yield{j, i};
        end
    end
end
h = heatmap(names,metabolites_of_interest,temp_1);
h.Title = 'Biomass yields on substrates';
h.XLabel = 'Organisms';
h.YLabel = 'Substrates';

%% check for accuracy

% create a confusion matrix
temp_bio_cell= zeros(length(biomass_yield), 1);
temp_bio = zeros(length(biomass_yield), 1);
bio_matrix= zeros(length(biomass_yield), length(models));

for i = 1: length(models)
     temp_bio_cell = temp_1(:, i);
    for k = 1:length(temp_bio_cell)
        temp_bio= temp_bio_cell(k);
        if temp_bio > 0
        bio_matrix(k, i)= 1;
        else
        bio_matrix(k, i)= 0; 
        end 
    end 
end 

% added entries for growth for manually curated model 
bio_matrix (:, 6)= [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0];
curated_model = bio_matrix(:, 6);

% Predicted models
predicted_models = bio_matrix(:, 1:5);

% Initialize variables
num_models = size(predicted_models, 2);
tp = zeros(1, num_models);
fp = zeros(1, num_models);
fn = zeros(1, num_models);
precision = zeros(1, num_models);
recall = zeros(1, num_models);
f1_score = zeros(1, num_models);

% Calculate TP, FP, FN for each model
for i = 1:num_models
    tp(i) = sum(predicted_models(:, i) & curated_model);
    fp(i) = sum(predicted_models(:, i) & ~curated_model);
    fn(i) = sum(~predicted_models(:, i) & curated_model);
    
    % Calculate precision, recall, and F1 score
    precision(i) = tp(i) / (tp(i) + fp(i));
    recall(i) = tp(i) / (tp(i) + fn(i));
    f1_score(i) = 2 * (precision(i) * recall(i)) / (precision(i) + recall(i));
end