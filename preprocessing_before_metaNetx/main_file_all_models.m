% this file just concerns with running the loading models and running the
% pipeline 
%% Initialize cobratoolbox and set solver
initCobraToolbox(false);
solverOK=changeCobraSolver('gurobi','LP');

%% Raven 
RE1A=importModel('Raven_edited_main.xml');
% remove bad reactions
[updated_RE1A, deleted_reactions_RE]= removeBadRxns(RE1A, 3);
% Deleted orphan gene VEJY3_oo520 from gene list to avoid rxnGene Mat
% mismatch 
writeCbModel(updated_RE1A, 'RE1A_from_pipeline.xml', 'format', 'sbml');

%% CarveMe 
CE1A=readCbModel('1A01_carveMe.xml');
[updated_CE1A, deleted_reactions_CE] = run_pipeline_for_rest(CE1A);
new_CE1A= update_reaction_information(CE1A, updated_CE1A);
writeCbModel(new_CE1A, 'CE1A_from_pipeline.xml', 'format', 'sbml');
%% Kbase 
KE1A=readCbModel('1A01_Kbase.sbml'); %%Warning only use readCbModel to import Kbase and ModelSEED
[updated_KE1A, deleted_reactions_KE] = run_pipeline_for_SEED(KE1A);
writeCbModel(updated_KE1A, 'KE1A_from_pipeline.xml', 'format', 'sbml');

%% ModelSEED
ME1A=readCbModel('ModelSeed_1A01.sbml'); %%Warning only use readCbModel to import Kbase and ModelSEEDK
[updated_ME1A, deleted_reactions_ME] = run_pipeline_for_SEED(ME1A);
writeCbModel(updated_ME1A, 'ME1A_from_pipeline.xml', 'format', 'sbml');

% The resulting models from modelSEED and Kbase needs to be validated by SBML validator. It will show
% a bunch of errors regarding "null" value in chemical species for "RNA
% transcription", "DNA replication", "protein biosynthesis" and "biomass".
% Enter the appropriate value for this and remove the
% <geneProduct></geneProducts> tags if GrRules do not exist in the model. 
% for Kbase--> the following was changed name="DNA_replication_c0"
% fbc:chemicalFormula="null"; RNA_transcription_c0;Protein_biosynthesis_c0;
% MPT_Synthases_c0;Acceptor_c0; Biomass_c0_b
