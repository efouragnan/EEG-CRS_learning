function [ npar ] = get_npar_ef( modelID)
% Lookup table to get number of free parameters per model
% EF 2023


%%%%%
if       strcmp(modelID,'ms_RL'),            npar = 3;
elseif   strcmp(modelID,'RLpred'),           npar = 3;
elseif   strcmp(modelID,'RLpred_FH'),        npar = 4;
elseif   strcmp(modelID,'RLpred_RW'),        npar = 4;
elseif   strcmp(modelID,'RLpred_DYNA'),      npar = 4;

end


end

