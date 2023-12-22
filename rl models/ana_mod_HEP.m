% Fits RL models and does model comparison
% MK Wittmann, October 2018
% updated by Elsa FOuragnan 2023
%
%%
%== -I) Prepare workspace: ============================================================================================

clearvars
cd('/Users/efouragnan/Library/CloudStorage/OneDrive-UniversityofPlymouth/root/data/HEP_EEG/rev_office/OneDrive_2023-03-08/RL guassian prior/')
addpath('models'); 
addpath('tools');
setFigDefaults; 

%== 0) Load and organise data: ==========================================================================================
% load data:
load('sorted_s.mat')

% define data set(s) of interest:
expids = {'beha'};

% how to fit RL:
M.dofit     = 1;                                                                                                     % whether to fit or not                                                           
M.doMC      = 1;                                                                                                     % whether to do model comparison or not  
M.quickfit  = 1;                                                                                                     % whether to lower the convergence criterion for model fitting (makes fitting quicker) (1=quick fit)
M.omitBMS   = 1;                                                                                                     % omit bayesian model comparison if you don't have SPM instaslled
M.modid     = {'ms_RL_HEP_simple2' 'ms_RL_HEP_simple' 'ms_RL_HEP_dyna','ms_RL_HEP_rw'}; % complete list of models to fit


%== I) RUN MODELS: ======================================================================================================

for iexp = 1:numel(expids)
   if M.dofit == 0,  break; end
   cur_exp = expids{iexp};                                                   
   
   %%% EM fit %%%
   for im = 1:numel(M.modid)
      dotry=1;
      while 1==dotry
         try
            close all;
            s.(cur_exp).em = EMfit_ms_HEP(s.(cur_exp),M.modid{im},M.quickfit);dotry=0;
         catch
            dotry=1; disp('caught');
         end
      end
   end

   %%% calc BICint for EM fit
   for im = 1:numel(M.modid)
      s.(cur_exp).em.(M.modid{im}).fit.bicint =  cal_BICint_ms(s.(cur_exp).em, M.modid{im});     
   end   
end

%== II) COMPARE MODELS: ================================================================================================
e=1
for iexp = 1:numel(expids)
   if M.doMC~=1, break; end
   cur_exp = expids{iexp};
   EMmc_ms(s.(cur_exp),M.modid,M.omitBMS);  
end



















