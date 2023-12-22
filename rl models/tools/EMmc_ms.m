function [ out ] = EMmc_ms(rootfile,modnames,omitBMS)
% Model comparison for EM fitted models
% MK Wittmann, Nov 2017
%
% INPUT:        - rootfile: root of a specific experiment
%               - modnames: which model to looks at
%               - omitBMS: you can omit plotting the bayesian model comparison, because it requres you to install spm
%
%%


expname = rootfile.expname;


%%%% get info
for imod=1:numel(modnames)
   lme(:,imod)    = rootfile.em.(modnames{imod}).fit.lme;
   bicint(imod)   = rootfile.em.(modnames{imod}).fit.bicint;
end


% I. Plot LME and BIC: --------------------------------------------------
% 1) LME sum
h=figure('name', expname);
suptitle('Bayesian Model Comparison');

subplot(2,2,1);   bar(sum(lme)); xtickrot = 25;
set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
ylabel('Summed log evidence (more is better)','FontWeight','bold');

% 2) BICint
subplot(2,2,2);
bar(bicint);
set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames);
set(gca,'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
ylabel('BICint (less is better)','FontWeight','bold');
%------------------------------------------------------------------------


% II Calculate exceedance probability and compare log model evidence ----
% Compare LME of two models:
compM   = {'ms_RL','ms_RL_cl'};                                    % models to compare specifically
subplot(2,2,3);
for imod=1:numel(compM)
   lmepair(:,imod)    = rootfile.em.(compM{imod}).fit.lme;
   bicint(imod)   = rootfile.em.(compM{imod}).fit.bicint;
end
lmediff = lmepair(:,2) - lmepair(:,1); 
barh(1:numel(lmediff), sort(lmediff),'Facecolor',[237 177 32]./255); hold all;
set(gca,'ytick',1:numel(lmediff));
ylabel('Sessions (sorted)');
xlabel([compM{1} ' < ' compM{2} ]);
title('Log Evidence Differences / Log Bayes factor');


% this is coming from SPM you need to have it somewhere in your paths
if omitBMS == 0
   [~,~,BMS.xp] = spm_BMS(lme);                                                   % log model evidence
   subplot(2,2,4)
   bar(BMS.xp);
   set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
   rl1 = refline(0,.95);     set(rl1,'linestyle','--','Color','r');
   ylabel('Exceedance Probability','FontWeight','bold');
end
setfp(gcf)
figname=['EM_BMC_' expname '_'  date];
saveas(gcf,figname);



