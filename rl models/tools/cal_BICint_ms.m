function [ bicint ] = cal_BICint_ms(mroot, modelID,Nsample )
% calculates BICint. 
% MK Wittmann, 2017
% 
% INPUT:    - mroot: model root with all relevant information
%           - modelID: the name of the model
% OUTPUT:   - bicint
%
%%


% define settings
if nargin <3, Nsample     = 2000; end                                                        % nr of times drawn from population distribution

% dont do prior, just get the NLL
fit.dofit   = 0;                                                                             % can be zero, because no fitted parameters needed, just nll
fit.doprior = 0;
fit.objfunc = str2func(['mod_' modelID]); 
fit.npar    = size(mroot.(modelID).q,2);
fit.ntrials = mroot.(modelID).ntrials;
fit.beh     = mroot.(modelID).behaviour;

% info for normpdf, and flip if it is the wrong orientation
mu          = mroot.(modelID).gauss.mu;            if size(mu,2)>size(mu,1),                 mu = mu'; end
sigmasqrt   = sqrt(mroot.(modelID).gauss.sigma);   if size(sigmasqrt,2)>size(sigmasqrt,1),   sigmasqrt = sigmasqrt'; end                       

% collect integrated nll
iLog        = nan(numel(fit.beh),1);                                                         % proxy for integrated log


%% start computing

%%% 1)get integrated nll by sampling nll from group gaussian
fprintf([modelID ' - BICint: ']);
for is = 1:numel(fit.beh)
   
   subnll = nan(1,Nsample);
   Gsamples    = normrnd(repmat(mu,1,Nsample),repmat(sigmasqrt,1,Nsample));                     % samples from gaussian distribution found during EM; draw anew for each subject

   % for each subject, get NLL for input params from gaussian
   for k=1:Nsample
      subnll(k)  = fit.objfunc(fit.beh{is},Gsamples(:,k),fit.doprior,fit.dofit); 
   end
   fprintf([ num2str(is) ',']);
   iLog(is) = log(sum(exp(-subnll))/Nsample);
end
     

%%% 2) Compute BICint

bicint  = -2*sum(iLog)   + fit.npar*log(sum(fit.ntrials));                     % integrated BIC 
fprintf(['\n']);






end

