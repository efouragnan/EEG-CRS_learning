function [fval,fit] = mod_ms_RL_HEP_dyna(behavData,q, doprior,dofit,varargin)
% runs standard RL model
% MK Wittmann, Oct 2018
%
% INPUT:    - behavData: behavioural input file
% OUTPUT:   - fval and fitted variables
%
%
%
%%
% -------------------------------------------------------------------------------------
% 1 ) Define free parameters
% -------------------------------------------------------------------------------------

if nargin > 4
    prior      = varargin{1};
end

qt = norm2par('ms_RL_HEP_dyna',q); % transform parameters from gaussian space to model space


% Define free parameters and set unused ones to zero
lrate        = qt(1);
beta         = qt(2);
choice_corr  = qt(3);
gamma        = qt(4);

% -------------------------------------------------------------------------------------
% 2-4) Middle code is the same for all models
% -------------------------------------------------------------------------------------

% Body of RL model that can run CL-trace, CS-trace and RT-trace. Compatible with all RL models.
% Ran directly from modms_RL_cl_cs_rt etc scripts
% MKW, Oct 2018
%%

% -------------------------------------------------------------------------------------
% 1 ) Define stuff, get variables
% -------------------------------------------------------------------------------------

% define behaviour on which model is fitted:
rt                = behavData.outcome;
offers            = behavData.offers;                                         % column 1 refers to right option, columm 2 refers to left option
choice            = behavData.choice_id;
choice(choice==2)=-1;
for it=1:length(rt)
    if choice(it)==1 && rt(it)==1
        rt(it)=1;
    elseif choice(it)==-1 && rt(it)==1
        rt(it)=-1;
    elseif choice(it)==-1 && rt(it)==0
        rt(it)=1;
    elseif choice(it)==1 && rt(it)==0
        rt(it)=-1;
    end
end
stimuli           = behavData.stimuli;
cue               = sum(stimuli,2)-1;
slope(1)=0;
a_t(1) = 0;

% prepare variables to collect information of interest
num_opt     = 2;
allvals     = nan(numel(rt),num_opt);
ChoiceProb  = nan(numel(rt),1);
stim1_val   = nan(numel(rt),1);
stim2_val   = nan(numel(rt),1);
pe          = nan(numel(rt),1);
pred_var    = nan(numel(rt),1);
prevch      = nan(numel(rt),1);
prevch(1)   = 0;

% define starting values:
allvals(1,:)    = zeros(1,num_opt);                                                    % order: [option1 option2 option3] % was 0 before sep18

%%
for it = 1:numel(rt)

    if it>1
        prevch(it) = choice(it-1);
    end

    % -------------------------------------------------------------------------------------
    % 2 ) Learning model:
    % -------------------------------------------------------------------------------------
    if cue(it)==1
        predch(it) = 0.5*allvals(it,1) + 0.5*allvals(it,1);
        pe(it) = rt(it) - predch(it);
        if it==1
            peabs(it) = abs(pe(it));
            a_t(it) = lrate;
        else
            peabs(it) = peabs(it-1)*(1-lrate) + abs(pe(it))*lrate;
            slope(it) = (peabs(it)-peabs(it-1)) / ((peabs(it)+peabs(it-1))/2);
            fm = sign(slope(it))*(1-exp(-(slope(it)/gamma)^2));
            if isnan(fm), fm=0; end;

            if slope(it)>=0, a_t(it) = a_t(it-1) + fm*(1-a_t(it-1));
            else a_t(it) = a_t(it-1) + fm*a_t(it-1);
            end
        end
        allvals(it+1,1) = allvals(it,1) + a_t(it)*pe(it);
        allvals(it+1,2) = allvals(it,2);
    elseif cue(it)==2
        predch(it) = 0.5*allvals(it,1) + 0.5*allvals(it,2);
        pe(it) = rt(it) - predch(it);
        if it==1
            peabs(it) = abs(pe(it));
            a_t(it) = lrate;
        else
            peabs(it) = peabs(it-1)*(1-lrate) + abs(pe(it))*lrate;
            slope(it) = (peabs(it)-peabs(it-1)) / ((peabs(it)+peabs(it-1))/2);
            fm = sign(slope(it))*(1-exp(-(slope(it)/gamma)^2));
            if isnan(fm), fm=0; end;
            if slope(it)>=0, a_t(it) = a_t(it-1) + fm*(1-a_t(it-1)); else a_t(it) = a_t(it-1) + fm*a_t(it-1); end
        end
        allvals(it+1,1) = allvals(it,1) + a_t(it)*pe(it);
        allvals(it+1,2) = allvals(it,2) + a_t(it)*pe(it);
    elseif cue(it)==3
        predch(it) = 0.5*allvals(it,2) + 0.5*allvals(it,2);
        pe(it) = rt(it) - predch(it);


        if it==1
            peabs(it) = abs(pe(it));
            a_t(it) = lrate;
        else
            peabs(it) = peabs(it-1)*(1-lrate) + abs(pe(it))*lrate;
            slope(it) = (peabs(it)-peabs(it-1)) / ((peabs(it)+peabs(it-1))/2);
            fm = sign(slope(it))*(1-exp(-(slope(it)/gamma)^2));
            if isnan(fm), fm=0; end;
            if slope(it)>=0, a_t(it) = a_t(it-1) + fm*(1-a_t(it-1)); else a_t(it) = a_t(it-1) + fm*a_t(it-1); end
        end
        allvals(it+1,2) = allvals(it,2) + a_t(it)*pe(it);
        allvals(it+1,1) = allvals(it,1);
    end

    % -------------------------------------------------------------------------------------
    % 3 ) Observation model:
    % -------------------------------------------------------------------------------------

    dv(it) = beta*(predch(it) - .5) + choice_corr*prevch(it);
    ChoiceProb(it) = 1 / (1 + exp(-dv(it)));

end

%%
allvals(end,:)=[];


% -------------------------------------------------------------------------------------
% 4 ) Calculate model fit:
% -------------------------------------------------------------------------------------

% find when 1 was actually choosen:
choice_bin = (choice == 1);
nll =-nansum(choice_bin.*log(ChoiceProb));   % the thing to minimize

if doprior == 0                                                               % NLL fit
    fval = nll;
elseif doprior == 1                                                           % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign
    fval = -(-nll + prior.logpdf(q));
end

if sum(isnan(ChoiceProb))>0,
    disp('ERROR');
    keyboard;
    return;
end                    % error check


%--------------------------------------------------

if isinf(fval)
    e=1;
else
    e=0;
end

% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save:
% -------------------------------------------------------------------------------------

if dofit ==1
    vsum = sum(allvals,2);

    fit         = struct;
    fit.xnames  = {'alpha'; 'beta';'Choice_corr'};

    fit.mat    = [allvals  ChoiceProb dv' ...
        predch' pe];
    fit.names  = {'V1';'V2'; 'ChoiceProb' ; 'DV'; ...
        'PredVariable';'PE'};
end




end

