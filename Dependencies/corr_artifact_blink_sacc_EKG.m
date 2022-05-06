function [all_artifacts] = corr_artifact_blink_sacc_EKG(cfg)

% unpack cfg
compdat = cfg.compdat; %components data
sigdat = cfg.sigdat;  % EEG data
sig_index = cfg.sig_index; nindex = length(sig_index);% number/index of EOG channel
p_thresh = cfg.p_thresh; % pvalue thres .05
r_thresh = cfg.r_thresh; % r value thres .1

blink =[]; saccade =[]; ECG =[]; all_artifacts =[];

for s = 1:nindex
      indx = sig_index(s);
    % create matrices for input data
    compmat = zeros(numel(compdat.trial),size(compdat.trial{1},1),size(compdat.trial{1},2));
    sigmat = zeros(numel(sigdat.trial),size(sigdat.trial{1},2));

    for trl = 1:numel(compdat.trial)
                compmat(trl,:,:) = compdat.trial{trl};
                sigmat(trl,:) = sigdat.trial{trl}(indx,:);
    end
       
   % trials to include in correlation
     ntrials = numel(compdat.trial);
        
    % correlate components with signal
    r = []; p = [];
    for i = 1:size(compdat.trial{1},1)
        [r(i),p(i)] = corr(reshape(squeeze(compmat(1:ntrials,i,:))',1,numel(compmat(1:ntrials,i,:)))',reshape(sigmat(1:ntrials,:)',1,numel(sigmat(1:ntrials,:)))');
    end
    
    % find highly correlated components
    sig_corr = find(p < p_thresh);
    corr_comp = find(abs(r) > r_thresh);
    
    % take comps above threshold
    if ~isempty(intersect(sig_corr,corr_comp))
        artefact_comp = intersect(sig_corr,corr_comp);
    else % take highest one
        artefact_comp = find(r == max(r));
    end
    
    % if s is 64 which is th VEOG...
    if s == 1, blink = artefact_comp; elseif s ==2, saccade = artefact_comp; else ECG = artefact_comp; end
    
    clearvars -except ECG blink cfg compdat nindex p_thresh r_thresh saccade sig_index sigdat artefact_comp
end


all_artifacts = union([blink, saccade], ECG);
end