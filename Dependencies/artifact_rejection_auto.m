function [data_clean_artifacts] = artifact_rejection_auto(fname, data_seg, trl_visual_only,rd)

% trl_visual_only= trl_visual_only(1:10,:);
% trl_short = trl(1:10,:);

%%%%%%%%%%%%%%%%%% detection of muscle artifacts in taks trials and rest trials %%%%%%%%%%%%%%%%%
% http://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection
cfg            = [];
cfg.dataset = fname;
cfg.continuous = 'no'; % changed from yes
% channel selection excluding ECG channel, cutoff and padding
cfg.artfctdef.zvalue.channel = rd.listChan([1:30,32:end]);
cfg.artfctdef.zvalue.cutoff      = 20; % changed from 8 
cfg.artfctdef.zvalue.trlpadding  = 0;% 0
cfg.artfctdef.zvalue.fltpadding  = 0;% 0
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

%the maximum value of the artifact within its range will be found and stored;
cfg.artfctdef.zvalue.artfctpeak = 'no';

% make the process interactive GUI WITH FEEDBACK
cfg.artfctdef.zvalue.interactive = 'no';

% Detecting artifact from trl time-locked to visual onset 
cfg.trl        = trl_visual_only;
[cfg, artifact_muscle] = ft_artifact_zvalue(cfg); %visual_muscle = cfg.artfctdef.zvalue.peaks;
save ([rd.dirSav,['artifact_muscle_' num2str(rd.subj) '.mat']],'artifact_muscle');%save ([rd.dirSav,['visual_muscle_' num2str(rd.subj) '.mat']],'visual_muscle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Remove muscle artifcats from the data %%%%%%%%%%%%%%%%%%%%%%%%
% not removing jumps (typically from MEG data), but check if this is necessary 
cfg=[];
cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
cfg.artfctdef.feedback = 'yes';
cfg.artfctdef.muscle.artifact = artifact_muscle;
data_no_artifacts = ft_rejectartifact(cfg,data_seg);
data_no_artifacts.cfg.previous =[]; % this reduces the file size
data_clean_artifacts = data_no_artifacts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% cfg.trl = trl_R_feed; 
% [cfg, artifact_muscle] = ft_artifact_zvalue(cfg); R_muscle = cfg.artfctdef.zvalue.peaks;
% art_R = artifact_muscle; save ([rd.dirData,['art_R_' num2str(rd.subj) '.mat']],'art_R');save ([rd.dirData,['R_muscle_' num2str(rd.subj) '.mat']],'R_muscle');
% 
% 
% all_artifact_muscle = [art_visual;art_R]; all_artifact_muscle = sortrows(all_artifact_muscle,1);
