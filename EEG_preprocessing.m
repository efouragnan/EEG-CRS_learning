%%% Alejandra Sel 25/06/2018 [BH Update 2021] %%%%
%%% Pre-processing EEG data HEP reward learning task

close all
clear all


% Configure eeglab, ERPlab, fieldtrip and dependencies
% Update with your paths.
addpath('\eeglab2020');
addpath('\CARE-rCortex');
addpath('\ERPLAB8.02');
addpath('\Dependencies');
addpath('\fieldtrip');


surrogate = 0; %Set to 1 if creating surrogates for data validity tests.

% Define data pathway.

rd.dirProj = ('\Data\'); %Full Project directory.
rd.dirData =  fullfile(rd.dirProj); %Used in place of dirProj for loading data.
rd.dirSavP =  ('\Results\'); % Output project directory.

% Usable Participants
rd.listSubj = {'02','04','05','06','07','08','10','11','12','13','14','15','16','17','18','19','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36'};
rd.nSubj = length(rd.listSubj);
% some participants have channels that need to be repaired, notably 6 & 12.
rd.repchan = [ 0,   0 ,  0,   1,   0,   0  ,  0,   0,   1,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ,  0,  0,  0,   0,   0,   0,    0,   0,   0,   0,  0,   0,  0];




% Define list of channels.
rd.listChan = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F7' 'F8' 'T7' 'T8' 'P7' 'P8' 'Fz' 'Cz' 'Pz' 'FCz' 'FC1' 'FC2' 'CP1' 'CP2' 'FC5' 'FC6' 'CP5' 'CP6' 'TP9' 'TP10'...
    'EXT1' 'F1' 'F2' 'C1' 'C2' 'P1' 'P2' 'AF3' 'AF4' 'FC3' 'FC4' 'CP3' 'CP4' 'PO3' 'PO4' 'F5' 'F6' 'C5' 'C6' 'P5' 'P6' 'AF7' 'AF8' 'FT7' 'FT8' 'TP7' 'TP8' 'PO7' 'PO8' 'Fpz' 'CPz' 'POz' 'Oz'};
nchan = length(rd.listChan);




% Loop over participants
data = [];
fprintf('\n');
% 92 conditions in total - including correct and incorrect trials
ncond = 92;
allppts_all_conds_ave(1,ncond) = {0}; allppts_all_conds(1,ncond) = {0}; all_ppt_trlnumb_feed = []; all_ppt_trlnumb_visual = []; all_ppt_trlnumb(1,ncond) = 0;
allppts_HEP_corr_conds(1,4) ={0}; % 4 conditions for HEP correct trials (average before and after reversal and across colours/items)

dwnsmple = 1; %Downsample or not. This is the default. Everything is set up around downsampled data.


for s = 1:rd.nSubj
    rd.subj = rd.listSubj{s};
    fprintf('\n\t%s',rd.subj);
    mkdir(fullfile(rd.dirProj,['sub' rd.subj]))
    rd.dirData = fullfile(rd.dirProj,['sub' rd.subj],'/');
    
    if dwnsmple == 1
        mkdir(fullfile(rd.dirSavP,['sub' rd.subj],'/Downsampled'));
        rd.dirSav = fullfile(rd.dirSavP,['sub' rd.subj],'/Downsampled/');
    else
        rd.dirSav = fullfile(rd.dirSavP,['sub' rd.subj],'/FullData/');
    end
    
    
    % For later switch. Currently requires full run, as there is no loading
    % within each case, just continual running.
    % todo = {'preproc','artifact_corr','ica'};
    todo = {'preproc','artifact_corr','ica','conditions'};
    
    
    
    
    for t = 1:length(todo)
        r = char(todo(t));
        
        
        
        switch r
            case 'preproc'
                
                % Define trials - they are different in some participants but
                % all feedback conditions can appear, removed the if statement
                % as it was losing trials.
                
                
                % Conditions time-locked to feedback, two groups of participants
                
                rd.listCond_feedback_g1 = {' 11','12',' 13','14',' 21','22','23','24','31','32','33','34','51','52','53','54'};
                rd.listCond_feedback_g2 = {' 11','12',' 13','14',' 21','22','23','24','41','42','43','44','51','52','53','54'};
                rd.listCond_visual_g1 = {'61',' 62',' 63','64','65','66','67','68','71','72','73','74','75','76','77','78','81',' 82',' 83',' 84','85','86','87','88','101','102','103','104'};
                rd.listCond_visual_g2 = {'61',' 62',' 63','64','65','66','67','68','71','72','73','74','75','76','77','78','91',' 92',' 93',' 94','95','96','97','98','101','102','103','104'};
                
                
                rd.nCond = length(rd.listCond_feedback_g1);
                rd.nCond = length(rd.listCond_visual_g1);
                
                % Define data file.
                fname = fullfile(rd.dirData, ['Pilot_' rd.subj '.eeg']);
                
                % Read in the data, selecting active channels
                cfg = []; data = [];
                cfg.dataset = fname;
                cfg.channel = [1:30 32:64]; % keep active chanels, excluding 31(LM).
                data = ft_preprocessing(cfg);
                
                
                
                
                
                
                % Repairing bad channels in the best way that would work.
                if rd.repchan(s)
                    cfg =[]; cfg.neighbours = struct; cfg.trials = 'all';
                    if strcmp (rd.subj,'06')
                        cfg.neighbours(1).label = 'T8';
                        cfg.neighbours(1).neighblabel = {'FT8'; 'C6'; 'TP8'};
                        cfg.badchannel = {'T8'};
                        
                    end
                    if strcmp (rd.subj,'12')
                        cfg.neighbours(1).label = 'AF7';
                        cfg.neighbours(1).neighblabel = {'Fp1'; 'AF3'; 'F7'};
                        cfg.badchannel = {'AF7'};
                    end
                    [channelmatch, channelidx] = intersect(data.label, cfg.neighbours(1).neighblabel);
                    [badchannelmatch, badchannelidx] = intersect(data.label, cfg.badchannel);
                    load ('layout.mat')
                    [~, badchanidx] = intersect(lay.label, cfg.badchannel);
                    [~, goodchanidx] = intersect(lay.label, cfg.neighbours(1).neighblabel);
                    distance = sqrt(sum((lay.pos(goodchanidx, :) - repmat(lay.pos(badchanidx, :), numel(goodchanidx), 1)).^2, 2));
                    %weighting by distance
                    weight = distance./sum(distance);
                    data.trial{1}(badchannelidx,:) = sum(weight.*data.trial{1}(channelidx,:)) / sum(weight);
                    
                end
                
                
                
                % Bandpass filter.
                cfg = [];
                cfg.bpfilter = 'yes';
                cfg.bpfreq        = [0.5 40];
                data_filtered = ft_preprocessing(cfg, data);
                
                % Downsample after filtering.
                if dwnsmple == 1
                    
                    cfg = [];
                    cfg.resamplefs = 500;
                    data_filtered  = ft_resampledata(cfg, data_filtered);
                end
                
                % To remove hardcoded variables.
                samplefactor = cfg.resamplefs/1000;
                
                
                % Extract the raw ecg channel, Detect the R-waves with pan_tompkin
                
                cfg=[];
                cfg.channel = data_filtered.label(31); % EKG channel = 31
                ecg_struct = ft_preprocessing(cfg,data_filtered);
                ecg = ecg_struct.trial{1}';
                fs = data_filtered.fsample;
                gr = 0;
                [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs,gr);
                
                
                % surrogacy creation for control analysis. 
                if surrogate
                    qrs_i_raw = qrs_i_raw(1,:)+ (200 * samplefactor);
                end % add 200ms to R onset to create surrogates
                
                R_trl = [qrs_i_raw(2:end)' - (200*samplefactor), qrs_i_raw(2:end)' + (600*samplefactor)];
                R_trl(:,3) = -(200*samplefactor); % trl for R waves, -200ms for baseline, 600ms epoch
                
                
                
                % Trial definition feedback, according to group variable predictive face vs house
                
                % Required in epoch definition.
                cfg = [];
                cfg.dataset = fname;
                cfg.trialdef.pre = 0.9; % 900ms before feedback onset
                cfg.trialdef.post = 4; % 4secs after feedback onset
                %cfg.trialdef.eventtype = 'trigger';
                %cfg.trialdef.eventtype = 'Stimulus';
                cfg.trialdef.vispre = 0.2; % 200ms before stimuli onset
                cfg.trialdef.vispost = 0.5; % 500ms after stimuli onset
                cfg.trialdef.decpre = 0.1;
                cfg.trialdef.decpost = 1;
                
                
                
                % Function - define epochs based on triggers
                % This requires the events.mat structure as behavioural
                % data is used to correct for weird recording of trigger
                % codes.
                if dwnsmple == 1
                    [all_trl_visual, all_trl_decision, all_trls_feed, trl_visual_feed, trl_R_feed, group] = define_epochs_all(cfg,R_trl,rd,data_filtered);
                else
                    % code is tailored towards downsampling, but should
                    % still work.
                    [all_trls_feed, trl_visual_feed, trl_R_feed,  event_feedback, group, trlnumb_feedback] = define_epochs_all(cfg,R_trl,rd,data_filtered);
                end
                
                save ([rd.dirSav,['all_trls_feed_' num2str(rd.subj) '.mat']],'all_trls_feed');
                save ([rd.dirSav,['trl_visual_feed_' num2str(rd.subj) '.mat']],'trl_visual_feed');
                save ([rd.dirSav,['trl_R_feed_' num2str(rd.subj) '.mat']], 'trl_R_feed');
                save ([rd.dirSav,['all_trl_visual_' num2str(rd.subj) '.mat']],'all_trl_visual');
                save ([rd.dirSav,['all_trl_decision_' num2str(rd.subj) '.mat']],'all_trl_decision');
                
                
                % Merge all trials from face,house presentation and feedback
                all_trl_visual (:,6:8)= NaN;
                all_trl_visual = all_trl_visual(:,[1:4 6:7 5 8]);
                all_trl_decision(:,6:8) = NaN;
                all_trl_decision = all_trl_decision(:,[1:4 6:7 5 8]);
                trl_all = [all_trls_feed;all_trl_visual;all_trl_decision];
                trl_all = sortrows(trl_all,1);
                save ([rd.dirSav,['trl_all_' num2str(rd.subj) '.mat']],'trl_all');
                
                
                % Segmentation
                trl_all_basic = trl_all(:,[1:5,7]);
                cfg = [];
                cfg.trl = trl_all_basic;
                data_seg = ft_redefinetrial(cfg,data_filtered); %I cannot use definetrial because I input preprocess data
                
                
                
                
                % Some participants had their recordings abruptly ended,
                % this removes and NaN's recorded to make the data workable.
                if any([17,31] == str2double(rd.subj))
                    [data_seg, trl_visual_feed] = NaNremove(data_seg, trl_visual_feed);
                    save ([rd.dirSav,['data_seg_' num2str(rd.subj) '.mat']],'data_seg','-v7.3');
                    save ([rd.dirSav,['trl_visual_feed_' num2str(rd.subj) '.mat']],'trl_visual_feed');
                else
                    save ([rd.dirSav,['data_seg_' num2str(rd.subj) '.mat']],'data_seg','-v7.3');
                end
                
                
                
                % Trying to save RAM by removing unused things. Should
                % probably use clearvars -except
                clear rd.listCond_feedback_g1 rd.listCond_feedback_g2 rd.listCond_visual_g1 rd.listCond_visual_g2 rd.nCond rd.subj rd.dirData fname data_filtered ecg fs gr ecg_struct event_visual trlnumb_visual trl_all_basic qrs_amp_raw qrs_i_raw event_feedback event_visual all_trls_feed trlnumb_feedback trlnumb_visual R_trl ecg fs gr ecg_struct
                %trl_all
                
                
            case 'artifact_corr'
                
                % Load raw data, segmented data, trl time-locked to house/face and to feedback onset,
                % and trl time-locked to R during feedback,
                fname = fullfile(rd.dirData, ['Pilot_' rd.subj '.eeg']);
                
                % Uncomment if not starting from the beginning.
                %load ([rd.dirSav,['data_seg_' num2str(rd.subj) '.mat']],'data_seg');
                %load ([rd.dirSav,['trl_visual_feed_' num2str(rd.subj) '.mat']],'trl_visual_feed'); load ([rd.dirSav,['all_trl_visual_' num2str(rd.subj) '.mat']],'all_trl_visual');
                %load ([rd.dirSav,['trl_R_feed_' num2str(rd.subj) '.mat']],'trl_R_feed');
                trl_visual_feed(:,4:end) =[];
                all_trl_visual(:,4:end) =[];
                trl_R_feed(:,4:end) = [];
                trl_visual_only = [trl_visual_feed;all_trl_visual;trl_R_feed];
                trl_visual_only = sortrows(trl_visual_only,1);
                
                
                % Clean and save data - only clean segments time-locked to
                % feedback which includes trl_R_feed.
                % Z threshold is very generous.
                % Currently set at 20. This keeps trial losses to around
                % 20% max.
                [data_clean_artifacts] = artifact_rejection_auto(fname, data_seg, trl_visual_only,rd);
                
                save ([rd.dirSav,['data_clean_artifacts_' num2str(rd.subj) '.mat']],'data_clean_artifacts','-v7.3');
                
                
                clear Ntrlsrej Ntrlsreject trl_visual_feed all_trl_visual data_seg
                
            case 'ica' 
                
                % Will do both per trial and per heartbeat, including
                % average heartbeat.
                
                %endcond = {'PTrial','PHEP','PStim'};
                %endcond = {'PTrial', 'PHEP'};
                
                % We will be using the events structure here to add
                % behavioural and model data to the final structures.
                
                % Load here if beginning the process here.
                % load ([rd.dirSav,['data_clean_artifacts_' num2str(rd.subj) '.mat']],'data_clean_artifacts');
                load ([rd.dirData 'BEHAV\events.mat']) % Behavioural data and learning models.
               
                
                trig_feed = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44,51,52,53,54]; all_ind_trig_feed =[];
                for jj = 1:length(trig_feed)
                    trigfeed = trig_feed(jj);
                    indtrigfeed = find(data_clean_artifacts.trialinfo(:,1) == trigfeed);
                    all_ind_trig_feed =[all_ind_trig_feed;indtrigfeed];
                end
                all_ind_trig_feed = sortrows(all_ind_trig_feed);
                clear trig_feed jj trigfeed indtrigfeed
                
                
              
                sub_data_clean_artifacts = data_clean_artifacts;
                sub_data_clean_artifacts.trial = sub_data_clean_artifacts.trial(1,all_ind_trig_feed);
                sub_data_clean_artifacts.time = sub_data_clean_artifacts.time(1,all_ind_trig_feed);
                sub_data_clean_artifacts.trialinfo = sub_data_clean_artifacts.trialinfo(all_ind_trig_feed,:);
                sub_data_clean_artifacts.sampleinfo = sub_data_clean_artifacts.sampleinfo(all_ind_trig_feed,:);
                
                % Perform the independent component analysis (i.e., decompose the data)
                % Ran on every channel in the EEG
                cfg = [];
                cfg.method = 'runica'; % ICA; this is the default and uses the implementation from EEGLAB
                comp = ft_componentanalysis(cfg, sub_data_clean_artifacts); % it doesn't have the layout because the plots aren't necessary
                
                
                % Remove the bad components and backproject the
                % data.
                
                cfg = [];
                cfg.component = []; % reset to zero, as it was defined in applyica
                
                % Variables needed for cleaning blinks and
                % saccades.
                cfg.compdat = comp;% independent components
                cfg.sigdat = sub_data_clean_artifacts; % EEG data
                cfg.sig_index = [1, 54, 31]; % FP1 = 1 and FT7 = 54 that correlate the most with ocular and ECG artifacts. 31 = ECG channel
                cfg.p_thresh = 0.05; % pvalue thres .05
                cfg.r_thresh = 0.25; % r value thres .25 - increased the r value, to only select the components that are noise, but not activity of interest
                
                % Function - detect independent components that represent blinks
                % and saccades.
                [all_artifacts] = corr_artifact_blink_sacc_EKG(cfg);
                
                
                % Remove the artifacts
                cfg = [];
                cfg.channel = [1:30 32:63]; % Old Channel 31 is missing, replaced with ECG channel as new 31.
                data_clean_ICA = ft_preprocessing(cfg, data_clean_artifacts);
                cfg = [];
                cfg.component = all_artifacts; % Takes the components whose signal most correlate with the eye channels
                data_seg_clean = ft_rejectcomponent(cfg, comp, data_clean_ICA); % extract components from signal
                data_seg_clean.cfg = [];
                kept_trials = unique(data_seg_clean.trialinfo(all_ind_trig_feed,3));
                
                save ([rd.dirSav, 'all_artifacts.mat'],'all_artifacts')
                save ([rd.dirSav,'group.mat'],'group');
                save ([rd.dirSav,'kept_trials.mat'], 'kept_trials');
                save ([rd.dirSav,['data_seg_clean_' num2str(rd.subj) '.mat']],'-v7.3', 'data_seg_clean');
                % Print a text document in the save folder of
                % how many trials kept.
                %fileID = fopen([rd.dirSav,[num2str(length(kept_trials)) '.txt']],'w');
                %fprintf(fileID,'Number of Trials');
                %fclose(fileID);
                
                clear comp all_ind_trig_feed all_artifacts sub_data_clean_artifacts clear keeptrials jjj alltrials trialnum trialnumb cfg
                
            case 'conditions'
                
                
                %load ([rd.dirSav 'data_seg_clean_' num2str(rd.subj) '.mat'])
                %load ([rd.dirSav,'kept_trials.mat'], 'kept_trials');
                %load ([rd.dirSav,'group.mat'],'group');
                %load ([rd.dirData 'BEHAV\events.mat'])

                cfg = [];
                cfg.reref         = 'yes';
                %cfg.refchannel    = rd.listChan([1:30,32:end]); 
                % average of all active channels excluding ECG 
                % Exluding channels 1 and 54 as they were used in the ICA
                % corrections.
                cfg.refchannel    = rd.listChan([2:30,32:53, 55:end]);
                cfg.refmethod     = 'avg';
                data_clean_reref  = ft_preprocessing(cfg, data_seg_clean);
                data_clean_reref.cfg.previous =[];
                % Doing this here because for some reason in the
                % define_epochs I decided to make Correct = 1 and Incorrect
                % = 2. This is confusing to me now, so I am changing it to
                % Correct = 1 and Incorrect = 0.
                tempincor = find(data_clean_reref.trialinfo(:,2) == 2);
                data_clean_reref.trialinfo(tempincor,2) = 0;
                
                %Data segmented solely by trial.
                
                trig_feed = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44,51,52,53,54];
                all_ind_trig_feed =[];
                for jj = 1:length(trig_feed)
                    trigfeed = trig_feed(jj);
                    indtrigfeed = find(data_clean_reref.trialinfo(:,1) == trigfeed);
                    all_ind_trig_feed =[all_ind_trig_feed;indtrigfeed];
                end
                % Sorting as a precaution, not sure if it affects anything
                % or not.
                all_ind_trig_feed = sortrows(all_ind_trig_feed);
                
                cfg = [];
                cfg.trials = all_ind_trig_feed;
                cfg.demean = 'yes';
                cfg.baselinewindow = [-0.1 -0.05]; %Never have to change this. It's based on time and not samples.
                data_seg_by_feedback  = ft_preprocessing(cfg,data_clean_reref);
                
                % Adding model data.
                data_seg_by_feedback.model.PEs = events.model.PEs_resigned_all(kept_trials)';
                data_seg_by_feedback.model.Value = events.model.V_all(kept_trials,:);
                data_seg_by_feedback.model.ProbaCh = events.model.ProbaCh_all(kept_trials);
                
                
                save ([rd.dirSav,['data_seg_by_feedback_' num2str(rd.subj) '.mat']],'-v7.3', 'data_seg_by_feedback');
                
                
                
                %HEP condition by trial
                all_good_trials_ind = [];
                for ji = 1:length(kept_trials)
                    good_trials = kept_trials(ji);
                    ind_good_trial = find(data_clean_reref.trialinfo(:,3) == good_trials);
                    all_good_trials_ind =[all_good_trials_ind;ind_good_trial];
                end
                
                all_ind_HEP_feed = [];
                HEP_feed = [611,612,613,614,621,622,623,624,631,632,633,634,641,642,643,644,651,652,653,654];
                for jj = 1:length(HEP_feed)
                    HEPfeed = HEP_feed(jj);
                    indHEPfeed = find(data_clean_reref.trialinfo(:,1) == HEPfeed);
                    all_ind_HEP_feed =[all_ind_HEP_feed;indHEPfeed];
                end
                % Intersect is kind enough to sort for us,
                % so we don't have to.
                [hepcondidx] = intersect(all_good_trials_ind, all_ind_HEP_feed);
                clear all_ind_trig_feed
                
                cfg = [];
                cfg.trials = hepcondidx;
                cfg.demean = 'yes';
                cfg.baselinewindow = [-0.2 -0.05];
                HEP_all_by_trial  = ft_preprocessing(cfg,data_clean_reref);
                [~, trlinfoidx] = intersect(HEP_all_by_trial.trialinfo(:,3), unique(HEP_all_by_trial.trialinfo(:,3)));
                
                for jkj = 1:length(HEP_all_by_trial.trialinfo)
                    HEP_all_by_trial.model.Value(jkj,:) = events.model.V_all(HEP_all_by_trial.trialinfo(jkj,3),:);
                    HEP_all_by_trial.model.PEs(jkj,1) = events.model.PEs_resigned_all(HEP_all_by_trial.trialinfo(jkj,3));
                    HEP_all_by_trial.model.ProbaCh(jkj,1) = events.model.ProbaCh_all(HEP_all_by_trial.trialinfo(jkj,3));
                end
                alltrials = unique(HEP_all_by_trial.trialinfo(:,3));
                for j = 1:length(alltrials)
                    curtrial = alltrials(j);
                    idx = find(HEP_all_by_trial.trialinfo(:,3) == curtrial);
                    HEPnum = 1;
                    for k = 1:length(idx)
                        HEP_all_by_trial.trialinfo(idx(k),4) = HEPnum;
                        HEPnum = HEPnum + 1;
                    end
                end
                
                save ([rd.dirSav,['HEP_all_by_trial_' num2str(rd.subj) '.mat']],'-v7.3', 'HEP_all_by_trial');
                
                
                
                % HEP Average.
                
                HEP_average_by_trial_temp = [];
                trials = unique(HEP_all_by_trial.trialinfo(:,3));
                for a = 1:length(trials)
                    cfg = [];
                    all_trig_con =[];
                    trial2avg = trials(a);
                    indx = find(HEP_all_by_trial.trialinfo(:,3) == trial2avg);
                    cfg.trials = indx;
                    data_cond  = ft_preprocessing(cfg,HEP_all_by_trial);
                    data_cond.cfg.previous = [];
                    cfg = [];
                    avg_data = ft_timelockanalysis(cfg, data_cond);
                    HEP_average_by_trial_temp{a} =  avg_data;
                end
                
                
                %This just restructures the output to be a lot more like
                %the normal output from Fieldtrip, as this is what is used
                %in all the other code. It's unneccassary if you want to
                %refactor everything for specific segments, but is better
                %to keep uniform code later.
                HEP_average_by_trial = [];
                HEP_average_by_trial.label = HEP_average_by_trial_temp{1}.label;
                HEP_average_by_trial.dimord = HEP_average_by_trial_temp{1}.dimord;
                HEP_average_by_trial.cfg = HEP_average_by_trial_temp{1}.cfg;
                HEP_average_by_trial.trialinfo = HEP_all_by_trial.trialinfo(trlinfoidx,:);
                HEP_average_by_trial.sampleinfo = HEP_all_by_trial.sampleinfo(trlinfoidx,:);
                HEP_average_by_trial.model.ProbaCh = HEP_all_by_trial.model.ProbaCh(trlinfoidx);
                HEP_average_by_trial.model.PEs = HEP_all_by_trial.model.PEs(trlinfoidx);
                HEP_average_by_trial.model.Value = HEP_all_by_trial.model.Value(trlinfoidx,:);
                HEP_average_by_trial.fsample = HEP_all_by_trial.fsample;
                HEP_average_by_trial.hdr = HEP_all_by_trial.hdr;
                for b = 1:length(HEP_average_by_trial_temp)
                    HEP_average_by_trial.time{b} = HEP_average_by_trial_temp{b}.time;
                    HEP_average_by_trial.trial{b} = HEP_average_by_trial_temp{b}.avg;
                    HEP_average_by_trial.var{b} = HEP_average_by_trial_temp{b}.var;
                    HEP_average_by_trial.dof{b} = HEP_average_by_trial_temp{b}.dof;
                end
                
                clear HEP_average_by_trial_temp
                save ([rd.dirSav,['HEP_average_by_trial_' num2str(rd.subj) '.mat']],'-v7.3', 'HEP_average_by_trial');
                
                
                % Now for segmented by stimulus.
                
                
                all_ind_stim_feed = [];
                stim_feed = [61,62,63,64,65,66,67,68,71,72,73,74,75,76,77,78,81,82,83,84,85,86,87,88,91,92,93,94,95,96,97,98,101,102,103,104];
                for jj = 1:length(stim_feed)
                    stimfeed = stim_feed(jj);
                    indstimfeed = find(data_clean_reref.trialinfo(:,1) == stimfeed);
                    all_ind_stim_feed =[all_ind_stim_feed;indstimfeed];
                end
                all_ind_stim_feed = sort(all_ind_stim_feed);
                
                cfg = [];
                cfg.trials = all_ind_stim_feed;
                cfg.demean = 'yes';
                cfg.baselinewindow = [-0.2 -0.05];
                data_seg_by_stim  = ft_preprocessing(cfg,data_clean_reref);
                % Bit late but need to correct the trialinfo so they have
                % the outcome column. Gotta use event.mat to do this.
                for jlj = 1:length(data_seg_by_stim.trialinfo)
                    data_seg_by_stim.trialinfo(jlj,2) = events.behavior.Reward(data_seg_by_stim.trialinfo(jlj,3));
                end

                
                
                % Can now do anything to this, average it by any condition,
                % etc. No model data for now.
                save ([rd.dirSav,['data_seg_by_stim_' num2str(rd.subj) '.mat']],'-v7.3', 'data_seg_by_stim');
                
        end
        
        
    end
end
clearvars -except r t s todo dwnsmple allppts_HEP_corr_conds allpts_all_conds_ave ncond nchan rd surrogate
