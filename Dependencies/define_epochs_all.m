% Alex Sel 18/10/18
% Outputs all the trials time-locked to the visual onset, correct and
% incorrect, as coded by last colunm (1 = correct, 2 = incorrect)

function [all_trl_visual, all_trl_decision, all_trls_feed, trl_visual_feed, trl_R_feed, trl_R_feed2, group] = define_epochs_all(cfg,R_trl,rd,data_filtered)


%define columns in trl_visual_feed
trlbegtm =1;
trialendtm =2;
trigcode = 4; % triger type
resp = 5; % ppt resp correct vs inco
trialltype = 6; % trial type corr vs incor
ordtrlblck =7; % order trial in the block
ordrtrlexp =8; % order trial in the experiment

%firstITI =110; % First ITI of the next trial, right after the feedback
%lasttirgblock=99;% trigger indicating end of the block
%lasttrrev =29; % last trl before reversal

% Distance between beggining of two consecutive trials in sampling points(?)
% dist_trl = 14000;

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);


% change trigger string to numbers
for i = 1:length(event)
    event(i).value = event(i).value(2:end);
    event(i).value = str2double(event(i).value);
end
% Define participant's group
if find([event.value] == 31)
    group = 1; 
else
    group = 2; 
end

if group == 1 
    cfg.trialdef.eventvalue = rd.listCond_feedback_g1; 
else
    cfg.trialdef.eventvalue = rd.listCond_feedback_g2; 
end

                if group == 1
                    cfg.trialdef.viseventvalue = rd.listCond_visual_g1;
                else
                    cfg.trialdef.viseventvalue = rd.listCond_visual_g2;
                end

% search for "trigger" events
value  = [event(find(strcmp('Stimulus', {event.type}))).value]';
sample = [event(find(strcmp('Stimulus', {event.type}))).sample]';

factor = hdr.Fs / data_filtered.fsample;

sample = round((sample-1)/factor +1);
% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * data_filtered.fsample);
posttrig =  round(cfg.trialdef.post * data_filtered.fsample);

vispretrig = -round(cfg.trialdef.vispre * data_filtered.fsample); 
visposttrig =  round(cfg.trialdef.vispost * data_filtered.fsample);

decpretrig = -round(cfg.trialdef.decpre * data_filtered.fsample); 
decposttrig =  round(cfg.trialdef.decpost * data_filtered.fsample);
% To sum up, this finds the first INI of each trial, the 110 code, and then
% shoves everything from there to the next 110 into a structure.
% After, it loops through the structure for each 'trial' to pull out the
% important information, the response, keyboard trig code, calculate
% trlbegin and trlend, etc.
% Unfortunately, for Pilot_02, keyboard trig code & response are coded the
% same (1 to 3), therefore the keyboard trig code is left as a column of 0s.

allvalues = value;
trial = [];
fulltrl = struct('idx',0,'data',0,'trlbegin',0,'trlend',0,'offset',0,'trialtrigger',0,'response',0,'trigcodes',0,'trialnumber',0);
trl = str2double(cfg.trialdef.eventvalue);
vistrl = str2double(cfg.trialdef.viseventvalue);
indextrial = 0;
responses = [1,2,3];
trigcodes = [10,20];
vistrials = [];
dectrials = [];
for j = 1:length(allvalues)
    trgcode = allvalues(j);
    idx = j;
    if trgcode == 110
        if indextrial == 0
            indextrial = indextrial + 1;
            fulltrl(indextrial).idx = idx;
            trial = trgcode;
            
        elseif indextrial > 0
            indextrial = indextrial + 1;
            fulltrl(indextrial).idx = idx;
            trial = trgcode;
            
        end
    elseif trgcode ~= 110
        trial = [trial, trgcode];
        for z = 1:length(trl)
            for i = 1:length(trial)
                if trial(i) == trl(z)
                    fulltrl(indextrial).trlbegin = (sample(j) + pretrig);
                    fulltrl(indextrial).trlend = (sample(j) + posttrig);
                    fulltrl(indextrial).offset = pretrig;
                end
            end
        end
        for y = 1:length(vistrl)
                if trgcode == vistrl(y)
                    vistrlbegin = (sample(j) + vispretrig);
                    vistrlend = (sample(j) + visposttrig);
                    visoffset = vispretrig;
                    vistrlnum = indextrial;
                    newvis = [vistrlbegin vistrlend visoffset trgcode vistrlnum];
                    vistrials = [vistrials; newvis];
                end
        end
            if trgcode == 100
                decbegin = (sample(j) + decpretrig);
                decend = (sample(j) + decposttrig);
                decoffset = decpretrig;
                dectrl = indextrial;
                newdec = [decbegin decend decoffset trgcode dectrl];
                dectrials = [dectrials; newdec];
            end
    end
    if indextrial == 0
        trial = [];
    else
        fulltrl(indextrial).data = [trial];
        
    end
end

% Loading behavioural data here because we need a better way to identify
% bad trials. Comparing to the behavioural seems to be the best way.
load(fullfile(rd.dirData,'\BEHAV\events.mat'))
badtrials = [];
for k = 1:length([events.behavior.Missed])
    if any(ismember(trl, [fulltrl(k).data])) == 0
        events.behavior.Missed(k) = 1;
        badtrials = [badtrials, k];
    end
    for p = 1:length(fulltrl(k).data)
        if ismember(fulltrl(k).data(p), trl) == 1
            fulltrl(k).trialtrigger = fulltrl(k).data(p);
        end
        
        if ismember(fulltrl(k).data(p), responses) == 1
            fulltrl(k).response = fulltrl(k).data(p);
        end
        if ismember(fulltrl(k).data(p), trigcodes) == 1
            fulltrl(k).trigcodes = fulltrl(k).data(p);
        end
    end
    if isempty(fulltrl(k).trigcodes) == 1
        fulltrl(k).trigcodes = 0;
    end
end

% Often the response is grabbed incorrectly from the EEG. The behavioural
% data is much more reliable so comparing it and replacing
% irregularities.
for k = 1:length([events.behavior.Reward])
    if events.behavior.Reward(k) == 0
        events.behavior.Reward(k) = 2;     %Recoding events.behavior.Reward as it is a 0/1 binary reverse not 1/2
    end
    if isempty(fulltrl(k).response) == 1
        fulltrl(k).response = 0;
    end
    if fulltrl(k).response ~= events.behavior.Reward(k)
        fulltrl(k).response = events.behavior.Reward(k);
    end
end

missingtrials = [events.behavior.Missed] == 1;




fulltrl(missingtrials) = [];
fulltrlNewunits = events.behavior.Newunit;

fulltrlNewunits(missingtrials) = [];

newunit = find(fulltrlNewunits == 1);

%firstppts = {'02','04','05','06','07','08','09','10','11','12'};
%if ismember(rd.subj,firstppts)
    
    if group ==1
    cond = str2double({'11','13','21','23','31','33','51','53'});
    cond2 = str2double({' 61',' 62',' 63','64',' 71',' 72',' 73',' 74',' 81',' 82',' 83',' 84'}); %
    else
        cond = str2double({'11','13','21','23','41','43','51','53'});
        cond2 = str2double({' 61',' 62',' 63','64',' 71',' 72',' 73',' 74',' 91',' 92',' 93',' 94'}); %
    end
    for a = 1:length(cond)
        trig = cond(a);
        trlindex = [];
        for b = 1:length(fulltrl)
            if fulltrl(b).trialtrigger == trig
                trlindex = [trlindex, b];
            end
        end
    [~, matchidx] = intersect(trlindex,newunit);
    for c = trlindex(matchidx(2)):trlindex(end)
        fulltrl(c).trialtrigger = fulltrl(c).trialtrigger + 1;
    end
    end

    
    % Someone failed the first trial of a block, which ruins this script.
    % This 'fixes' it.
if str2num(rd.subj) == 31
events.behavior.Newunit(195) = 0;
elseif str2num(rd.subj) == 17
    events.behavior.Newunit(260) = 0;
end
    
newunit = find(events.behavior.Newunit == 1);



for a = 1:(length(cond2)/2)
    b = 1:2:length(cond2);
    trig = [cond2(b(a)) cond2(b(a)+1)];
  vistrialnumbs = vistrials(ismember(vistrials(:,4),trig),5);
  visindex = find(ismember(vistrials(:,4),trig) == 1);
  [~, matchidx] = intersect(vistrialnumbs,newunit);
  vistrials(visindex(matchidx(2)):visindex(end),4) = vistrials(visindex(matchidx(2)):visindex(end),4) + 4;
end


cond2 = str2double({'101','103'});
for a = 1:length(cond2)
    trig = cond2(a);
    vistrialnumbs = vistrials(ismember(vistrials(:,4),trig),5);
    visindex = find(ismember(vistrials(:,4),trig) == 1);
    [~, matchidx] = intersect(vistrialnumbs,newunit);
    vistrials(visindex(matchidx(2)):visindex(end),4) = vistrials(visindex(matchidx(2)):visindex(end),4) + 1;
end
% for c = 1:length(newunit)
%trialblock = find(vistrials(:,5) == newunit(c));
%   if events.behavior.Newblock(newunit(c))


%   elseif events.behavior.Newblock(newunit(c)) == 1 &&
%   vistrials(trialblock(1),5) == any([101,103]) && newunit(c+2) == any([101,103])
           
%   elseif newunit(c) == newunit(end) && vistrials(newunit(c)) ~= any([101,103])
%         
%         for d = newunit(c):newunit(end)
%         vistrials(d,4) = vistrials(d,4) + 4;
%         end
%   elseif vistrials(newunit(c)) == any([101,103])
%         for d = newunit(c):newunit(c+2)-1    
%             vistrials(d,4) = vistrials(d,4) + 1;
%         end
%         events.behavior.Newblock(newunit(c+4)) = 1;
%   else
%         
%         for d = newunit(c):newunit(c+1)-1    
%         vistrials(d,4) = vistrials(d,4) + 4;
%         end
%         
%    end
% end
    

% for c = 1:length(newunit)
%     newunit(c);
%     if events.behavior.Newblock(newunit(c)) == 1
%     elseif newunit(c) == newunit(end)
%         for d = newunit(c):newunit(end)
%             
%     fulltrl(d).trialtrigger = fulltrl(d).trialtrigger + 1; 
%     
%         end
%     else
%         for d = newunit(c):newunit(c+1)-1
%             
%     fulltrl(d).trialtrigger = fulltrl(d).trialtrigger + 1; 
%     
%         end
%     end
% end
% badppt = {'02'}; %Manually updating anyone who has random erroneous 'correct' trial defs, which are made incorrect by the code
% if ismember(rd.subj,badppt)
% fulltrl(183).trialtrigger = 24;
% end 

%End of 'if participants'
%end

% missingtrials = [events.behavior.Missed] == 1;
% 
% 
% 
% 
% fulltrl(missingtrials) = [];



% deltrials = find(missingtrials == 1);
% dectrials(deltrials,:) = [];
% 
% for e = 1:length(deltrials)
%     d = deltrials(e);
% vistrials(find(vistrials(:,5) == d),:) = [];
% end
% %vistrials([find(missingtrials ==1)],:) = [];
 for k = 1:length([fulltrl.idx])
     fulltrl(k).trialnumber = k;
     %dectrials(k,5)= k;
 end



 
% 
%     for q = 1:2:length(vistrials)
%         fulltrl((q+1)/2).trialnumber;
%         vistrials(q,6) = fulltrl((q+1)/2).trialnumber;
%     end
%     for q = 2:2:length(vistrials)
%         fulltrl((q)/2).trialnumber;
%         vistrials(q,6) = fulltrl((q)/2).trialnumber;
%     end



cond = unique([fulltrl(1:end).trialtrigger]);
for c = 1:length(cond)
    trig = cond(c);
    numbtrls = length(find([fulltrl(:).trialtrigger] == trig)); disp(['There are ' num2str(numbtrls)  ' trials in condition '  num2str(trig)]);
end


cond = unique(vistrials(:,4)); trlnumb = [];
for cc = 1:length(cond)
    condi = cond(cc);
    numbtrls = length(find(vistrials(:,4)==condi)); disp(['There are ' num2str(numbtrls)  ' trials in condition '  num2str(condi)]);
    trlnumb = [trlnumb,numbtrls];
end



% Creating trl_visual_feed in the smame format as before and is used
% throughout the scripts.

clear trl
trl_visual_feed = [];
newtrl = [];
cond = unique([fulltrl(1:end).trialtrigger]);
for k = 1:length(cond)
    trl = cond(k);
    trltype = [];
    for p = 1:length([fulltrl.trialtrigger])
        if fulltrl(p).trialtrigger == trl
            newtrl = [fulltrl(p).trlbegin, fulltrl(p).trlend, fulltrl(p).offset, fulltrl(p).trialtrigger, fulltrl(p).response, fulltrl(p).trigcodes];
            trltype = [trltype; newtrl];
        end
    end
    trlnumb = [1:1:size(trltype,1)];
    trltype(:,end+1) = trlnumb;
    trl_visual_feed = [trl_visual_feed; trltype];
    
    clear trltype
end

trl_visual_feed = sortrows(trl_visual_feed);
index = (1:length(trl_visual_feed)).';
trl_visual_feed(:,ordrtrlexp)= index;
trl_visual_feed = trl_visual_feed(:, [1:6 flip([7,8])]);
%Now we have trl_visual_feed we can create the rest using it. Code from old
%define_epocks_feedback

% SUBJ 31 HAD RECORDING CUT EARLY, SO NO HEARTBEATS FOR
% LAST TRIAL
if str2num(rd.subj) == 31
    trl_visual_feed(480,:) = [];
end

cond = unique(trl_visual_feed(:,trigcode)); trl_R_feed =[]; trl_R_feed2 = [];
for d = 1:length(cond)
    trig = cond(d);
    ind = find([trl_visual_feed(:,trigcode)] == trig);
    %trlnumbR = trl_visual_feed(ind,end);
    trbg =[]; trend =[]; ind_all_R =[]; trl_NR_all =[]; all_resptrial = []; all_trialtype =[]; ind_all_R2 = []; trl_NR_all2 =[]; all_resptrial2 = []; all_trialtype2 =[];
    for k = 1:length(ind)
        t = ind(k);
        trbg = trl_visual_feed(t,trlbegtm)+(1100/factor); trend = trl_visual_feed(t,trialendtm);% added 1100ms as we only want to take R-waves from +200 after visual onset(because of the baseline of -200 from the R-wave), taking into account the -900ms baseline
        indR = find([R_trl(:,1)] > trbg & [R_trl(:,1)] < trend);
        indR2 = [indR(1) - 1; indR];
        ind_all_R =[ind_all_R;indR];
        ind_all_R2 = [ind_all_R2;indR2];
        trl_N_R = zeros(length(indR),1) + ind(k); % add trial number
        trl_N_R2 = zeros(length(indR2),1) + ind(k);
        %         resptrial = trl_visual_feed(trl_N_R,resp); trialtype = trl_visual_feed(trl_N_R,trialltype);
        resptrial = trl_visual_feed(trl_N_R,resp);
        resptrial2 = trl_visual_feed(trl_N_R2,resp);
        trialtype = trl_visual_feed(trl_N_R,trialltype);
        trialtype2 = trl_visual_feed(trl_N_R2,trialltype);
        trl_NR_all =[trl_NR_all;trl_N_R]; 
        trl_NR_all2 =[trl_NR_all2;trl_N_R2]; 
        all_resptrial = [all_resptrial;resptrial]; 
        all_resptrial2 = [all_resptrial2;resptrial2]; 
        all_trialtype =[all_trialtype;trialtype];
        all_trialtype2 =[all_trialtype2;trialtype2];
    end
    clear indR trl_N_R
    r_trl_cond = R_trl(ind_all_R,:); r_trl_cond(:,trigcode)= trig + 600; r_trl_cond(:,resp) = all_resptrial; r_trl_cond(:,trialltype) = all_trialtype; r_trl_cond(:,ordtrlblck) = trl_NR_all;
    trl_R_feed =[trl_R_feed; r_trl_cond];
    
    r_trl_cond2 = R_trl(ind_all_R2,:); r_trl_cond2(:,trigcode)= trig + 600; r_trl_cond2(:,resp) = all_resptrial2; r_trl_cond2(:,trialltype) = all_trialtype2; r_trl_cond2(:,ordtrlblck) = trl_NR_all2;
    trl_R_feed2 =[trl_R_feed2; r_trl_cond2];
end
trl_R_feed= sortrows(trl_R_feed);
index = (1:length(trl_R_feed)).';
trl_R_feed(:,ordrtrlexp)= index ;

trl_R_feed2= sortrows(trl_R_feed2);
index2 = (1:length(trl_R_feed2)).';
trl_R_feed2(:,ordrtrlexp)= index2 ;


all_trls = [trl_visual_feed; trl_R_feed];
all_trls_feed = sortrows(all_trls,1);

all_trl_visual = vistrials;
all_trl_decision = dectrials;


end