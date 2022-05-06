
% Potential Solution

% Isnan search on data_seg
% loop through each and every cell
% if find (1) = true
% store idx number in a variable




% Specifically for 31 & 17 as both stopped recording mid trial. 

function [data_seg, trl_visual_feed, ECG_seg] = NaNremove(data_seg, trl_visual_feed, ECG_seg)
    EEG_data = cellfun(@isnan,data_seg.trial,'UniformOutput',false);
    idx = [];
    for i = 1:length(data_seg.trial)
        if isempty(find([EEG_data{i}] == 1)) == 0
            idx = [idx; i];
        end
    end

    %idx; % Uncomment if you need to see where the NaN's are.
    data_seg.trial(:,idx:end) = [];
    data_seg.time(:,idx:end) = [];
    data_seg.trialinfo(idx:end,:) = [];
    data_seg.sampleinfo(idx:end,:) = [];
    data_seg.cfg.trl(idx:end,:) = [];
    
    trl_visual_feed(end,:) = [];

clear junkdata idx i 
end



