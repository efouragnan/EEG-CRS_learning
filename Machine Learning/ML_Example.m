% Template for creating ML documents. Using the HEP average by trial data.
% X, truth, duration, skipLOO, perm, eigratio,gamma, method)
% Sample length = 1000 Hz Downsampled to 500 Hz


tau     = 5;                        % 10 ms of increment for discrimination
duration = 30;                      % 60 ms of windows length for discrimination
windInt = [1 401-30];               % feedback loocked onsets times from -200 to 600ms
stepWin = 1:tau:sum(abs(windInt));  % number of steps in slicing window approach
timex = ((windInt) - 101);
timexStep = 1:tau:timex(end);
eigratio = 0; gamma = 0;            % Not doing this right now
method = 1;                         % Loop through fisher discriminant--
skipLOO = 0;                        % Don't skip LOO
opt = 1;
perm = [0 100 0.05];


rd.dirProj =  ('\project_directory');
rd.listSubj = {'02','04','05','06','07','08','10','11','12','13','14','15','16','17','18','19','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36'};


rd.nSubj = length(rd.listSubj);

        for i = 1:rd.nSubj
            rd.subj = rd.listSubj{i};
            rd.dirData = fullfile(rd.dirProj,['sub' rd.subj],'/Downsampled/');
            rd.dirSav = fullfile(rd.dirProj,['sub' rd.subj],'/ML_HEP/');
            mkdir(rd.dirSav);
            load ([rd.dirData,['HEP_average_by_trial_' num2str(rd.subj) '.mat']]);
            condPEvNeg = find(abs(HEP_average_by_trial.model.PEs) < prctile(abs(HEP_average_by_trial.model.PEs), 25)); % 25th percentile and below PEs
            condPEvPos = find(abs(HEP_average_by_trial.model.PEs) > prctile(abs(HEP_average_by_trial.model.PEs), 75)); % 75th percentile and above PEs
            EEG.full = cat(3, HEP_average_by_trial.trial{:}); % EEG is stored in a cell array, converting it to a 3D matrix
            EEG.full = baselinecorr(EEG.full,75,1); %Baseline correction is done per HEP, but not when averaging. Do it again here, but only for averaged data.
            EEG.vPos = EEG.full(:,:,condPEvPos);
            truth = [zeros(1,length(condPEvPos)) ones(1,length(condPEvNeg))]; %B = 1's, A = 0's in single_trial_analysis.EEG
            idx = 1;
            for j = stepWin
                X_vNeg = EEG.vNeg(:, j:duration+j, :);
                X_vPos = EEG.vPos(:, j:duration+j, :);
                
                X = [reshape(X_vNeg, [size(X_vNeg,1), size(X_vNeg,2) * size(X_vNeg,3)]) reshape(X_vPos, [size(X_vPos,1), size(X_vPos,2) * size(X_vPos,3)])];
                if opt %0 0.1 1
                    
                    idx2 = 1;
                    gamma_v = 0:0.1:1;
                    for f = gamma_v
                        [TempAz(idx2)] = single_trial_analysis_EEG(X,truth,duration,skipLOO,perm,eigratio,f,method);
                        idx2 = idx2 + 1;
                    end
                    [~, idxmax] = max(TempAz);
                    igamma = gamma_v(idxmax);
                    tuned_igamma = [igamma-0.1:0.01:igamma+0.1];
                    tuned_igamma = tuned_igamma(tuned_igamma >= 0 & tuned_igamma <= 1);
                    clear TempAz
                    idx2 = 1;
                    for g = tuned_igamma
                        [TempAz(idx2)] = single_trial_analysis_EEG(X,truth,duration,skipLOO,perm,eigratio,g,method);
                        idx2 = idx2 + 1;
                    end
                    [~, idxmax] = max(TempAz);
                    gamma = tuned_igamma(idxmax);
                    clear TempAz
                end
                [Azloo(idx),~,Y(idx),a(:,idx),v(:,idx)] = single_trial_analysis_EEG(X,truth,duration,skipLOO,perm,eigratio,gamma,method);
                gamma_all(idx) = gamma;
                idx = idx + 1;
            end
            analysis.Azloo = Azloo;
            analysis.Y = Y;
            analysis.Ya = Y.condA;
            analysis.Yb = Y.condB;
            analysis.a = a;
            analysis.v = v;
            analysis.gamma = gamma_all;
            
            save([rd.dirSav,['analysis_HEP_PE_Pos_vs_Neg' num2str(rd.subj) '.mat']], 'analysis', '-v7.3');

          
            clear analysis EEG condCorr condIncorr truth
        end
        
