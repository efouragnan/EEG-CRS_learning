function [Azloo,Azsig,Y,a,v,D,Azloomedian] = single_trial_analysis_EEG(X, truth, duration, skipLOO, perm, eigratio,gamma, method)

% INPUT: 
%       X (data matrix) =
%       [Ch1: Sample1toM (Trial1) Sample1toM (Trial2) .... Sample1toM (TrialN)]
%       [Ch2: Sample1toM (Trial1) Sample1toM (Trial2) .... Sample1toM (TrialN)]
%       [....
%       [ChD: Sample1toM (Trial1) Sample1toM (Trial2) .... Sample1toM (TrialN)]
%
%       truth (trial labels - this is a vector of zeros and ones - dimensions [1 x Trials])
%       duration (length of window of interest in samples M - e.g., 60ms)
%       skipLOO (flag (0|1) for running the leave-one-cross validation
%       perm (three parameter vector - param 1: set flag (0|1) for running permutation test for Az significance 
%                                      param 2: number of permutation runs, e.g. 100)
%                                      param 3: significance threshold, p-value (e.g., 0.01) 
%       eigratio (scalar - for pca subspace reduction)
%       gamma   (scalar - regularization parameter for LDA)
%       method  (LR or LDA - 0|1)
%
% OUTPUT: 
%       Azloo (discriminator performance - this is a scalar between 0.5 and 1) 
%       Azsig (Az significance threshold - this is a scalar between 0.5 and 1) - Azsig = -1 if permutation test is turned off
%       Y (discriminator output structure - Y.condA and Y.condB single trial estimates for conditions A and B respectively)
%       a (forward mode - scalp map of discriminating activity - dimensions [1 X Channels])
%       v (spatial filter used for discrimination)
%       D (# dimensions used for discrimination)


if nargin < 3
    fprintf('Data matrix X, true labels and traing window length are required. Exiting...\n')
    return;
elseif nargin < 4
    skipLOO = 1;
    perm = [0 100 0.01];
    eigratio = 0.01;
    gamma = 0;
    method = 0;
elseif nargin < 5
    perm = [0 100 0.01];
    eigratio = 0.01;
    gamma = 0;
    method = 0;
elseif nargin < 6
    eigratio = 0.01;
    gamma = 0;
    method = 0;
elseif nargin < 7
    gamma = 0;
    method = 0;
elseif nargin < 8
    method = 0;
end

% augment truth label vector to match samples going into data matrix
truth_aug = [];
for i=1:length(truth)
    truth_aug = [truth_aug; truth(i)*ones(duration,1)];
end

% get trial numbers for condition A, B, and for all trials
Na = length(find(truth==0));
Nb = length(find(truth==1));
N = Na + Nb;

% re-sort data matrix based on trial labels
[truth_aug truthi] = sort(truth_aug);
X = X(:,truthi);

if ~method
    % run logistic regression to get weight vector
    [v,D] = logistpca(X',truth_aug,[],[],[],eigratio);
else
    % run lda to get weight vector
    [v,c] = ldareg(X',truth_aug,gamma);
    v(end+1)=-c; clear c;
    D = size(X,1);
end

% apply weight vector to get discriminator output (y)
y = X'*v(1:end-1) + v(end);
% compute forward model (scalp map for discriminating activity)
a = y \ X';

p=bernoull(1,y);
[Az] = rocarea(p,truth_aug);

y = reshape(y,[duration,N]);
Y.condA = y(:,1:Na);
Y.condB = y(:,Na+1:end);
clear yy

% number of permutation runs
nperms = perm(2);

% determine if we are bootstraping
if perm(1)
    azsigimax=nperms+1;
else
    azsigimax=1;
end

% perform leave-one-out cross validation
clear y;
ploo=zeros(N,1);

if ~skipLOO
    %fprintf('Az Significance loop: ');
    for azsigi=1:azsigimax,
        
        %fprintf('%d.',azsigi);
        if azsigi>1,
            truth_aug=truth_aug(randperm(length(truth_aug)));
        end
        
        for i=1:N
                        
            % dummy indexing for all trials
            indx = ones(N*duration,1);
            % remove one trial from the mix
            indx((i-1)*duration+1:i*duration)=0;
            
            % get remaining data
            tmp = X(:,find(indx))';
            % get truth labels for remaining trials
            tmpt = [truth_aug(find(indx))];
            
            if ~method
                % estimate weight vector with one trial removed
                [vloo D(i)] = logistpca(tmp,tmpt,[],[],[],eigratio);
            else
                % run lda to get weight vector
                [vloo,cloo] = ldareg(tmp,tmpt,gamma); 
                vloo(end+1)=-cloo;
            end
            
            % apply vector to get discriminator output for the trial we left out
            y = [X(:,(i-1)*duration+1:i*duration)' ones(duration,1)]*vloo;
            % build probability distribution for mean y (to be used in ROC)
            ploo(i)=bernoull(1,mean(y));
            ploomedian(i)=bernoull(1,median(y));
            
            clear tmp tmpt y vloo cloo indx
          
        end
        
        % run ROC analysis to get leave-one-one discriminator performance
        truthmean = ([zeros(Na,1); ones(Nb,1)]);
        if azsigi > 1
            [Azsig(azsigi-1)] = rocarea(ploo,truthmean);
        else
            [Azloo] = rocarea(ploo,truthmean);
            [Azloomedian] = rocarea(ploomedian,truthmean);
        end
        
    end
    
    if perm(1),
        psig = perm(3);
        tmpaz = sort(Azsig);
        Azsig = tmpaz(end-round(nperms*psig));
    else
        Azsig = -1;
    end
    
    if Azloo < 0.5, Azloo = 0.5; end
    if Azloomedian < 0.5, Azloomedian = 0.5; end
    %fprintf('\n');
else  
    Azloomedian = Az; 
    Azloo = Az; 
    Azsig = -1;
end
