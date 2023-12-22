function [V,pe,pch,beta,lrate,chcorr,err,bic,choice,predch,dv]=rl_fits2_prediction(X)

%% inputs:
% X matrix as N x 5 (N = trials)
% 1st column: Choices for option A
% 2nd column: Choices for option B
% 3th column: Reward/Feedback Received 
% 4th column: Start of new block

%% outputs:
% V         : Matrix of value for the two options
% pe        : Prediction error on trial-by-trial basis
% pa        : probability of choice for option A
% pb        : probability of choice for option B
% alpha     : indecision point 
% beta      : temperature parameter
% lrate     : learning rate
% err       : returning error
% bic       : BIC as BIC = -2*logL + M*log(N)/N - M: fitted parameters, N: trials

%% References:
% 1) Fouragnan et al. 2019
% 2) Fouragnan et al. 2017
% 3) Fouragnan et al. 2015

%% get fitted parameters:
[beta,lrate,chcorr,err] = quickfit_lr2_prediction(X);

%% Some initializations
n = size(X,1);
V = zeros(n,2); 
V(1,1)=0; V(1,2)=0; 
pe = zeros(n,1);
pch = zeros(n,1);

ca = X(:,1);	 % binary vector with item 1 choices
choice = ca;
choice(choice==0)=-1;
rt = X(:,3);  % binary vector with rewarded choices
reward = X(:,3);
% change this to actual color:
for i=1:length(rt)
   if choice(i)==1 && rt(i)==1
       rt(i)=1;
   elseif choice(i)==-1 && rt(i)==1
       rt(i)=-1;
   elseif choice(i)==-1 && rt(i)==0
       rt(i)=1;
   elseif choice(i)==1 && rt(i)==0
       rt(i)=-1;
   end
end
bk = X(:,4);  % binary vector with start of new block
stim1 = X(:,5);  % binary vector with start of new block
stim2 = X(:,6);  % binary vector with start of new block
cue=(sum([stim1 stim2],2)/10)-1;
previousch(1) = 0;

n = length(ca);

for i=1:n

    if i>1
        previousch(i) = choice(i-1);
    end
    
    if cue(i)==1
        predch(i) = 0.5*V(i,1) + 0.5*V(i,1);
        pe(i) = rt(i) - predch(i);
        V(i+1,1) = V(i,1) + lrate*pe(i);
        V(i+1,2) = V(i,2);
    elseif cue(i)==2
        predch(i) = 0.5*V(i,1) + 0.5*V(i,2);
        pe(i) = rt(i) - predch(i);
        V(i+1,1) = V(i,1) + lrate*pe(i);
        V(i+1,2) = V(i,2) + lrate*pe(i);
    elseif cue(i)==3
        predch(i) = 0.5*V(i,2) + 0.5*V(i,2);
        pe(i) = rt(i) - predch(i);
        V(i+1,2) = V(i,2) + lrate*pe(i);
        V(i+1,1) = V(i,1);
    end
    
    dv(i) = beta*(predch(i) - .5) + chcorr*previousch(i);
	pch(i) = 1 / (1 + exp(-dv(i)));
	
end

% BIC 
V = V(1:n,:);
m = 3; 

[aic,bic]= aicbic(err,m,n);