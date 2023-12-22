function err = fitfunrl2_prediction(q)
%FITFUNRL Used by QUICKFITRL.

global Data
TINY = 1e-100;

ca = Data(:,1);	 % binary vector with item 1 choices
cb = Data(:,2);	 % binary vector with item 2 choices
choice = ca;
choice(choice==0)=-1;
rt = Data(:,3);  % binary vector with rewarded choices
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
bk = Data(:,4);  % binary vector with start of new block
stim1 = Data(:,5);  % binary vector with start of new block
stim2 = Data(:,6);  % binary vector with start of new block
cue=(sum([stim1 stim2],2)/10)-1;

n = length(ca);
na = sum(ca);
nb = sum(cb);

% Some initializations
V = zeros(n,2);
V(1,1)=0; V(1,2)=0;
previousch(1) = 0;
pe = zeros(n,1);
pch = zeros(n,1);

beta  = q(1);
lrate = q(2);
chcorr = q(3);


for i=1:n

    if i>1
        previousch(i) = choice(i-1);
    end

    %     if bk(i)
    %         pe(i)=0;
    %         V(i,1)=0; V(i,2)=0;
    %         pch(i,1)=0;
    %         previousch(1) = 0;
    %     end

    if cue(i)==1
        predch(i) = 0.5*V(i,1) + 0.5*V(i,1);
        pe(i) = rt(i) - predch(i);
        V(i+1,1) = V(i,1) + lrate*pe(i) ;
        V(i+1,2) = V(i,2) ;

    elseif cue(i)==2
        predch(i) = 0.5*V(i,1) + 0.5*V(i,2);
        pe(i) = rt(i) - predch(i);
        V(i+1,1) = V(i,1) + lrate*pe(i) ;
        V(i+1,2) = V(i,2) + lrate*pe(i) ;

    elseif cue(i)==3
        predch(i) = 0.5*V(i,2) + 0.5*V(i,2);
        pe(i) = rt(i) - predch(i);
        V(i+1,2) = V(i,2) + lrate*pe(i) ;
        V(i+1,1) = V(i,1) ;

    end

    dv(i) = beta*(predch(i) - .5) + chcorr*previousch(i);
    pch(i) = 1 / (1 + exp(-dv(i)));

end

pch = pch - (pch > .999999)*TINY + (pch < .0000001)*TINY;

ca(ca==-1)=0;
err = -(sum(ca.*log(pch)))/na ;
