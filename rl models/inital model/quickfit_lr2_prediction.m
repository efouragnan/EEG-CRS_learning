function [beta,lrate,chcorr,err]=quickfit_lr2_prediction(data)

global Data
Data = data;

% initialisation
q(1,1) = 3;    % beta
q(2,1) = 0.4;  % learning rate
q(3,1) = 0.1;

% search for hidden variables
quick = fmincon(@(q) fitfunrl2_prediction(q), q, [],[],[],[],[0; 0; -1],[10; 1; 1],[],optimset('Display','off','TolX',.0001,'Algorithm','interior-point'));

% quick
beta  = (quick(1,1));
lrate = (quick(2,1));
chcorr = (quick(3,1));
err = fitfunrl2_prediction(quick);
