function [mu,sigma,flagsigma,covmat ] = compGauss_ms(m,h,vargin)
% compute group-level gaussian from fminunc computed parameters and their covariances
% MKW, 2017
%
% INPUT:    - m:  fitted parameters (npar x nsub matrix)
%           - h:  individual-level hessians(npar x npar x nsub)
%           - vargin: if set to 2, computes covariance matrix in addition
%
% OUTPUT:   - group mu and group sigma 
%           - flagcov: flag indicating whether model variance was calculated successfully 
%           - covmat: full covariance matrix; is [] if no vargin specified
%%

% get info
nsub = size(m,2);
npar = size(h,1);
covmat = [];

% ------ 1) compute mean: -------------------------------------------------
mu =  mean(m,2);


% ------2) Compute sigma: -------------------------------------------------

sigma   = zeros(size(h,1),1);

for is = 1:nsub
   sigma = sigma + m(:,is).^2 + diag(pinv(h(:,:,is)));
end
sigma = sigma./nsub  - mu.^2;

% give error message in case:
if min(sigma) < 0, flagsigma = 0; disp('CovError!'); else flagsigma=1;  end % negative values for the variance cannot be



% ----- 3) Optional: Get full covariance matrix----------------------------


if nargin < 3, return; end % go on only if vargin is defined
covmat   = zeros(npar,npar);
if vargin == 2 
   for is=1:nsub
      covmat = covmat + m(:,is)*m(:,is)' - m(:,is)*mu' - mu*m(:,is)' + mu*mu' + pinv(h(:,:,is));
   end
   covmat = covmat./nsub;
end
if det(covmat)<=0
   fprintf('negative/zero determinant - prior covariance not updated');
end

end

