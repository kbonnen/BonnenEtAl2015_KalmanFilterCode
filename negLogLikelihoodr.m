function[nLL] = negLogLikelihoodr(rr,Q,X,Xhat)
% function[nLL] = negLogLikelihoodr(rr,Q,X,Xhat)
% 
%
%   Takes Stimulus variance (Q), log observation noise variance(rr), X (true
%   target position), Xhat (position estimates) and reports the negative
%   log likelihood calculated via the Kalman Filter model.  
%

numTrials = size(X,2);          % number of trials
N = size(X,1);                  % number of samples
nLL = 0;            
rr = exp(rr);                   % exponentiate observation noise variance                  
pp = Q/2 * (sqrt(1+4*rr/Q)-1);  % posterior variance - Eq B8

k = (pp+Q)/(pp+Q+rr);           % Kalman Gain - eq B10
d = spdiags([ones(N,1),(k-1)*ones(N,1)],[0,-1],N,N);    % Eq B12

% NOTE: The appendix describes maximum likelihood fitting of R using 
% p(xhat|x) (see Equation B14).  It is simpler (and significantly faster
% for the MATLAB implementation) to maximize the likelihood that 
% (D*xhat - K*x) ~ N(0,K^2*R)

for i=1:numTrials
    x = X(:,i);
    xhat = Xhat(:,i);
    temp = d*xhat - k*x; % eq B11
    nLL = nLL - (-1/(2*k^2*rr)*temp'*temp - N/2*log(rr) - N*log(k)); % the negative log likelihood that (D*xhat - K*x) ~ N(0,K^2*R) 
end