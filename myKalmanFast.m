function[xhat, y,K,P] = myKalmanFast(x,Q,R)
%  function[xhat, y,K,P] = myKalmanFast(x,Q,R)
%
%   Takes as input target position (x), target displacement variance (Q),
%   and observation noise variance (R).
%
%   Returns position estimates (xhat) for an observer with observation 
%   noise variance R, the associated noisy sensory observations (y), the
%   Kalman gain (K) and the posterior variance (P). See figure 5 and
%   appendix B.
% 

N = length(x);                          % Number of samples

v = randn(N,1)*sqrt(R);
y = x + v;                              % Eq. B2

P = Q/2 * (sqrt(1+4*R/Q)-1);            % Eq. B8
K = (P+Q)/(P+Q+R);                      % Eq. B10
D = spdiags([ones(N,1),(K-1)*ones(N,1)],[0,-1],N,N);    % Eq. B12
xhat = D\(K*y);                         % Eq. B13



