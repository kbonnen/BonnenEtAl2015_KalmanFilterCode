This package contains matlab files which can perform the Kalman filter analysis and the
ccg analysis described in Bonnen et al 2015 (JOV).

tracking2graphs.m performs the Kalman filter analysis and the CCG analysis on the data in
data.mat and produces example.fig.  The data used is for subject LKC in the paper.

negloglikelihoodr.m is a function used by tracking2graphs.m to fit positional uncertainty 
as described in appendix B.

myKalmanFast.m is a function which can be used to produce Kalman filter position estimates
for a set of target positions.