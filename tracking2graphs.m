%% EXAMPLE - The following script loads the responses from a single subject
%   (subject LKC in Bonnen et al 2015) and performs a kalman filter analysis 
%   as well as a ccg analysis (see figures 7 and 4 in Bonnen et al 2015).

clear;
load data
% Set the starting point to zero
target = bsxfun(@minus, target, target(:,1));  
response =bsxfun(@minus, response, response(:,1));

colors = [    0.6980    0.0941    0.1686;...
    0.9373    0.5412    0.3843;...
    0.9922    0.8588    0.7804;...
    0.8196    0.8980    0.9412;...
    0.4039    0.6627    0.8118;...
    0.1294    0.4000    0.6745];

%% Kalman Filter Analysis

Q=1;        
lag=12;     
clip =60;   % Clipping the first second (i.e. the first 60 frames)

% Take off the first 'clip' frames.  Offset the target and response by the 'lag'
targetc = bsxfun(@minus, target(:,1+clip:end-lag), mean(target(:,1+clip:end-lag),2));  
responsec = bsxfun(@minus,response(:,lag+1+clip:end), mean(response(:,lag+1+clip:end),2));

% Set up the minimization
r0 = log(100);
warning off
opt.Display = 'Off';


% Find the conditions (sigma - width of the gaussian blob).
sigmas = flip(unique(sigma));

% Find Positional uncertainty with minimization and plot for each condition
f=figure;
subplot(121);
for i=1:length(sigmas)
    ind = sigma==sigmas(i);
    [r(i),fval(i)] = fminunc(@negLogLikelihoodr,r0,opt,Q,targetc(ind,:)',responsec(ind,:)');
    semilogx(sigmas(i),sqrt(r(i)),'.','LineWidth',3,'MarkerSize',45,'Color',colors(i,:))
    hold on;
end

% The fitting code returns the log value of the positional uncertainty
% estimate in order to force it to be positive.  We exponentiate to return
% it to (-Inf,Inf).
r = exp(r);

xlabel('blob width (sigma in arcmin)');
ylabel('positional uncertainty estimate (arcmin)');
title('LKC - Kalman Filter Analysis');

%% CCG Analysis

% Find velocities and zscore in preparation for calculating CCGs
targVel = zscore(diff(target'))';
respVel = zscore(diff(response'))';

% Calculate and plot CCGs
subplot(122);
for i=1:length(sigmas)
    ind = find(sigma==sigmas(i));
    for j=1:length(ind)
        ccg(j,:) = xcorr(respVel(ind(j),:),targVel(ind(j),:),60,'coeff');
    end
    
    avg(i,:) = mean(ccg);
    h(i) = plot((-60:60)/60,avg(i,:),'LineWidth',3,'Color',colors(i,:));
    hold on;
end

xlim([0,1])
xlabel('time (s)');
ylabel('correlation value');
title('LKC - Average CCGs');
legend(flip(h),num2str(flip(sigmas)))

%%
saveas(f,'example.fig');
