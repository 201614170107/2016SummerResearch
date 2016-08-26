function [sgtErrAvg, empErrAvg, lapErrAvg] = normalL1Comp(sigma)
% Given mu and sigma, generate a random sample under normal distribution of 
% size 100 and estimate the true distributio from the sample using Simple 
% Good-Turing estimator, empirical estimator, and Laplace estimator. 
% Return the l1-distance averaged over 1000 trials as the error for each 
% estimator.

mu = 0;
min = -3*sigma; max = 3*sigma;

% Run 1000 trials
sgtErrs = zeros(1,1000);
empErrs = zeros(1,1000);
lapErrs = zeros(1,1000);
for trial = 1:1000
    % generate 100 integers under normal distribution
    sample = normrnd(mu, sigma, 1, 100);
    sample = round(sample);
    
    % calculate sgtProb and sgtErr
    [species, speciesR, sgtProb, ~, ~, ~, ~, ~] = mat2prob(sample, min, max, false);
    trueProb = normpdf(species, mu, sigma); % calculate the true probability
    sgtErr = sum(abs(sgtProb - trueProb));
    sgtErrs(trial) = sgtErr;
    
    N = length(sample); % sample size, should be 100
    
    % calculate empProb and empErr
    empProb = speciesR/N;
    empErr = sum(abs(empProb - trueProb));
    empErrs(trial) = empErr;
    
    % calculate lapProb and lapErr
    lapProb = (speciesR+1)/(N+length(trueProb));
    lapErr = sum(abs(lapProb - trueProb));
    lapErrs(trial) = lapErr;
end

% take the avg error of 100 trials
sgtErrAvg = mean(sgtErrs);
empErrAvg = mean(empErrs);
lapErrAvg = mean(lapErrs);

end