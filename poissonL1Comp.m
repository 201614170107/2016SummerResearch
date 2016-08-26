function [sgtErrAvg, empErrAvg, lapErrAvg] = poissonL1Comp(lambda)
% Given lambda, generate a random sample under Poisson distribution of size
% 100 and estimate the true distributio from the sample using Simple 
% Good-Turing estimator, empirical estimator, and Laplace estimator. 
% Return the l1-distance averaged over 1000 trials as the error for each 
% estimator.

min = lambda*3/5; max = lambda*5/3;

% what if not narrow done the range?
% min = 1; max = 100; % in this case emprical estimator does the best

% Run 1000 trials
sgtErrs = zeros(1,1000);
empErrs = zeros(1,1000);
lapErrs = zeros(1,1000);
for trial = 1:1000
    % generate 100 integers under Poisson distribution
    sample = poissrnd(lambda, 1, 100);
    
    % calculate sgtProb and sgtErr
    [species, speciesR, sgtProb, ~, ~, ~, ~, ~] = mat2prob(sample, min, max, false);
    trueProb = poisspdf(species, lambda); % calculate the true probability
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