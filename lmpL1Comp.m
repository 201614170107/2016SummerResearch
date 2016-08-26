function [sgtErrAvg, empErrAvg, lapErrAvg] = lmpL1Comp(number)
% Given number denoting one of four representative cases, use the
% corresponding LMP data of size 10000 and empirical estimator to calculate
% the near-truth-probability; then partition the data into 100 samples of
% size 100 each and use Simple Good-Turing estimator, empirical estimator, 
% and Laplace estimator to estimate the probability distribution. Return the
% l1-distance averaged over 100 samples as the error for each estimator.

load('LMP_small.mat')
switch number
    case 1
        data = LMP_small(1).LMP(500,:);
        min = 76000; max = 76500; % TODO
    case 2
        data = LMP_small(1).LMP(3120,:);
        min = 61000; max = 61500; % TODO
    case 3
        data = LMP_small(2).LMP(1,:);
        min = 81800; max = 82100; % TODO
    case 4
        data = LMP_small(4).LMP(2000,:);
        min = 94200; max = 96200; % TODO
end

% what if we cannot narrow down the range?
% min = -156355; max = 210810; % in this case emprical estimator does the best

% rescale and round data to get integer samples
scale = 4;
samples = round(data.*10^scale);

% Compute near-truth probability using empirical estimator
[~, ~, ~, speciesR] = sample2rNr(samples, min, max, false);
N = length(samples);
trueProb = speciesR/N;

% Run 100 trials
sgtErrs = zeros(1,100);
empErrs = zeros(1,100);
lapErrs = zeros(1,100);
for trial = 1:100
    % take 100 out of the 10000 as a sample (sparse data)
    sample = samples(1+100*(trial-1):100*trial);
    
    % calculate sgtProb and sgtErr
    [~, speciesR, sgtProb, ~, ~, ~, ~, ~] = mat2prob(sample, min, max, false);
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

% de-rescale in the end (if necessary)
% species = species./(10^scale);

end