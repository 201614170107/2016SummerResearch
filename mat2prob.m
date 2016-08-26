function [species, speciesR, speciesProb, r, rstar, Nr, Zr, P] = mat2prob(sample, min, max, plot)
% assume all species are integers
% Given sample, a vector of species of different occurrences, compute
% and return species, speciesR, speciesProb by SGT, as well as middle
% products r, rstar, Nr, and P.
% Species not within the range from min to max are ignored.
% If plot evaluates to true, then plot SGT speciesProb against species.

% call helpers
[r, Nr, species, speciesR] = sample2rNr(sample, min, max, plot);
[rstar, Zr, P] = rNr2rstar(r, Nr, plot);

len = length(r);
N = sum(r(2:len) .* Nr(2:len));
Nstar = sum(rstar .* Nr(2:len));
if length(r) <= 1 || ~(r(2) == 1)
    error('Bad data to run sgt!!!')
end
p0 = Nr(2)/N;
rProb = (1-p0) .* rstar ./ Nstar;
speciesProb = zeros(size(species));

for i = 1:length(speciesProb)
    if speciesR(i) == 0
        speciesProb(i) = p0/Nr(1);
    else
        speciesProb(i) = sum( (r(2:len) == speciesR(i)) .* rProb);
    end
end

if plot
    % plot sgt probability of each species
    figure
    plot(species, speciesProb, '.')
    xlabel('species')
    ylabel('sgt probability')
    title('sgt probability of each species')
end

end