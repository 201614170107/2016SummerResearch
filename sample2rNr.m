function [r, Nr, species, speciesR] = sample2rNr(sample, min, max, plot)
% Given sample, a vector of species of different occurrences, compute
% and return species, a vector of non-repeating species, and speciesR,
% their corresponding numbers of occurrences (counts), as well as r, 
% a vector of ascending counts, and Nr, the counts of counts (number of 
% speicies that occur r times in the sample).
% Species not within the range from min to max are ignored.
% If plot evaluates to true, then plot Nr against r and speciesR against
% species.

% Sample input and output:
% v = [0.1 0.4, 0.25, 0.75, 0.76, 0.77, 0.78];
% min = 1; max = 8; scale = 10;
% v = [1 4 3 8 8 8 8]
% species =  [1 2 3 4 5 6 7 8]
% speciesR = [1 0 1 1 0 0 0 4]
% r =  [0 1 4]
% Nr = [4 3 1]

% species-speciesR pairs.
species = min:max;
speciesR = zeros(size(species));
r = 0;
for k = 1:length(species)
    s = sum(sample == species(k));
    if s ~= 0
        speciesR(k) = s;
    end
    if sum(r==s) == 0
        r = [r s];
    end
end

% r-Nr pairs.
r = sort(r);
Nr = zeros(size(r));
for k = 1:length(r)
    Nr(k) = sum(speciesR == r(k));
end

if plot
    % plot species count against species
    figure
    hold on
    bar(species, speciesR)
    grid on
    xlabel('species')
    ylabel('species count')
    title('species count vs. species')
    hold off
    
    % plot Nr against r
    figure
    hold on
    bar(r, Nr)
    grid on
    xlabel('r')
    ylabel('Nr')
    title('Nr (count of count) vs. r (count)')
    hold off
end

end