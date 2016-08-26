function [rstar, Zr, P] = rNr2rstar(rIN, NrIN, plot)
% Given a vector of counts rIN, and a corresponding vector of counts
% of counts (how many species occur r times) NrIN, compute and return the
% smoothed counts using the Simple Good-Turing (SGT) estimator as rstar,
% where rstar(k) is the SGT estimate of r(k) for all indices
% 1 <= k <= length(r). Also returns the middle products Zr (after
% truncating Nr for large r) and P (best fit for lnZr-lnr).
% If plot evaluates to true, then plot Nr and Zr against r in log-log scale, 
% and plot the line for P as well.
%
% Implementation follows the methods described in:
% Gale, W. and Sampson, G. (1995). "Good-Turing Frequency Estimation
% Without Tears," J. Quant. Linguistics 2, 24.
%
% Sample input and output:
% r = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 ...
%     27 28 29 31 33 39 41 46 47 50 52 53 55 57 60 74 84 108 109 177 400 1918];
% Nr = [-1 268 112 70 41 24 14 15 14 8 11 9 6 6 3 7 9 4 4 8 2 4 2 2 3 4 4 4 1 1 ...
%     2 1 3 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1];
% rstar = [0.84,1.35,2.27,3.23,4.19,5.17,6.15,7.14,8.13,9.12,10.11,11.11,12.1,...
%     13.1,14.09,15.09,16.09,17.09,18.08,19.08,20.08,21.08,22.07,23.07,24.07,...
%     25.07,26.07,27.07,28.07,30.07,32.06,38.06,40.06,45.06,46.06,49.05,51.05,...
%     52.05,54.05,56.05,59.05,73.05,83.05,107.04,108.04,176.04,399.04,1917.04];
% Zr = [268,112,70,41,24,14,15,14,8,11,9,6,6,3,7,9,4,4,8,2,4,2,2,3,4,4,4,1,0.7,...
%     1,0.3,0.8,0.3,0.3,0.5,0.4,1.3,0.7,0.5,0.4,0.1,0.1,0.1,0.1,0,0,0,0];
% lnr = [0,0.69,1.1,1.39,1.61,1.79,1.95,2.08,2.2,2.3,2.4,2.48,2.56,2.64,2.71,...
%     2.77,2.83,2.89,2.94,3,3.04,3.09,3.14,3.18,3.22,3.26,3.3,3.33,3.37,3.43,...
%     3.5,3.66,3.71,3.83,3.85,3.91,3.95,3.97,4.01,4.04,4.09,4.3,4.43,4.68,4.69,...
%     5.18,5.99,7.56];
% lnZr = [5.59,4.72,4.25,3.71,3.18,2.64,2.71,2.64,2.08,2.4,2.2,1.79,1.79,1.1,...
%     1.95,2.2,1.39,1.39,2.08,0.69,1.39,0.69,0.69,1.1,1.39,1.39,1.39,0,-0.41,...
%     0,-1.39,-0.29,-1.25,-1.1,-0.69,-0.92,0.29,-0.41,-0.69,-0.92,-2.14,-2.48,...
%     -2.83,-2.53,-3.54,-4.98,-6.77,-7.33];
% P = [-1.9646,6.6834];

len = length(rIN);
% pop the first term off
r = rIN(2:len);
Nr = NrIN(2:len);
len = len - 1;

if len >= 2
    % Compute Zr
    k = [r(2:len) 2*r(len)-r(len-1)];
    i = [0 r(1:len-1)];
    diffki = k - i;
    Zr = 2 .* Nr ./ diffki; % vectorized operation
    
    % Compute lnr
    lnr = log(r);
    
    % Compute lnZr
    lnZr = log(Zr);
    
    % Compute P for linear regression where
    % lnZr =(best fit)= P(1).*lnr + P(2) = lnSr
    P = polyfit(lnr, lnZr, 1); % polynomial with degree 1
    
    % Compute rstar
    rstar = zeros(1, len);
    useY = false;
    for j = 1:len
        % Compute y using Sr
        y = (r(j)+1)*S(r(j)+1, P)/S(r(j), P); % exp(polyval(P, log(r(j)+1)))/exp(polyval(P, log(r(j))));
        if useY
            rstar(j) = y;
        elseif j+1 > len || r(j+1) ~= r(j)+1 % x is incomputable
            rstar(j) = y;
            useY = true;
        else
            % Compute x using Nr
            x = (r(j)+1)*Nr(j+1)/Nr(j);
            % Check the difference between x and y
            d = abs(x - y);
            threshold = 1.96*(r(j)+1)/Nr(j)*sqrt(Nr(j+1)*(1+Nr(j+1)/Nr(j))); %TODO
            if d > threshold
                rstar(j) = x;
            else
                rstar(j) = y;
                useY = true;
            end
        end
    end
else % len <= 1
    Zr = Nr; P = [0 0];
    rstar = r;
end

if plot
    % plot Nr against r in log-log scale
    figure
    ax = gca; ax.XLimMode = 'manual'; ax.XLim = [-10^(-1) 10^2]; ax.YLim = [0, 10^3];
    loglog(r, Nr, 'g*')
    xlabel('r')
    title('Nr vs. r log-log scale')
    
    % plot Nr and Zr against r in log-log scale, plot P
    figure
    x1 = linspace(1,r(len)+1);
    y1 = S(x1, P);
    ax = gca; ax.XLimMode = 'manual'; ax.XLim = [-10^(-1) 10^2]; ax.YLim = [0, 10^3];
    loglog(r, Nr, 'g*', r, Zr, 'bo', x1, y1, 'r')
    xlabel('r')
    title('Nr & Zr vs. r log-log scale')
    legend('Nr','Zr')
end

end

function Sr = S(r, P)
lnr = log(r);
lnSr = polyval(P, lnr); % P(1).*lnr + P(2)
Sr = exp(lnSr);
end