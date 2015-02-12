Kappa = zeros(3,3);
Sigma = zeros(3,3);
for i = 1:length(Ls)
    L = Ls(i);
        
    for j = 1:length(Hs)
        h = Hs(j);       
        
        % Repeat hist command to get the data
        [counts, centers] = hist(x(h,(1:L)),sqrt(L));
        
        % Calculate initial estimates of sigma and kappa
        sigma0 = std(x(h,(1:L)));
        kappa0 = counts(round(length(counts)/2));
        
        % Define the objective function which we'll try to minimize
        % p is a vector with fit parameters:
        %     p(1) standard deviation (sigma)
        %     p(2) scale factor (kappa)
        fobj = @(p) sum( (counts - p(2)*exp(-centers.^2/2/p(1)^2)).^2 );
        
        % Fit the Gaussian through the histogram data
        p_opt = fminsearch(fobj, [sigma0, kappa0]);
        sigma = p_opt(1);
        kappa = p_opt(2);
        Sigma(i,j) = sigma;
        Kappa(i,j) = kappa;
        % Select appropriate subplot
        subplot(length(Ls), length(Hs), (i-1)*length(Hs)+j);
                 
        % Plot the Gaussian on top of the histogram
        hold all;
        xg = linspace(min(min(x)), max(max(x)), 1000);
        yg = kappa.*exp(-(xg.^2)./(2.*sigma.^2));
        plot(xg, yg, 'r');
        xlabel(sprintf('sigma = %.2e (m)', sigma));
        hold off;        
        
    end
end    