% Initialize the number of simulation (vector of L values)
% and the number of samples in each simulation (vector of H values)
Ls = [30,300,3000];
Hs = [10^3,10^4,10^5];

% Take L the maximum number of simulations we need to compare
% (we can take subsets for the lower values of L)
L = max(Ls);

% Note that N (defined above) is equal to the maximum h

% Initialize signal and noise matrices (x and w)
x = zeros(N,L);
w = randn(N,L);

% Simulate the difference equation L times
for l = 1:L
    for k = 3:N
        x(k,l) = - beta1*x(k-1,l) - beta2*x(k-2,l) + beta3*w(k,l);
    end
end

% Take different subsets of the data for each combination of the 
% L and h parameters and plot as a histogram 
for i = 1:length(Ls)
    L = Ls(i);
        
    for j = 1:length(Hs)
        h = Hs(j);
        
        % Select appropriate subplot
        subplot(length(Ls), length(Hs), (i-1)*length(Hs)+j);
        
        % Plot the histogram
        hist(x(h,1:L),sqrt(L));        
        xlim([min(min(x)) max(max(x))]);
        xlabel('displacement (m)');
        title(sprintf('L=%d, h=%d', L, h));
        
    end
end