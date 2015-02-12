% SC3011TN - Stochastische Signaal Analyse
% Matlab Assignment 1
% January 2015

%% Constants
clear all
close all
clc

N = 1e5;
dt = 1e-8; % s
R = 1e-6; % m
kB = 1.38e-23; % J/K
T = 300; % K
eta = 1e-3; % Pa s
rho = 2.6e3; % kg/m^3
gamma = 6*pi*R*eta; % Pa m s

% Compute particle mass in kg (nb. particles are spheres)
m = 4/3 * pi * R^3 * rho; %[kg]

%% Exercise 1

% Expressions for the coefficients in terms of the constants given above
beta1 = -1*(2*m + gamma*dt)/(m + gamma*dt);
beta2 = m / (m + gamma*dt);
beta3 = sqrt(2*kB*T*gamma)*(dt)^(3/2)/(m + gamma*dt);

%% Exercise 2

% Initialize signal vector (x) and generate white noise samples vector (w)
x = zeros(N,1);
w = randn(N,1);

% Simulate the difference equation
for k = 3:N
    x(k) = beta3.*w(k) - beta1.*x(k-1) - beta2.*x(k-2);
end

% Plot the result as a function of time
plot(x);
xlabel('time (s)');
ylabel('displacement (m)');

%% Exercise 3

% Initialize number of simulation (L) and the signal and noise matrices (x and w)
L = 30;
x = zeros(N,L);
w = randn(N,L);

% Simulate the difference equation
for l = 1:L
    for k = 3:N
        x(k,l) = beta3.*w(k,l) - beta1.*x(k-1,l) - beta2.*x(k-2,l);
    end
end

% Plot the results
plot(x);
xlabel('time (s)');
ylabel('displacement (m)');


%% Exercise 4

% NB: Make sure your computer has enough memory / try small L values first!

% Initialize the number of simulation (vector of L values)
% and the number of samples in each simulation (vector of H values)
Ls = [30 300 3000];
Hs = [1e3 1e4 1e5];

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
        x(k,l) = beta3.*w(k,l) - beta1.*x(k-1,l) - beta2.*x(k-2,l);
    end
end
figure
% Take different subsets of the data for each combination of the 
% L and h parameters and plot as a histogram 
for i = 1:length(Ls)
    L = Ls(i);
        
    for j = 1:length(Hs)
        h = Hs(j);
        
        % Select appropriate subplot
        subplot(length(Ls), length(Hs), (i-1)*length(Hs)+j);
        
        % Plot the histogram
        hist(x(h,(1:L)),sqrt(L));        
        xlim([min(min(x)) max(max(x))]);
        xlabel('displacement (m)');
        title(sprintf('L=%d, h=%d', L, h));
        
    end
end

%% Exercise 5
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

