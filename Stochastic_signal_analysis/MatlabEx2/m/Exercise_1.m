%% Initialization

close all
clear all

% Load uncorrupted signal d
% N is the number of samples (signal length)
load gong.mat;
[N,k]=size(y);
d = y;

% Generate zero-mean white noise sequence g with standard deviation 0.35
sg = 0.35;
g = sg*randn(N,1);
g = g - mean(g);

% Generate noise sequences v1 and v2
a1 = [1 -0.90]; b1 = [1 -.2];
a2 = [1 -0.95]; b2 = [1 -.3];
v1 = filter(b1,a1,g);
v2 = filter(b2,a2,g);

% Generate the corrupted signal x
x = d + v1;

%% Exercise 1: Determining the optimal FIR Wiener Filter
% Goal: reconstruct d from x and v2 by estimating v1 from v2

% Let n vary between the desired filter orders
n = [1 2 4 6];
Stdd = zeros(length(n),1);
W_tot = zeros(max(n),length(n));
for k = 1:length(n)
% First we determine Rv2 and Rv1v2 needed to set up the Wiener-Hopf

% Calculate the first two values of rv2 (i.e. rv2(0) and rv2(1))
% using eq (3.116) from Hayes
h = dimpulse(b2, a2, 20);
c(1,1) = b2(1)*conj(h(1)) + b2(2)*conj(h(2));
c(2,1) = b2(2)*conj(h(1));

 rv2 = zeros(200,1);
rv2(1) = (sg^2*c(1,1) - a2(2)*sg^(2)*c(2,1))/(1-a2(2)^2);
rv2(2) = sg^2*c(2,1) - a2(2)*rv2(1);

% Calculate the rest of rv2 until it becomes (almost) zero
% We only determine one side of the auto-correlation function and then 
% mirror, to get the double sided ACF (centered at index 200)
for i=3:200,
    rv2(i) = -1*a2(2)*rv2(i-1);
end
rv2_ds = [rv2(end:-1:1); rv2(2:end)];

% Next we determine the the cross-correlation function between v1 and v2
bb = conv(b1,a2);
aa = conv(b2,a1);
rv1v2_ds = filter(bb,aa,rv2_ds);

% Put rv2 and rv1v2 into matrix form for the Wiener-Hopf equations
Rv2 = zeros(n(k),n(k));
for i = 1:n(k)
    for j = 1:n(k)
        Rv2(i,j) = rv2_ds(200+j-i);
    end
end
Rv1v2 = zeros(n(k),1);
for i=1:n(k),
    Rv1v2(i,1) = rv1v2_ds(200+i-1);
end

% Solve for the optimal filter
W = Rv2\Rv1v2;
v1e = filter(W,1,v2);
de  = x - v1e;

% Save the necessary data necessary to evaluate the sound
W_tot(1:length(W),k) = W;
Stdd(k) = std(d-de);
end