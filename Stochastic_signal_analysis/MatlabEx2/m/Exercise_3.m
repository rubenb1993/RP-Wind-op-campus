%% Exercise 3: Filter by approximating the correlation functions
% Goal: Approximate the auto- and cross correlation functions used in 
%       the Wiener-Hopf equations from the signals

%n = input('Order filter: ');
% Let n vary between the desired filter orders
n = [1 2 4 6];
Stdd2 = zeros(4,1);
W_tot2 = zeros(6,4);
for k = 1:length(n)
% Construct an n by (N-n+1) matrix V2 containing shifted versions of v2
V2 = zeros(n(k), N-n(k)+1);
for i=1:n(k)
    V2(i,:) = v2(n(k)+1-i:N+1-i);
end

% Use V2 and v1 to estimate Rv2 and Rv1v2
rv2e = zeros(n(k),1);
for i = 1:n(k)
    rv2e(i) = sum(V2(1,:).*V2(i,:));
end

% Put rv2e in a Toeplitz matrix like in Exercise 2
rv2e_ds = [rv2e(end:-1:1);rv2e(2:end)];
Rv2e = zeros(n(k),n(k));
for i = 1:n(k)
    for j = 1:n(k)
        Rv2e(i,j) = rv2e_ds(n(k)+j-i);
    end
end

Rv1v2e = zeros(n(k),1);
for i=1:n(k),
    Rv1v2e(i,1) = sum(x(n(k):end).*(V2(i,:).'));
end

% Calculate filter using Wiener-Hopf equations and reconstruct the signal
w = Rv2e\Rv1v2e;
v1e = filter(w,1,v2);
de = x - v1e;

% Save all the variables for the different values of n
std(d-de);
Stdd2(k) = std(d-de);
W_tot2(1:length(w),k) = w;
end