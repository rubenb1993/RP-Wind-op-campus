% Initialize number of simulation (L) and the signal and noise matrices (x and w)
L = 30;
x = zeros(N,L);
w = randn(N,L);

% Simulate the difference equation
for l = 1:L
    for k = 3:N
        x(k,l) = - beta1*x(k-1,l) - beta2*x(k-2,l) + beta3*w(k,l);
    end
end

% Plot the results
time = (0:dt:dt*(N-1));
plot(time,x);
xlabel('time (s)');
ylabel('displacement (m)');