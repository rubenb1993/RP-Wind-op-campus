%% Constants

N = 1e5;
dt = 1e-8; % s
R = 1e-6; % m
kB = 1.38e-23; % J/K
T = 300; % K
eta = 1e-3; % Pa s
rho = 2.6e3; % kg/m^3
gamma = 6*pi*R*eta; % Pa m s

% Compute particle mass in kg (nb. particles are spheres)
m = (4/3)*pi*R^3*rho;

% Expressions for the coefficients in terms of the constants given above
beta1 = -(2*m + gamma*dt)/(m + gamma*dt);
beta2 = m/(m + gamma*dt);
beta3 = sqrt(2*kB*T*gamma/dt)/(m/dt^2 + gamma/dt);

% Initialize signal vector (x) and generate white noise samples vector (w)
N2 = 1e3;
x = zeros(N2,1);
w = randn(N2,1);

% Simulate the difference equation
for k = 3:N2
    x(k) = - beta1*x(k-1) - beta2*x(k-2) + beta3*w(k);
end

% Plot the result as a function of time
time = (0:dt:dt*(N2-1));
plot(time,x);
xlabel('time (s)');
ylabel('displacement (m)');