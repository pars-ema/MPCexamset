% Scalar Standard Brownian Motion = Standard Wiener Process
Ns = 10; % Number of realizations
T = 10; % Final time
N = 1000; % Number of time steps
seed = 100; % Seed for reproducibility

% Realization of Ns Standard Brownian Motions
rng(seed);
dt = T/N;
dW = sqrt(dt)*randn(Ns,N);
W = cumsum(dW,2);


% Time vector
t = linspace(0, T, N);

% Plot all realizations
figure;
plot(t, W'); % transpose so each row is one trajectory
xlabel('t');
ylabel('W(t)');
title('Standard Brownian Motion Realizations');
grid on;

[W,Tw,dW] = ScalarStdWienerProcess(10,1000,2,100)
figure
plot(Tw,W, "linewidth",2);
hold on

[Wmean,sW,Wmeanp2sW,Wmeanm2sW]=ScalarSampleMeanStdVar(W)
plot(t, sW(:, 1:end-1))

figure
T = 10;
N = 1000;
Ns = 10;
seed = 100;

x0 = 10;
lambda = -0.5;
sigma = 1.0;

[W,Tw,dW]=ScalarStdWienerProcess(T,N,Ns,seed);

X = zeros(size(W));
for i=1:Ns;
    X(i,1) = x0;
    for k=1:N
    dt = Tw(k+1)-Tw(k);
    X(i,k+1) = X(i,k) + lambda*X(i,k)*dt + sigma*dW(i,k);
    end
end

plot(Tw, X)