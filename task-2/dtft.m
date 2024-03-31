function X = dtft(x)
% Compute the DTFT
N = length(x); % Length of the sequence
n = 0:N-1; % Time indices
omega = linspace(-pi, pi, 1000); % Frequency vector (adjust as needed)
% Compute the DTFT using the formula
X = zeros(size(omega)); % Initialize DTFT
for k = 1:length(omega)
    X(k) = sum(x .* exp(-1i * omega(k) * n));
end
% Plot the magnitude and phase of the DTFT
figure;
subplot(2,1,1);
plot(omega, abs(X));
title('Magnitude of DTFT');
xlabel('\omega');
ylabel('|X(e^{j\omega})|');
subplot(2,1,2);
plot(omega, angle(X));
title('Phase of DTFT');
xlabel('\omega');
ylabel('\angle X(e^{j\omega})');
