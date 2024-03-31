function [x_n_enc, min_x_n, max_x_n] = a2d(signal,Fs,bit)
%A2D Analog to Digital Convertor using Sampling, Quantisation and Encoding
%   Detailed explanation goes here

% Sampling(t --> n/Fs)
B = Fs/2; % Nyquist Rate Sampling
n = -B:1/Fs:B;
x_n = signal(n);
figure
subplot(131), stem(n,x_n); % Plotting the sampled signal
title('Sampling')
xlabel('n');
ylabel('x(n)');
grid on;

% Quantisation
k = min(x_n);
l = max(x_n);
x_n_shift = x_n - k;
x_n_quant = floor((x_n_shift)/((l-k)/(2.^(bit)-1)));
subplot(132), bar(n,x_n_quant, "w", "EdgeColor", "#0072BD");
title('Quantisation')
xlabel('n');
ylabel('x(n) Discrete Levels');
grid on;

% Encoding
enc = (int2bit(x_n_quant,8)).';
r = 1;
for p = 1:length(enc)
    for q = 1:8
        b =  enc(p,q);
        enc_out(1,r) = b;
        r = r + 1;
    end
end
subplot(133),stairs([enc_out,enc_out(end)]);
title('Encoding')
xlabel('n');
ylabel('x(n) Binary Encoding');

x_n_enc = enc_out;
min_x_n = k;
max_x_n = l;
end