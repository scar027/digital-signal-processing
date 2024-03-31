function x_t = d2a(x_n_enc,Fs,bit,min_x,max_x)
%D2A Digital to Analog Convertor using Decoding, Inverse Quantisation
%and Interpolation.
%   Detailed explanation goes here

% Decoding
s = 0;
r = 1;
for p = 1:length(x_n_enc)/bit
    for q = 1:8
        x_n_enc_mat(r,q) = x_n_enc(s + q);
    end
    r = r + 1;
    s = s + bit;
end
x_n_dec = bit2int((x_n_enc_mat).',bit);
figure
subplot(131),bar(x_n_dec, "w", "EdgeColor", "#0072BD");
title('Decoding')
xlabel('n');
ylabel('x(n)');
grid on;

% Inverse Quantisation
k = min_x;
l = max_x;
x_n_dec_inv_quant = x_n_dec.*((l-k)/(2.^(bit))) + k;
% Plotting the quantised signal
B = Fs/2; % Nyquist Rate Sampling
n = -B:1/Fs:B;
subplot(132), stem(n,x_n_dec_inv_quant); 
title('Inverse Quantisation')
xlabel('n');
ylabel('x(n)');
grid on;

% Interpolation
xq = -B:1/Fs:B;
vq1 = interp1(n,x_n_dec_inv_quant,xq);
subplot(133),plot(xq,vq1);
title('Interpolation')
xlabel('t');
ylabel('y(t)');
grid on;

x_t = vq1;
end

