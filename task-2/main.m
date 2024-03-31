% Program for Sampling, Reconstruction & DTFT Analysis 

%%
clc; clear; close all;
%% TASK 1-A-1)

Fs = 10;
B = Fs/2;
t = -B:1/Fs:B;
X_a = t - t + 1; % X_a(F)

figure
subplot(121), plot(t, X_a) % Plotting the given continuous signal
title('X_a(F)')
xlabel('F');
ylabel('X_a(F)');
grid on;

% Generating the signal x_a(t) from X_a(F) using IDFT
x_a = sinc(t);
subplot(122), plot(t, x_a) % Plotting the given continuous signal
title('x_a(t)')
xlabel('t');
ylabel('x_a(t)');
grid on;

% A/D
[x_n_enc, min_x, max_x] = a2d(@sinc,Fs,8);

% x^2(n)
x_n_enc_sq = x_n_enc.^2;
figure
stairs([x_n_enc_sq,x_n_enc_sq(end)]);
title('x^2(n)')
xlabel('n');
ylabel('x^2(n)');


% D/A
x_a_dec = d2a(x_n_enc_sq,Fs,8,min_x,max_x);

%% TASK-1-A-2)

Fs = 10;
B = Fs/2;
t = -B:1/Fs:B;

% x^2(n)
figure
s_a = (sinc(t)).^2;
plot(t, s_a) % Plotting the continuous signal after squaring
title('s_a(t)')
xlabel('t');
ylabel('s_a(t)');
grid on;

% A/D
[s_n_enc, min_s, max_s] = a2d(@(x) (sinc(x)).^2,Fs,8);

% D/A
s_a_dec = d2a(s_n_enc,Fs,8,min_s,max_s);

%% TASK 1-B-1)

Fs = 50; % Fs = 50Hz
B = Fs/2;
Fo = 20;

t = -B:1/Fs:B;
x2_a = cos(2*pi*Fo*t*(1/Fs));

figure
plot(t, x2_a) % Plotting the given continuous signal
title('x_2_a(t)')
xlabel('t');
ylabel('x_2_a(t)');
grid on;

% A/D
[x2_a_n_enc, min_x2a, max_x2a] = a2d(@(x) cos(2*pi*Fo/Fs*x),Fs,8);

% x^2(n)
x2_a_n_enc_sq = x2_a_n_enc.^2;
figure
stairs([x2_a_n_enc_sq,x2_a_n_enc_sq(end)]);
title('x^2(n)')
xlabel('n');
ylabel('x^2(n)');


% D/A
x2_a_dec = d2a(x2_a_n_enc_sq,Fs,8,min_x2a,max_x2a);

% x^2(n)
figure
s2_a_a = (x2_a).^2;
plot(t, s2_a_a) % Plotting the continuous signal after squaring
title('s_2_a(t)')
xlabel('t');
ylabel('s_2_a(t)');
grid on;

% A/D
[s2_a_n_enc, min_s2_a, max_s2_a] = a2d(@(x) (cos(2*pi*Fo/Fs*x)).^2,Fs,8);

% D/A
s2_a_dec = d2a(s2_a_n_enc,Fs,8,min_s2_a,max_s2_a);

%% TASK 1-B-2)

Fs = 30; % Fs = 30Hz
B = Fs/2;
Fo = 20;

t = -B:1/Fs:B;
x2_b = cos(2*pi*Fo*t*(1/Fs));

figure
plot(t, x2_b) % Plotting the given continuous signal
title('x_2_b_a(t)')
xlabel('t');
ylabel('x_2_b_a(t)');
grid on;

% A/D
[x2_b_n_enc, min_x2b, max_x2b] = a2d(@(x) cos(2*pi*Fo/Fs*x),Fs,8);

% x^2(n)
x2_b_n_enc_sq = x2_b_n_enc.^2;
figure
stairs([x2_b_n_enc_sq,x2_b_n_enc_sq(end)]);
title('x^2(n)')
xlabel('n');
ylabel('x^2(n)');


% D/A
x2_b_dec = d2a(x2_b_n_enc_sq,Fs,8,min_x2b,max_x2b);

% x^2(n)
figure
s2_b_a = (x2_b).^2;
plot(t, s2_b_a) % Plotting the continuous signal after squaring
title('s_2_b_a(t)')
xlabel('t');
ylabel('s_2_b_a(t)');
grid on;

% A/D
[s2_b_n_enc, min_s2_b, max_s2_b] = a2d(@(x) (cos(2*pi*Fo/Fs*x)).^2,Fs,8);

% D/A
x2_b_dec = d2a(s2_b_n_enc,Fs,8,min_s2_b,max_s2_b);
