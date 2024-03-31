% Program for generation of elementary signals

%% Generation of Elementary Signals
clc; clear; close all;
% Unit impulse response
n = 1:10;
u = [zeros(1,4) 1 zeros(1,5)];
figure
subplot(331),stem(n,u);
title('Unit Impulse Response');
xlabel('Time (n)---->');
ylabel('u[n]');
grid;

% code for unit sample response
s = [ones(1,10)];
subplot(332), stem(n,s);
title('Unit Sample Response');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% code for ramp signal
s = [ones(1,10)];
subplot(333),stem(n,n.*s);
title('Ramp Response');
xlabel('Time (n)---->');
ylabel('r[n]');
grid;

% code for decaying exponential signal
a = n*(-0.25); e = exp(a);
subplot(334),stem(n,e);
title('Decaying Exponential Signal');
xlabel('Time (n)---->');
ylabel('e[n]');
grid;

% code for sinewave signal
a = 1:10; pi = 3.14; s = sin(2*pi*a/10);
subplot(335),stem(n,s);
title('Sine-wave Signal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% code for sinewave signal - time reversal
a = 1:10; pi = 3.14; s = sin(2*pi*(-a)/10);
subplot(336),stem(n,s);
title('Sine-wave Signal - Time Reversal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% code for cosine signal
a = 1:10; pi = 3.14; s = cos(2*pi*a/10);
subplot(337),stem(n,s);
title('Cosine Signal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% code for cosine signal - time reversal
a = 1:10; pi = 3.14; s = cos(2*pi*(-a)/10);
subplot(338),stem(n,s);
title('Cosine Signal - Time Reversal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

%% TASK 1

% 1) code for sinewave signal
figure
n = 1:10;
a = 1:10; pi = 3.14; s = sin(2*pi*a/10);
subplot(331),stem(n,s);
title('Q.1) Sinusoidal Signal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% 1.a) code for sinewave time shifting
a = 1:10; pi = 3.14; k = 2; s = sin(2*pi*(a-k)/10); %k is the delay
subplot(332),stem(n,s);
title('Q.1)a) Time Shifting');
xlabel('Time (n)---->');
ylabel('s[n-k]');
grid;

% 1.b) code for sinewave time reversal
a = 1:10; pi = 3.14; s = sin(2*pi*(-a)/10);
subplot(333),stem(n,s);
title('Q.1)b) Time Reversal');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% 1.c) code for verifying sinewave TD[FD] = FD[TD]
a = 1:10; pi = 3.14; k = 2; s = sin(2*pi*((-a-k)/10)); %TD[FD]
subplot(334),stem(n,s);
title('Q.1)c) TD[FD]');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

s = sin(2*pi*((-(a)-k)/10)); %FD[TD]
subplot(335),stem(n,s);
title('Q.1)c) FD[TD]');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% 1.d) code for convolution
a = 1:10; pi = 3.14; s1 = sin(2*pi*(-a)/10);
s2 = sin(2*pi*a/10);
for i = 1:length(s1)
    for j = 1:length(s2)
        y(i,i+j-1) = s1(i)*s2(j);
    end
end
disp(y)
s = sum(y);
subplot(335);plot(s);
title('Q.1)d) Convolution');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% 1.e) code for correlation
a = 1:10; pi = 3.14; s1 = sin(2*pi*a/10);
s2 = sin(2*pi*a/10);
for i = 1:length(s1)
    for j = 1:length(s2)
        y(i,i+j-1) = s1(i)*s2(j);
    end
end
disp(y)
s = sum(y);
subplot(336);plot(s);
title('Q.1)e) Correlation');
xlabel('Time (n)---->');
ylabel('s[n]');
grid;

% 2)Verifying whether the given system is linear/non-linear stable/unstable
% System Input
figure
n = 1: 100;
u1 = [ones(1,100)];
u2 = [zeros(1,2) ones(1,98)];
u3 = u1 + u2;
subplot(331),stem(n,u3);
title('Q.2) System Input');
xlabel('Time (n)---->');
ylabel('u3[n]');
grid;

% System Output
a = 0;
for i = 1:100
    y1(i) = a.*a + u3(i);
    a = y1(i);
end
subplot(332),stem(n,y1);
title('Q.2) System Output');
xlabel('Time (n)---->');
ylabel('y1[n]');
grid;

% The output is unbounded if the input is unbounded. Hence Unstable

% Q.3)
% a)
figure
n = -5:5;
dir1 = [zeros(1,3) 1 zeros(1,7)];
dir2 = [zeros(1,9) 1 0];
x3a = 2*dir1 - dir2;
subplot(221),stem(n,x3a);
title('Q.3)a)');
xlabel('Time (n)---->');
ylabel('x3a[n]');
grid;

% b)
n = 1:20;
u_1 = [ones(1,20)];
u_2 = [zeros(1,9) ones(1,11)];
u_3 = [zeros(1,19) 1];
for i = 1:20
    x3b(i) = i*(u_1(i)-u_2(i))+10*exp(-0.3*(i-10))*(u_2(i)-u_3(i));
end
subplot(222),stem(n,x3b);
title('Q.3)b)');
xlabel('Time (n)---->');
ylabel('x3b[n]');
grid;

% c)
n = 1:50;
w = randn(1, 50);
for i = 1:50
    x3c(i) = cos(0.04*pi*i)+0.2*w(i);
end
subplot(223),stem(n,x3c);
title('Q.3)c)');
xlabel('Time (n)---->');
ylabel('x3c[n]');
grid;

% Q.4)
% a)
figure
a = [1 -1 0.9];
b = [1 0 1];
h = impz(b, a, 100);
subplot(221),stem(h);
title('Q.4)a)');
xlabel('Time (n)---->')
ylabel('Impulse Response h(n)')
grid on;

% b)
a = [1 -1 0.9];
b = [1 0 1];
step_res = stepz(b, a, 100);
subplot(222),plot(step_res);
title('Q.4)b)');
xlabel('Time (n)---->')
ylabel('s[n]')
grid on;

% c)
% The graph is converging hence the system is stable

% 5) Linear and Time Invariant System defined by Difference equation
% a)
% for the given difference equation, we can take Z-Transform on both sides 
% and get the Transfer function as:
% H(z) = (1 + 2z^-1 + z^-3)/(1-0.5z^-1+0.25z^-2)
% if the roots of the equation in the denominator lie within the unit
% circle in the z-plane then the system is said to be stable.
figure
coeff = [1, -0.5, 0.25];
roots = roots(coeff);
magnitudes = abs(roots);
if all(magnitudes < 1)
    disp('Q.5)a) The system is stable')
else
    disp('Q.5)b) The system is unstable')
end
% The system is stable as this method returns stable in the command window 
% output.

%b)
% Transfer function: H(z) = (1 + 2z^-1 + z^-3)/(1-0.5z^-1+0.25z^-2)
b = [1 2 0 1];
a = [1 -0.5 0.25];
h = impz(b, a, 100);
subplot(111),stem(h);
title('Q.5)b)');
xlabel('Time (n)---->')
ylabel('Impulse Response h(n)')
grid on;
% The impulse response decays over time hence the system is stable.

% 6) Simple digital differentiator
% a) Rectangular Pulse
figure
n = 1:20;
u_n = [ones(1,20)];
u_n10 = [zeros(1,9) ones(1,11)];
u_n20 = [zeros(1,19) 1];
x6a = 5*(u_n - u_n20);
subplot(331),stem(n,x6a);
title('Q.6)a)');
xlabel('Time (n)---->')
ylabel('Q.6)a)')
grid on;

% Implementing the differentiator
for i = 2:20
    y6a(i) = x6a(i)-x6a(i-1);
end
subplot(332),stem(n,y6a);
title('Q.6)a) Differentiator Output');
xlabel('Time (n)---->');
ylabel('y6a[n]');
grid;

% b) Triangular Pulse
x6b = n.*(u_n - u_n10)+(20-n).*(u_n10 - u_n20);
subplot(334),stem(n,x6b);
title('Q.6)b)');
xlabel('Time (n)---->')
ylabel('Q.6)b)')
grid on;

% Implementing the differentiator
for i = 2:20
    y6b(i) = x6b(i)-x6b(i-1);
end
subplot(335),stem(n,y6b);
title('Q.6)b) Differentiator Output');
xlabel('Time (n)---->');
ylabel('y6b[n]');
grid on;

% c) Sinusoidal Pulse
n = 1:100;
u_n = [ones(1,100)];
u_n100 = [zeros(1,99) 1];
x6c = sin(pi*n/25).*(u_n - u_n100);
subplot(337),stem(n,x6c);
title('Q.6)c)');
xlabel('Time (n)---->')
ylabel('Q.6)c)')
grid on;

% Implementing the differentiator
for i = 2:100
    y6c(i) = x6c(i)-x6c(i-1);
end
subplot(338),stem(n,y6c);
title('Q.6)c) Differentiator Output');
xlabel('Time (n)---->');
ylabel('y6c[n]');
grid on;
