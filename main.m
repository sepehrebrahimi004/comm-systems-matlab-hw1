%MAIN  Communication Systems - MATLAB homework 1: Fourier transforms with the
%      FFT, direct DFT vs. FFT, the Hilbert transform and Fourier series.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

% Continuous-time Fourier transform approximated by the FFT:
%   X(f_k) ~ dt * sum_n x(t_n) e^{-j 2 pi f_k t_n}
ctft = @(x, dt) fftshift(fft(x)) * dt;
faxis = @(N, dt) (-floor(N/2):ceil(N/2)-1) / (N*dt);

%% Q1  Spectra of three signals
dt = 0.01;
sigs = {
    'x_1(t) = \Sigma_{n=0}^{5} e^{-(t-2n)} rect((t-2n)/2)', -5:dt:15-dt, ...
        @(t) sum(cell2mat(arrayfun(@(n) exp(-(t-2*n)) .* (abs(t-2*n) <= 1), (0:5)', 'UniformOutput', false)), 1)
    'x_2(t) = ln(t) rect((t-10)/10)', 0:dt:20-dt, @(t) log(max(t, eps)) .* (abs(t-10) <= 5)
    'x_3(t) = sinc^4(t) rect(t/10)', -20:dt:20-dt, @(t) sinc(t).^4 .* (abs(t) <= 5)
};
fig = figure('Visible', 'off', 'Position', [100 100 1100 750]);
for k = 1:3
    [name, t, fx] = sigs{k, :};
    x = fx(t);
    X = ctft(x, dt);  f = faxis(numel(t), dt);
    E_t = sum(abs(x).^2)*dt;  E_f = sum(abs(X).^2)/(numel(t)*dt);
    fprintf('Q1  %-48s energy %.4f (time) vs %.4f (frequency, Parseval)\n', name, E_t, E_f);
    subplot(3, 2, 2*k-1); plot(t, x); grid on; title(name); xlabel('t (s)');
    subplot(3, 2, 2*k);   plot(f, abs(X)); grid on; xlim([-5 5]); title('|X(f)|'); xlabel('f (Hz)');
end
% x3: sinc^4 has a compactly supported transform (Lambda * Lambda, width 4)
save_figure(fig, 'q1_spectra');

%% Q2  Direct DFT (O(N^2)) vs. FFT (O(N log N))
t = 0:dt:20-dt;  x = sigs{2, 3}(t);  N = numel(t);
f = faxis(N, dt);
tic; X1 = zeros(1, N);
for i = 1:N
    X1(i) = sum(x .* exp(-1j*2*pi*f(i)*t)) * dt;           % note the 2*pi
end
t_direct = toc;
tic; X2 = ctft(x, dt); t_fft = toc;
fprintf('Q2  N = %d: direct %.3f s, FFT %.5f s (%.0fx faster), max |difference| = %.1e\n', ...
        N, t_direct, t_fft, t_direct/t_fft, max(abs(X1 - X2)));

%% Q3  Hilbert transform: hilbert() vs. convolution with 1/(pi t)
M = 2000;
t = (-M/2:M/2-1)*(2/M);                                   % one period of length 2, no sample at t = 0
x = sin(pi*t).^3 .* cos(pi*t).^2 + sin(2*pi*t).^2 .* cos(3*pi*t);
xh1 = imag(hilbert(x));                                   % FFT-based (exact for periodic signals)
% periodic convolution with the Hilbert kernel 1/(pi t) (principal value):
dtt = t(2) - t(1);
kernel = @(tau) 1 ./ (2*tan(pi*tau/2));                   % periodised 1/(pi t) for period 2
xh2 = zeros(size(x));
for i = 1:M
    tau = t(i) - t;  tau(i) = NaN;                        % skip the singular point
    xh2(i) = sum(x .* kernel(tau), 'omitnan') * dtt;
end
fprintf('Q3  max |hilbert() - convolution| = %.2e;  <x, x_hat> = %.2e (orthogonal)\n', ...
        max(abs(xh1 - xh2)), sum(x .* xh1)*dtt);
fig = figure('Visible', 'off', 'Position', [100 100 900 320]);
plot(t, x, 'k', t, xh1, 'b', t, xh2, 'r--', 'LineWidth', 1.2); grid on;
legend('x(t)', 'Hilbert via hilbert()', 'Hilbert via convolution with 1/(\pit)');
xlabel('t'); title('Hilbert transform computed two ways');
save_figure(fig, 'q3_hilbert');

%% Q4  Exponential vs. trigonometric Fourier series of a period-4 signal
N = 1000;  t = linspace(-2, 2, N+1);  t = t(1:end-1);
x = [zeros(1, N/4) ones(1, N/4) linspace(0, 1, N/4) linspace(1, 0, N/4)];
K = 50;  dt = t(2) - t(1);
c = arrayfun(@(n) sum(x .* exp(-1j*n*pi/2*t))*dt/4, -K:K);
x_exp = real(sum(c.' .* exp(1j*pi/2*(-K:K).'*t), 1));
a0 = sum(x)*dt/4;
an = arrayfun(@(n) sum(x .* cos(n*pi/2*t))*dt/2, 1:K);
bn = arrayfun(@(n) sum(x .* sin(n*pi/2*t))*dt/2, 1:K);
x_trig = a0 + an*cos(pi/2*(1:K).'*t) + bn*sin(pi/2*(1:K).'*t);
fprintf('Q4  max |exponential - trigonometric series| = %.1e  (c_n = (a_n - j b_n)/2)\n', max(abs(x_exp - x_trig)));
fig = figure('Visible', 'off', 'Position', [100 100 900 320]);
plot(t, x, 'k', t, x_exp, 'b', t, x_trig, 'r--', 'LineWidth', 1.2); grid on;
legend('x(t)', 'exponential series (K = 50)', 'trigonometric series (K = 50)', 'Location', 'northwest');
title('The two forms of the Fourier series are identical'); xlabel('t');
save_figure(fig, 'q4_fourier_series');
