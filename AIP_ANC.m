% Anex to the Paper:
% Low Delay and Low Cost Sigma-Delta Adaptive Controller for Active Noise Control
% Paulo Lopes

% simulation parameters
L = 10000;        % main simulation samples (at 44 kHz)
Li = 400;         % initial wait to fill the filter buffers
NW = 32;          % controler filter size
NS = 32;          % secondary path filter size
fs = 44100;       % lower sampling frequency
K = 32;           % oversample

rng(7283723);

load AAF.mat
SAAF = AAF;      % still uses the sharp anti-alising filter in some cases

load AAF3.mat
% AAF3.mat: transition 1 to 22050 kHz, attenuation 50 dB riplie 1 dB
% order 114
N_AAF = length(AAF);
AAF1 = AAF(1:K:end);
AAF1 = conv(AAF1,AAF1);

w_delays = [3.7, 2.2, 5.3];
w_amplitudes = [1.1, -0.7, 0.5]/2;
wop = sinc((0:NW*K-1)'-w_delays*K)*w_amplitudes';

s_delays = [5.2, 7.2, 3.7];
s_amplitudes = [-1.3, 0.9, -0.5];
sp = sinc((0:NS*K-1)'-s_delays*K)*s_amplitudes';
sp1 = conv(sp, ones(K,1)/K);
sp1 = sp1(K/2:end);
sh = conv(conv(sp1, AAF),AAF);
sh = [zeros(K,1); sh]; % there is an additional delay
sh = K*sh(1:K:end);
sh = sh(1:NS);

% simulation signals
up = randn(K*L,1); % up=conv(up, ones(4,1)/4); up=up(1:end-3);
yp = zeros(K*L,1);
dp = zeros(K*L,1);
ep = zeros(K*L,1);

% algorithm signals
u0 = zeros(K*L,1);
yq = zeros(K*L,1);
y0 = zeros(K*L,1);
e0 = zeros(K*L,1);

d = zeros(L,1);
u = zeros(L,1);
u1 = zeros(L,1);
e = zeros(L,1);
y = zeros(L,1);
dh = zeros(L,1);
eh = zeros(L,1);

% algorithm parameters
mu = 0.2;

% algorithm inicilizations
w = zeros(NW,1);
n = Li;

for n0=Li*K:L*K
    % discrete
    u0(n0) = up(n0);
    y0(n0) = y(n-1);

    % physical
    yq(n0) = wop'*up(n0:-1:n0-NW*K+1);
    yp(n0) = AAF*y0(n0:-1:n0-length(AAF)+1);
    dp(n0) = sp'*yq(n0:-1:n0-NW*K+1);
    ep(n0) = dp(n0) - sp'*yp(n0:-1:n0-NW*K+1);
    
    % discrete
    e0(n0) = ep(n0);
    
    if mod(n0, K)==0
        n = n0/K;
        u(n) = AAF*u0(n0:-1:n0-length(AAF)+1);
        e(n) = AAF*e0(n0:-1:n0-length(AAF)+1);
        
        u1(n) = sh'*u(n:-1:n-NS+1);

        y(n) = w'*u(n:-1:n-NW+1);
        dh(n) = e(n) + sh'*y(n:-1:n-NS+1);
        eh(n) = dh(n) - w'*u1(n:-1:n-NW+1);

        u1v = u1(n:-1:n-NW+1);
        w = w + mu*eh(n)*u1v/(u1v'*u1v+1e-1);
    end
end

figure(1);
ep_filtered = conv(ep, SAAF);
ep_down_sampled = ep_filtered(K:K:end);
plot((0:length(e)-1)/fs*1e3,e); hold on;
plot((0:length(ep_down_sampled)-1)/fs*1e3,ep_down_sampled); hold off;
legend('e','e_p')
set(gca,'XLim',[0 L/fs*1e3]);
xlabel('Time (ms)'); ylabel('Level');
set(gcf,'Name','External and internal residual error in time');
saveas(gcf, '../results/fig4.png')

figure(2);
[Pd, f_pd] = pwelch(dp(end/2:end), [], [], K*128, K*fs);
[Pe, f_pe] = pwelch(ep(end/2:end), [], [], K*182, K*fs);
plot(f_pd/1e3, 10*log10(Pd)); hold on;
plot(f_pe/1e3, 10*log10(Pe)); hold off;
set(gca, 'XLim', [0 fs/1e3]);
grid on;
legend('ANC off', 'ANC on','Location','southwest');
ylabel('Noise power spectral density (dB)');
xlabel('Frequency (kHz)');
set(gcf,'Name','Residual error PSD');

e1 = conv(ep, SAAF); e1=e1(1:K:end);
d1 = conv(dp, SAAF); d1=d1(1:K:end);
fprintf('Residual Noise Power: %f\n', mean(e1(end/2:end).^2));
fprintf('Noise Power: %f\n', mean(d1(end/2:end).^2));
fprintf('Atenuation (dB): %f\n', ...
    10*log10(mean(d1(end/2:end).^2)) - 10*log10(mean(e1(end/2:end).^2)));
saveas(gcf, '../results/fig5.png')
