% Anex to the Paper:
% Low Delay and Low Cost Sigma-Delta Adaptive Controller for Active Noise Control
% Paulo Lopes

global P M;

% simulation parameters
L = 3000;         % main simulation samples (at 44 kHz)
Li = 400;         % initial wait to fill the filter buffers
NW = 32;          % controler filter size
NS = 32;          % secondary path filter size
fs = 44100;       % lower sampling frequency
K = 32;           % oversample
P = 2;            % sigma delta order
M = 2^P;          % number of levels do the quantitizer

rng(7283723);

load AAF.mat
AAF = conv(AAF,AAF);
N_AAF = length(AAF);

w_delays = [3.7, 2.2, 5.3];
w_amplitudes = [1.1, -0.7, 0.5]/2;
wop = sinc((0:NW*K-1)'-w_delays*K)*w_amplitudes';

s_delays = [5.2, 7.2, 3.7];
s_amplitudes = [-1.3, 0.9, -0.5];
sp = sinc((0:NS*K-1)'-s_delays*K)*s_amplitudes';
sp1 = conv(sp, ones(K,1)/K);
sp1 = sp1(K/2+K/4:end); % with secondary path modeling error
%sp1 = sp1(K/2:end); % no secondary path modeling error
sh = conv(sp1, AAF);
sh = K*sh(1+(N_AAF-1)/2:K:end);
sh = sh(1:NS)';

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

NW0 = NW*K-K/2;
w0 = zeros(NW0,1);

NQNS = 4;
QNSx = cell(NQNS,1);
Levels = [2 1/K 4 10];
for i=1:NQNS
    QNSx{i} = QNS(Levels(i));
end
n1 = 1;

% logs
xqns = zeros(L*K, NQNS);
yqns = zeros(L*K, NQNS);
q = zeros(L*K, NQNS);

for n0=Li*K:L*K
    % discrete
    u0(n0) = step(QNSx{1}, up(n0));
    
    yq(n0) = - w0'*u0(n0:-1:n0-NW0+1);
    y0(n0) = step(QNSx{3}, yq(n0));

    % physical
    yp(n0) = wop'*up(n0:-1:n0-NW*K+1);
    dp(n0) = sp'*yp(n0:-1:n0-NW*K+1);
    ep(n0) = dp(n0) + sp'*y0(n0:-1:n0-NW*K+1);
    
    % discrete
    e0(n0) = step(QNSx{4}, ep(n0));
    
    if mod(n0, K)==0
        n = n0/K;
        u(n) = AAF*u0(n0:-1:n0-length(AAF)+1);
        e(n) = AAF*e0(n0:-1:n0-length(AAF)+1);
        d(n) = AAF*dp(n0:-1:n0-length(AAF)+1); % just for debug
        
        u1(n) = sh'*u(n:-1:n-NS+1);

        y(n) = w'*u(n:-1:n-NW+1);
        dh(n) = e(n) + sh'*y(n:-1:n-NS+1);
        eh(n) = dh(n) - w'*u1(n:-1:n-NW+1);

        u1v = u1(n:-1:n-NW+1);
        w = w + mu*eh(n)*u1v/(u1v'*u1v+1e-2);
    end
    
    if n1 > NW0
        n1 = 1;
        reset(QNSx{2});
    end
    w0(n1) = step(QNSx{2}, w(floor((n1-1)/K+0.5)+1)/K);
    n1 = n1 + 1;

    for i=1:NQNS
      xqns(n0,i) = QNSx{i}.x;
      yqns(n0,i) = QNSx{i}.yd;
      q(n0,i) = QNSx{i}.q;
    end
end

fprintf(1,'QNS: u w y e\n');
for i=1:NQNS
    m = sqrt(max(conv(xqns(:,i).^2,ones(10,1))/10));
    fprintf(1, 'QNS(%d) -- max input signal RMS: %f\n', i, m);
    
    x = conv(AAF,xqns(:,i)); x=x(length(AAF):end);
    y = conv(AAF,yqns(:,i)); y=y(length(AAF):end);
    normalized_mse = mean((y-x).^2)/(QNSx{i}.L)^2;
    fprintf(1,'QNS(%d) -- nmse: %f dB bits: %f\n', i, ...
        10*log10(normalized_mse),(-10*log10(3)-10*log10(normalized_mse))/(20*log10(2)));
end

figure(1);
ep_filtered = conv(ep, AAF);
ep_down_sampled = ep_filtered(K:K:end);
plot((0:length(e)-1)/fs*1e3-Li/fs*1e3, e); hold on;
plot((0:length(ep_down_sampled)-1)/fs*1e3-Li/fs*1e3, ep_down_sampled); hold off;
legend('e','e_p')
set(gca,'XLim',[0 (L-Li)/fs*1e3]);
xlabel('Time (ms)'); ylabel('Level');
saveas(gcf, '../results/fig1.png')

figure(2);
[hw0, f_hw0] = freqz(w0,1,1024*K,fs*K);
[hw1, f_hw1] = freqz(w,1,1024,fs);
hsh = freqz(ones(K,1)/K,1,1024*K,fs*K)./...
    freqz([zeros(1,K/2),1],1,1024*K,fs*K);

plot(f_hw0/1e3, 20*log10(abs(hw0./hsh))); hold on;
plot(f_hw1/1e3, 20*log10(abs(hw1))); hold off;
set(gca,'XLim', [0, fs/2/1e3]);
set(gca,'YLim', [-15, 5]);
set(gcf,'Name','Frequency Responce of W1 and W0');
grid on;
legend('w_0 filter', 'w filter','Location','southwest');
ylabel('Frequency Response Ampltitude (dB)');
xlabel('Frequency (kHz)');
set(gcf,'Name','Residual error PSD');
saveas(gcf, '../results/fig2.png')

figure(3);
[Pd, f_pd] = pwelch(dp(end/2:end), [], [], K*128, K*fs);
[Pe, f_pe] = pwelch(ep(end/2:end), [], [], K*128, K*fs);
plot(f_pd/1e3, 10*log10(Pd)); hold on;
plot(f_pe/1e3, 10*log10(Pe)); hold off;
set(gca, 'XLim', [0 fs/1e3]);
grid on;
legend('ANC off', 'ANC on','Location','southeast');
ylabel('Noise power spectral density (dB)');
xlabel('Frequency (kHz)');
set(gcf,'Name','Residual error PSD');
saveas(gcf, '../results/fig3.png')

fprintf('Residual Noise Power: %f\n', mean(e(end*3/4:end).^2));
fprintf('Noise Power: %f\n', mean(d(end*3/4:end).^2));
fprintf('Atenuation (dB): %f\n', ...
    10*log10(mean(d(end*3/4:end).^2)) - 10*log10(mean(e(end*3/4:end).^2)));