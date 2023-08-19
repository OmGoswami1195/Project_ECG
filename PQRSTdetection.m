close all;clear;clc;
sig=load("C:\Users\user1\Downloads\Telegram Desktop\ecg3.dat");
N=length(sig);
fs=200;
t=[0:N-1]/fs;
figure(1);subplot(4,2,1);plot(sig)
title('Original Signal')

%%
     %           Low Pass Filter

b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];
sigL=filter(b,a,sig);
subplot(4,2,3);plot(sigL)
title('Low Pass Filter')
subplot(4,2,4);zplane(b,a)

%%
     %           High Pass Filter

b=[-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32];
a=[1 -1];
sigH=filter(b,a,sigL);
subplot(4,2,5);plot(sigH)
title('High Pass Filter')
subplot(4,2,6);zplane(b,a)

%%
     %          Derivative Base Filter

b=[1/4 1/8 0 -1/8 -1/4];
a=[1];
sigD=filter(b,a,sigH);
subplot(4,2,7);plot(sigD)
title('Derivative Base Filter')
subplot(4,2,8);zplane(b,a)

%%
     %      be tavane 2 miresanim
sigD2=sigD.^2;

%%
     %      normalization
signorm=sigD2/max(abs(sigD2));

%%
  
h=ones(1,31)/31;
sigAV=conv(signorm,h);
sigAV=sigAV(15+[1:N]);
sigAV=sigAV/max(abs(sigAV));
figure(2);plot(sigAV)
title('Moving Average filter')

%%
treshold=mean(sigAV);
P_G= (sigAV>0.01);
figure(3);plot(P_G)
title('treshold Signal')
figure;plot(sigL)
%%
difsig=diff(P_G);
left=find(difsig==1);
raight=find(difsig==-1);

%%
     %      run cancel delay
     %      6 sample delay because of LowPass filtering
     %      16 sample delay because of HighPass filtering
left=left-(6+16);
raight=raight-(6+16);

%%
    % P-QRS-t
for i=1:length(left);
   
    [R_A(i) R_t(i)]=max(sigL(left(i):raight(i)));
    R_t(i)=R_t(i)-1+left(i) %add offset
   
    [Q_A(i) Q_t(i)]=min(sigL(left(i):R_t(i)));
    Q_t(i)=Q_t(i)-1+left(i)
  
    [S_A(i) S_t(i)]=min(sigL(left(i):raight(i)));
    S_t(i)=S_t(i)-1+left(i)
    
    [P_A(i) P_t(i)]=max(sigL(left(i):Q_t(i)));
    P_t(i)=P_t(i)-1+left(i)
    
    [T_A(i) T_t(i)]=max(sigL(S_t(i):raight(i)));
    T_t(i)=T_t(i)-1+left(i)+47
    
   
end

%%



figure;plot(t,sigL,t(Q_t),Q_A,'*g',t(S_t),S_A,'^k',t(R_t),R_A,'ob',t(P_t),P_A,'+b',t(T_t),T_A,'+r');
for i=1:((length(P_t))-1)
    
    HRV=P_t(i+1)-P_t(i)
end
% Hyperkalaemia Tented T waves detection
hyperkalaemia_detected = false(length(P_t), 1);
for i = 1:length(P_t)
    T_window_start = P_t(i) + round(0.15 * fs); % Window starts 0.15 seconds after P wave
    T_window_end = P_t(i) + round(0.4 * fs);    % Window ends 0.4 seconds after P wave
    
    % Extract a segment of the ECG signal around the T wave
    T_wave_segment = sigL(T_window_start:T_window_end);
    
    % Calculate T wave amplitude
    T_amplitude = max(T_wave_segment) - min(T_wave_segment);
    
    % Determine if T wave meets Hyperkalaemia criteria (tenting)
    if T_amplitude >= 0.5 % Adjust threshold as needed
        hyperkalaemia_detected(i) = true;
    end
end

% Display detected Hyperkalaemia Tented T waves
figure;
plot(t, sigL, t(Q_t), Q_A, '*g', t(S_t), S_A, '^k', t(R_t), R_A, 'ob', t(P_t), P_A, '+b', t(T_t), T_A, '+r');
hold on;
plot(t(P_t(hyperkalaemia_detected)), P_A(hyperkalaemia_detected), 'or', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG with Detected Hyperkalaemia Tented T Waves');
legend('ECG', 'Q Waves', 'S Waves', 'R Waves', 'P Waves', 'T Waves', 'Detected Hyperkalaemia Tented T Waves');
grid on;
% Calculate Heart Rate (HR) and Heart Rate Variability (HRV)
RR_intervals = diff(P_t); % Calculate RR intervals (in samples)
RR_intervals_ms = RR_intervals / fs * 1000; % Convert RR intervals to milliseconds
mean_RR = mean(RR_intervals_ms); % Calculate mean RR interval
std_RR = std(RR_intervals_ms); % Calculate standard deviation of RR intervals

% Calculate heart rate (bpm)
heart_rate_bpm = 60 / (mean_RR / 1000);

% Determine if heart rate is normal or abnormal
if heart_rate_bpm >= 60 && heart_rate_bpm <= 100 % Normal heart rate range (60-100 bpm)
    heart_rate_status = 'Normal';
else
    heart_rate_status = 'Abnormal';
end

% Display heart rate and heart rate variability
fprintf('Mean RR interval: %.2f ms\n', mean_RR);
fprintf('Standard deviation of RR intervals: %.2f ms\n', std_RR);
fprintf('Heart rate: %.2f bpm\n', heart_rate_bpm);
fprintf('Heart rate status: %s\n', heart_rate_status);


    