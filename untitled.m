clc;
clear all;
close all;

% Reading the ECG data from file
data = load("C:\Users\user1\OneDrive\Desktop\ECG signals (1000 fragments)\MLII\1 NSR\100m (2).mat");
signal = data.val; % Extract the signal data from the loaded structure

fs = 360;
filterorder = 10;
cutofflow = 18;
cutoffhigh = 4;

% Design Low-pass Butterworth filter
[b, a] = butter(filterorder, cutofflow/(fs/2), 'low');
[h, w] = freqz(b, a);

% Apply Low-pass filter
eeglowpass = filter(b, a, signal);

% Design High-pass Butterworth filter
[b1, a1] = butter(filterorder, cutoffhigh/(fs/2), 'high');
[h1, w1] = freqz(b1, a1);

% Apply High-pass filter
eeghighpass = filter(b1, a1, eeglowpass);

% Find R-wave peaks
[pks, locs] = findpeaks(eeghighpass(1801:end));
threshold = 2 * (rms(eeghighpass(1801:end)));
rwave = pks > threshold;
rwaveform = zeros(size(eeghighpass(1801:end)));
rwaveform(locs(rwave)) = pks(rwave);

% Calculate Heart Rate
beat = sum(rwave);
hr = (beat / 5) * 60; % Assuming a 5-second window

% Calculate RR intervals
rr_intervals = diff(locs) / fs;

% Calculate SDNN (Standard Deviation of RR intervals)
sdnn = std(rr_intervals);

% Plotting
figure;

subplot(3, 1, 1);
plot(eeglowpass);
title('Low-pass Filtered ECG Signal');
xlabel('Sample');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(eeghighpass);
title('High-pass Filtered ECG Signal');
xlabel('Sample');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(eeghighpass(1801:end));
hold on;
plot(rwaveform);
title('R-wave Detection');
xlabel('Sample');
ylabel('Amplitude');
legend('Filtered ECG', 'Detected R-wave');

% Display Results
fprintf("Heart Rate: %.2f bpm\n", hr);
fprintf("Heart Rate Variability (SDNN): %.2f ms\n", sdnn * 1000);
