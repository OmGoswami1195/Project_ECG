%basic syntax
clc;
clear all;
close all;

% Reading the ECG data from file
data = load("C:\Users\user1\Downloads\a01m.mat");% The load function is used to load data from .mat files
signal = data.val; % Extract the signal data from the loaded structure


% Signal pre-processing
processed_signal = signal - mean(signal);
% Centering the Data: Subtracting the mean of the signal from each data point effectively centers the signal around zero.
% Statistical Analysis: Centering the data is a common practice in statistical analysis. 
% centering the data can make it easier to interpret regression coefficients or other statistical parameters
% Noise Reduction: Subtracting the mean can also help remove constant or low-frequency noise from the signal. 
%Feature Extraction: Depending on the specific analysis being performed, centering the data might be a necessary step before extracting certain features or characteristics from the signal.



% Initializing parameters
fs = 200; % sampling frequency ,it means that 200 samples are taken from the continuous signal every second
% sampling frequency -It refers to the number of samples taken per unit of time from a continuous signal to convert it into a discrete signal suitable for digital processing and storage
%
t = (0:length(signal) - 1) / fs; % matlab can plot only the dicrete time signals , the selection of the t like this will aide to analyse the signal properly
% t represents the time axis before fourier transform
N = length(signal);

k = 0:N - 1;
F = k * fs / N; %  reprsents the frequency axis after fourier transform 

% FIR Notch Filter Design
f0 = 60; % Center frequency for the notch , or the cut off frequency

% Magnitude spectrum of the signal or performing fourier transform
processed_sig_fft = fft(processed_signal);

% Powerline noise removal using a notch filter or creation of the filter
w1 = 2 * pi * (f0 / fs);
% w1 = 2 * pi * (f0 / fs);
% This line calculates the digital radian frequency w1 corresponding to the desired notch frequency f0. 
% f0 represents the frequency that you want to attenuate using the notch filter, and fs is the sampling frequency of the signal. 
% Multiplying by 2 * pi converts the frequency from Hertz to radians.


% diffrence between poles and zeros
% Poles: Poles are the values of the complex variable at which the denominator of a transfer function becomes zero. 
% In simpler terms, they are the values that make the system's response tend to infinity. 
% Poles influence the stability and transient response of a system. In a filter, the positions of the poles affect how the filter attenuates or amplifies different frequencies. 
% Poles can also determine the decay rate of oscillations in a system.

% Zeros: Zeros are the values of the complex variable at which the numerator of a transfer function becomes zero. 
% Zeros correspond to frequencies where the system's response is attenuated. 
% In a filter, zeros indicate where specific frequencies are canceled or attenuated.
% They can be used to shape the frequency response of the filter.

z1 = [exp(1j * w1); exp(-1j * w1)];
% This creates a complex pair of zeros (z1) for the notch filter. 
% The zeros are placed at frequencies that correspond to f0 in the positive and negative directions on the complex plane. 
% The use of complex numbers helps define the frequency and phase characteristics of the filter

[b_notch, a_notch] = zp2tf(z1, [0; 0], 1 / sum(poly(z1)));
% This line converts the zeros z1 and two poles (at the origin) into filter coefficients for a discrete-time transfer function using the zp2tf function. 
% b_notch represents the feedforward coefficients (numerator coefficients) of the filter, and a_notch represents the feedback coefficients (denominator coefficients).

% The specific filter coefficients are calculated such that the notch filter has a gain of 1 at DC (frequency 0 Hz). 
% The sum(poly(z1)) calculates the polynomial coefficients for the zeros, and 1 / sum(poly(z1)) scales the filter to achieve the desired DC gain.

filtered_signal = filtfilt(b_notch, a_notch, processed_signal);
% This line applies the notch filter defined by the b_notch and a_notch coefficients to the processed_signal. 
% The filtfilt function performs zero-phase filtering, which filters the signal forward and then backward to eliminate phase distortion.
% Phase distortion refers to an alteration of the phase relationships between different frequency components of a signal, leading to a distortion of the signal's temporal characteristics. 
% In other words, phase distortion occurs when the relative timing or alignment of different frequencies in a signal is changed, often due to the nonlinear behavior of a system or component.
% The result is stored in filtered_signal, which now represents the processed signal with the specified notch frequency attenuated.

% Magnitude spectrum of the filter
[h, f] = freqz(b_notch, a_notch, 1024, fs);
% freqz(b_notch, a_notch, 1024, fs) - This function call calculates the frequency response of the filter defined by the coefficients b_notch (feedforward) and a_notch (feedback). Let's break down the parameters:
% b_notch and a_notch: These are the feedforward and feedback coefficients of the filter that define its transfer function.
% 1024: This is the number of points at which the frequency response will be evaluated. More points provide a finer-grained view of the frequency response.
% fs: The sampling frequency of the signal. This parameter is used to correctly scale the frequency values in the output.
% The output of freqz provides the complex frequency response of the filter. The variable h will contain this complex response, which includes both magnitude and phase information for different frequencies. 
% The variable f will contain the corresponding frequency values in Hertz.


% Filtered ECG signal
filtered_ecg = filtfilt(b_notch, a_notch, filtered_signal);

processed_filtered_signal = filtered_ecg - mean(filtered_ecg);

filtered_signal_fft = fft(processed_filtered_signal);

% Plots
% Question 1
figure;
subplot(2, 1, 1);
plot(t, processed_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original ECG Signal in Time Domain');

subplot(2, 1, 2);
plot(F, abs(processed_sig_fft));
title('Magnitude Spectrum of Original ECG Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Question 2
figure;
plot(f, abs(h));
title('Magnitude Spectrum of Notch Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Question 3
figure;
zplane(b_notch, a_notch);
title('Pole-Zero Plot of Notch FIR Filter');

% Question 4
figure;
subplot(2, 1, 1);
plot(t, processed_filtered_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered ECG Signal in Time Domain');

subplot(2, 1, 2);
plot(F, abs(filtered_signal_fft));
title('Magnitude Spectrum of Filtered ECG Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');