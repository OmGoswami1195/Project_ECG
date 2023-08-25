clc;                % Clear the command window
clear all;          % Clear all workspace variables
close all;          % Close all open figures/windows

delay = 0;          % Initialize a variable 'delay' with a value of 0
skip = 0;           % Initialize a variable 'skip' with a value of 0
% 'skip' is used to indicate whether a T wave has been detected

m_selected_RR = 0;  % Initialize a variable 'm_selected_RR' with a value of 0
% This variable seems to be related to RR interval analysis

mean_RR = 0;        % Initialize a variable 'mean_RR' with a value of 0
% This variable might be used to store the mean RR interval

ser_back = 0;       % Initialize a variable 'ser_back' with a value of 0
% It's unclear what 'ser_back' stands for without further context

ax = zeros(1,6);    % Initialize an array 'ax' with 6 zeros
% 'ax' seems to be used for storing plot axes information

fs = 360;           % Assign the value 200 to the variable 'fs'
% 'fs' likely represents the sampling frequency of the ECG signal

% Load ECG data from the specified file path
ecg_data = load("C:\Users\user1\OneDrive\Desktop\ECG signals (1000 fragments)\MLII\1 NSR\100m (0).mat");

% Extract the 'val' field from the loaded data and convert it to double
ecg = double(ecg_data.val);

% Remove the mean from the ECG signal
ecg = ecg - mean(ecg);

%It has come to my attention the original filter doesn't achieve 12 Hz
% so some modifications are made 

% Calculate the normalized cutoff frequency for the low-pass filter
Wn = 12*2/fs;

% Wn = 12*2/fs;: This calculates the normalized cutoff frequency (Wn) for the low-pass    filter. The value 12 is the desired cutoff frequency in Hz, and fs is the sampling %frequency of 
% the ECG signal. The multiplication by 2 is a common practice to convert the cutoff frequency from Hz to the normalized frequency range of [0, 1] that the Butterworth filter function
% expects




% Choose the filter order (higher order leads to sharper roll-off)
N = 3;

% Design a Butterworth low-pass filter
[a,b] = butter(N, Wn, 'low');

% Apply zero-phase forward and reverse filtering to the ECG signal
ecg_l = filtfilt(a, b, ecg);

% This uses the filtfilt function to apply the designed filter coefficients a and b to the ECG signal ecg. The filtfilt function performs zero-phase forward and reverse filtering, 
% which   reduces phase distortion.

% Normalize the filtered ECG signal
ecg_l = ecg_l / max(abs(ecg_l));

% After filtering, the filtered ECG signal is normalized by dividing it by the maximum absolute value of the filtered signal. This step scales the signal to have values in the 
% range [-1, 1].

figure ;
plot(ecg);
axis tight;
xlabel("samples");
ylabel("amplitude");
title('Raw signal');
figure ;
plot(ecg_l);
xlabel("samples");
ylabel("amplitude");
axis tight;
%The command axis tight; is a MATLAB command used to adjust the axes limits of a plot or figure to tightly fit the data being plotted.
title('Low pass filtered');

% original code  doesnt achieve  achieve 5 Hz frequecy cut off
% high pass filter
Wn = 5*2/fs;
N = 3;                                                                  % order of 3 less processing
[a,b] = butter(N,Wn,'high');                                            % bandpass filtering
ecg_h = filtfilt(a,b,ecg_l); 
ecg_h = ecg_h/ max(abs(ecg_h));
figure;
plot(ecg_h);axis tight;xlabel("samples");
ylabel("amplitude");title('High Pass Filtered');
% if fs is not equal to 200 then direct implentation of bandpass filter
f1=5;                                                                      % cuttoff low frequency to get rid of baseline wander
 f2=15;                                                                     % cuttoff frequency to discard high frequency noise
 Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
 N = 3;                                                                     % order of 3 less processing
 [a,b] = butter(N,Wn);                                                      % bandpass filtering
 ecg_h = filtfilt(a,b,ecg);
 ecg_h = ecg_h/ max( abs(ecg_h));
 figure
 %plot(ecg);axis tight;title('Raw Signal');
 plot(ecg_h);axis tight;
xlabel("samples");
ylabel("amplitude");
title('Raw signal');
 title('Band Pass Filtered');



%% ==================== derivative filter ========================== %%
% ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- %
if fs ~= 200
    % Calculate interpolation coefficient for resampling
 int_c = (5-1)/(fs*1/40);
  % Generate the derivative filter coefficients using interpolation
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
     % Use the default derivative filter coefficients for fs = 200
 b = [1 2 0 -2 -1].*(1/8)*fs;   
end
 ecg_d = filtfilt(b,1,ecg_h);
 ecg_d = ecg_d/max(ecg_d);
 figure;
 plot(ecg_d);
xlabel("samples");
ylabel("amplitude");
title('derivative');
 axis tight;
 % Interpolation for Resampling: If the sampling frequency is not 200, the code calculates an interpolation coefficient int_c based on the desired resampling rate (fs * 1/40). This 
% coefficient is used to interpolate the derivative filter coefficients.

%% ========== Squaring nonlinearly enhance the dominant peaks ========== %%
 ecg_s = ecg_d.^2;
 figure
 plot(ecg_s);
  axis tight;
  title('Squared');
  xlabel("samples");
ylabel("amplitude");

%% ============  Moving average ================== %%
%-------Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]---------%

% Define the length of the moving average window (in seconds)
window_length = 0.150;  % 150 milliseconds

% Calculate the number of samples in the moving average window
window_size = round(window_length * fs);

% Apply the moving average filter using convolution
ecg_m = conv(ecg_s, ones(1, window_size) / window_size);

% Adjust the delay variable to account for the filter's delay
delay = delay + window_size / 2;
% Applying Moving Average Filter: The conv function is used to apply the moving average filter. It convolves the ECG signal ecg_s with a filter kernel consisting of ones divided by the
% window size. This averaging operation effectively calculates the mean of the signal within the window.
% Applying Moving Average Filter: The conv function is used to apply the moving average filter. It convolves the ECG signal ecg_s with a filter kernel consisting of ones divided by the window size. This averaging operation effectively calculates the mean of the signal within the window.
figure
plot(ecg_m);
 xlabel("samples");
ylabel("amplitude");
  axis tight;
  title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
  axis tight;

  %% ===================== Fiducial Marks ============================== %% 
% Note : a minimum distance of 40 samples is considered between each R wave
% since in physiological point of view no RR wave can occur in less than
% 200 msec distance
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs));
%Peak Detection using findpeaks: The findpeaks function is employed to detect peaks (QRS complexes) in the smoothed ECG signal ecg_m. The 'MINPEAKDISTANCE' option enforces a minimum distance of 200 ms (40 samples at 200 Hz) between consecutive R waves. The detected peak amplitudes are stored in pks, and their locations in samples are stored in locs.

%% =================== Initialize Some Other Parameters =============== %%
LLp = length(pks);
% ---------------- Stores QRS wrt Sig and Filtered Sig ------------------%
qrs_c = zeros(1,LLp);           % amplitude of R
qrs_i = zeros(1,LLp);           % index
qrs_i_raw = zeros(1,LLp);       % amplitude of R
qrs_amp_raw= zeros(1,LLp);      % Index
% ------------------- Noise Buffers ---------------------------------%
nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);
% ------------------- Buffers for Signal and Noise ----------------- %
SIGL_buf = zeros(1,LLp);
NOISL_buf = zeros(1,LLp);
SIGL_buf1 = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf1 = zeros(1,LLp);
THRS_buf = zeros(1,LLp);
% Initialize Parameters: Several arrays are initialized to store information about QRS complexes and related thresholds and levels. LLp is calculated as the length of the detected peaks pks.
%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(ecg_m(1:2*fs))*1/3;                                          % 0.25 of the max amplitude 
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;                                       % 0.5 of the mean signal is considered to be noise
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;
%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
%Threshold Initialization: Initial thresholds and signal/noise levels are determined based on the maximum and mean amplitudes of the first 2 seconds of the smoothed ECG signal.
%Initialize Bandpass Filter Threshold: Similarly, initial thresholds and levels are determined for the bandpass-filtered signal.
THR_SIG1 = max(ecg_h(1:2*fs))*1/3;                                          % 0.25 of the max amplitude 
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2; 
SIG_LEV1 = THR_SIG1;                                                        % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1;                                                    % Noise level in Bandpassed filter
%% ============ Thresholding and desicion rule ============= %%
Beat_C = 0;                                                                 % Raw Beats
Beat_C1 = 0;                                                                % Filtered Beats
Noise_Count = 0;  % Noise Counter
%QRS Detection Loop: This loop iterates through each detected peak (pks) and its corresponding location (locs).

%Locate Corresponding Peak in Filtered Signal: The code checks if the corresponding R wave peak exists in the filtered signal ecg_h. It finds the maximum amplitude within a specified window around the peak's location.
for i = 1 : LLp  
   %% ===== locate the corresponding peak in the filtered signal === %%
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)
          [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
       else
          if i == 1
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
          elseif locs(i)>= length(ecg_h)
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
          end       
    end       
  %% ================= update the heart_rate ==================== %% 
    if Beat_C >= 9        
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));                                   % calculate RR interval
        mean_RR = mean(diffRR);                                            % calculate the mean of 8 previous R waves interval
        comp =qrs_i(Beat_C)-qrs_i(Beat_C-1);                                     % latest RR
    
        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
     % ------ lower down thresholds to detect better in MVI -------- %
                THR_SIG = 0.5*(THR_SIG);
                THR_SIG1 = 0.5*(THR_SIG1);               
        else
            m_selected_RR = mean_RR;                                       % The latest regular beats mean
        end 
          
    end
    
 %% Calculate Mean of Last 8 RR Waves: The mean of the last 8 RR intervals is calculated to ensure that the current QRS complex is not missed due to irregularities in heart rate.== calculate the mean last 8 R waves to ensure that QRS is not ==== %%
       if m_selected_RR
           test_m = m_selected_RR;                                         %if the regular RR availabe use it   
       elseif mean_RR && m_selected_RR == 0
           test_m = mean_RR;   
       else
           test_m = 0;
       end
        
    if test_m
          if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)                  % it shows a QRS is missed 
              [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs))); % search back and locate the max in this interval
              locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;      % location 
             
              if pks_temp > THR_NOISE
               Beat_C = Beat_C + 1;
               qrs_c(Beat_C) = pks_temp;
               qrs_i(Beat_C) = locs_temp;      
              % ------------- Locate in Filtered Sig ------------- %
               if locs_temp <= length(ecg_h)
                  [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):locs_temp));
               else
                  [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):end));
               end
              % ----------- Band pass Sig Threshold ------------------%
               if y_i_t > THR_NOISE1 
                  Beat_C1 = Beat_C1 + 1;
                  qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs)+ (x_i_t - 1);% save index of bandpass 
                  qrs_amp_raw(Beat_C1) = y_i_t;                               % save amplitude of bandpass 
                  SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;                      % when found with the second thres 
               end
               
               not_nois = 1;
               SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;                       % when found with the second threshold             
             end             
          else
              not_nois = 0;         
          end
    end
  
    %% ===================  find noise and QRS peaks ================== %%
    if pks(i) >= THR_SIG      
      % ------ if No QRS in 360ms of the previous QRS See if T wave ------%
       if Beat_C >= 3
          if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
              Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));       % mean slope of the waveform at that position
              Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C)))); % mean slope of previous R wave
              if abs(Slope1) <= abs(0.5*(Slope2))                              % slope less then 0.5 of previous R
                 Noise_Count = Noise_Count + 1;
                 nois_c(Noise_Count) = pks(i);
                 nois_i(Noise_Count) = locs(i);
                 skip = 1;                                                 % T wave identification
                 % ----- adjust noise levels ------ %
                 NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                 NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
              else
                 skip = 0;
              end
            
           end
        end
        %---------- skip is 1 when a T wave is detected -------------- %
        if skip == 0    
          Beat_C = Beat_C + 1;
          qrs_c(Beat_C) = pks(i);
          qrs_i(Beat_C) = locs(i);
        
        %--------------- bandpass filter check threshold --------------- %
          if y_i >= THR_SIG1  
              Beat_C1 = Beat_C1 + 1;
              if ser_back 
                 qrs_i_raw(Beat_C1) = x_i;                                 % save index of bandpass 
              else
                 qrs_i_raw(Beat_C1)= locs(i)-round(0.150*fs)+ (x_i - 1);   % save index of bandpass 
              end
              qrs_amp_raw(Beat_C1) =  y_i;                                 % save amplitude of bandpass 
              SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;                       % adjust threshold for bandpass filtered sig
          end
         SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;                          % adjust Signal level
        end
              
    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                        % adjust Noise level in filtered sig
         NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                       % adjust Noise level in MVI       
    elseif pks(i) < THR_NOISE
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);    
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                         % noise level in filtered signal    
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                        % adjust Noise level in MVI     
    end
               
    %% ================== adjust the threshold with SNR ============= %%
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    %------ adjust the threshold with SNR for bandpassed signal -------- %
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    
%--------- take a track of thresholds of smoothed signal -------------%
SIGL_buf(i) = SIG_LEV;
NOISL_buf(i) = NOISE_LEV;
THRS_buf(i) = THR_SIG;
%-------- take a track of thresholds of filtered signal ----------- %
SIGL_buf1(i) = SIG_LEV1;
NOISL_buf1(i) = NOISE_LEV1;
THRS_buf1(i) = THR_SIG1;
% ----------------------- reset parameters -------------------------- % 
skip = 0;                                                   
not_nois = 0; 
ser_back = 0;    
end
%% ======================= Adjust Lengths ============================ %%
qrs_i_raw = qrs_i_raw(1:Beat_C1);
qrs_amp_raw = qrs_amp_raw(1:Beat_C1);
qrs_c = qrs_c(1:Beat_C);
qrs_i = qrs_i(1:Beat_C);
figure
hold on,scatter(qrs_i,qrs_c,'m');
xlabel("samples");
ylabel("amplitude");
  hold on,plot(locs,NOISL_buf,'--k','LineWidth',2);
  hold on,plot(locs,SIGL_buf,'--r','LineWidth',2);
  hold on,plot(locs,THRS_buf,'--g','LineWidth',2);
 if any(ax)
  ax(~ax) = []; 
  linkaxes(ax,'x');
  zoom on;
 end
figure
plot(ecg_h);xlabel("samples");
ylabel("amplitude");
   title('QRS on Filtered Signal');
   axis tight;
   hold on,scatter(qrs_i_raw,qrs_amp_raw,'m');
   hold on,plot(locs,NOISL_buf1,'LineWidth',2,'Linestyle','--','color','k');
   hold on,plot(locs,SIGL_buf1,'LineWidth',2,'Linestyle','-.','color','r');
   hold on,plot(locs,THRS_buf1,'LineWidth',2,'Linestyle','-.','color','g');
 figure;
 plot(ecg_m);xlabel("samples");
ylabel("amplitude");
   title('QRS on MVI signal and Noise level(black),Signal Level (red) and Adaptive Threshold(green)');axis tight;
   hold on,scatter(qrs_i,qrs_c,'m');
   hold on,plot(locs,NOISL_buf,'LineWidth',2,'Linestyle','--','color','k');
   hold on,plot(locs,SIGL_buf,'LineWidth',2,'Linestyle','-.','color','r');
   hold on,plot(locs,THRS_buf,'LineWidth',2,'Linestyle','-.','color','g');
    figure;xlabel("samples");
ylabel("amplitude");
   plot(ecg-mean(ecg));
   xlabel("samples");
ylabel("amplitude");
   title('Pulse train of the found QRS on ECG signal');
   axis tight;
   line(repmat(qrs_i_raw,[2 1]),...
       repmat([min(ecg-mean(ecg))/2; max(ecg-mean(ecg))/2],size(qrs_i_raw)),...
       'LineWidth',2.5,'LineStyle','-.','Color','r');

RR_intervals = diff(qrs_i) / fs;  % Calculate RR intervals in seconds
mean_RR = mean(RR_intervals);      % Calculate mean RR interval

heart_rate = 60 / mean_RR;         % Calculate heart rate in beats per minute

fprintf('Calculated Heart Rate: %.2f BPM\n', heart_rate);
RR_intervals = diff(qrs_i) / fs;  % Calculate RR intervals in seconds
mean_RR = mean(RR_intervals);      % Calculate mean RR interval

heart_rate = 60 / mean_RR;         % Calculate heart rate in beats per minute

fprintf('Calculated Heart Rate: %.2f BPM\n', heart_rate);

RR_intervals = diff(qrs_i) / fs;

% Calculate SDNN (Standard Deviation of RR intervals)
SDNN = std(mean_RR);

% Calculate RMSSD (Root Mean Square of Successive RR Interval Differences)
RMSSD = sqrt(mean(diff(RR_intervals).^2));

% Calculate SDSD (Standard Deviation of Successive RR Interval Differences)
SDSD = std(diff(RR_intervals));

% Calculate NNSD (Mean of Successive RR Interval Differences)
NNSD = mean(abs(diff(RR_intervals)));

fprintf('SDNN: %.2f ms\n', SDNN);
fprintf('RMSSD: %.2f ms\n', RMSSD);
fprintf('SDSD: %.2f ms\n', SDSD);
fprintf('NNSD: %.2f ms\n', NNSD);

