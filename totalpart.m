%% DSP part B project
%% professor Doctor rajabi
%% student Mastooreh sadeghi
%% Part 2 - extracting first 1 mins data and Plotting it
clear all
clc
load('100m.mat')
fs = 360;    % sampling frequency
[tm,signal,Fs,siginfo]=rdmat('100m');   % 'wfdb-app-toolbox-0-9-9' Add to path
given_data_in_mins = (length(signal)/fs)/60;
given_data_in_secs = round(length(signal)/(given_data_in_mins*60));
S = (given_data_in_secs * 60);
data= signal(1:S,:);               % first 1 mins data
figure(1);
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('1 mins ECG data'); % plotting the data                                               % small section of the data

%% part 3 - get fft and fftshift from ECG signal and plot data in frequency domain
NFFT = 2 ^ nextpow2(length(data));  %compute FFT length depends on the signal length
Y = fft(data,NFFT);  %compute the fft
Y2 = fftshift(Y);
Y = Y(1:NFFT/2);  %we only need a one sided fft plot
Y_abs = 1/NFFT*abs(Y); %calculate the magnitude and normalize the spectrum
f_fft = (0:NFFT/2-1)*fs/NFFT; %scale the frequency axe and calculate the coresponding frequencys
figure(2)
plot(f_fft,Y_abs);
title('Amplitude Spectrum of ECG Signal');
xlabel('f (Hz)');
ylabel('Magnitude |P1(f)|');

%% part4 - use filter to remove noise method1

Pass_edge_freq =50; %pass edge freq
Trans_band =10; %transition band
N= round((3.3*Fs)/Trans_band);%N according to hamming window 
Mid_point= floor(N/2); % finding the midpoint
fc=(Pass_edge_freq+(Trans_band/Fs))/Fs; 
Hd(1)= 2*fc;
Wd(1)= 1; %according to Hamming window
H(1) = Hd(1) * Wd(1);
for x=2:N
    
Hd(x)=(2*fc*sin(2*pi*fc*(x-1)))/(2*fc*pi*(x-1)); 
Wd(x)=0.54 + (0.46*cos(((x-1)*2*pi)/N));

H(x)=Hd(x)*Wd(x);

end

figure(3)
plot(H); title('filter')
New_data=filter(H',1,data); %filtered data
figure(4)
plot(data'); hold on; plot(New_data)
legend('Input Data','Filtered Data');xlabel('samples'); ylabel('Amplitude(mv)');
figure(5)
subplot(2,1,1)
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('1 mins ECG data'); % plotting the data
subplot(2,1,2) 
plot(New_data); xlabel('samples'); ylabel('Amplitude(mv)'); title('1 mins filtered ECG data');


%% part4 - use filter to remove noise method2

[b,a] = butter(3,[0.0196 0.0204],'stop'); % Bandstop = 49Hz-51Hz, 6th Order
z = filter(b,a,data);

% butterworth filter
fc = 80;    % cuttoff frequency
fs = 5000;  % sampling frequency
N = 3;      % 3rd order filter
Wn = fc/(fs/2); % freq in rad/sample
[b,a] = butter(N, Wn);
New_data=filter(b,a,data); %filtered data
figure(6)
plot(data'); hold on; plot(New_data)
legend('Input Data','Filtered Data');xlabel('samples'); ylabel('Amplitude(mv)');

% frequency response
[H,W]= freqz(b, a);
magH = abs(H);
angH = angle(H);
% plot mangnitude and phase response
figure(7);
subplot(2,1,1); plot(W/pi, magH); grid
title('Magnitude response');
ylabel('Magnitude');
subplot(2,1,2); plot(W/pi, angH); grid
title('Phase response');
xlabel('Frequency in pi units');
ylabel('Phase in pi Radians');
figure(8)
subplot(2,1,1)
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('1 mins ECG data'); % plotting the data
subplot(2,1,2) 
plot(New_data); xlabel('samples'); ylabel('Amplitude(mv)'); title('1 mins filtered ECG data');


%% Part 5 - baseline correction the signal defention1
dataneww = detrend(data,1);
figure(9);
subplot(2,1,1)
plot(dataneww); xlabel('samples'); ylabel('Amplitude(mv)'); title('after  baseline correction'); % plotting the data                                               % small section of the data
subplot(2,1,2)
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('before  baseline correction'); % plotting the data                                               % small section of the data



%% Part 5 - baseline correction the signal defention2 method1
Ts = 1/Fs;
Fn = Fs/2;  
Wp = [2.1  15.0]/Fn;                                                % Passband Frequency (Normalised)
Ws = [1.8  18.0]/Fn;                                                % Stopband Frequency (Normalised)
Rp =   1;                                                           % Passband Ripple (dB)
Rs =  50;                                                           % Stopband Ripple (dB)
[n,Ws]  = cheb2ord(Wp,Ws,Rp,Rs);                                    % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                                          % Filter Design, Sepcify Bandpass
[sos,g] = zp2sos(z,p,k);                                            % Convert To Second-Order-Section For Stability
figure(10)
freqz(sos, 2^16, Fs)                                                % Filter Bode Plot
D_Filtered = filtfilt(sos, g, data);                                   % Filter Signal
figure(11)
subplot(2,1,1);
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('BEFORE  baseline correction'); % plotting the data
subplot(2,1,2);
plot(D_Filtered);  xlabel('samples'); ylabel('Amplitude(mv)'); title('AFTER  baseline correction '); % plotting the data
grid
%% Part 5 - baseline correction the signal   method2  
y = data
figure(12)
subplot(4,1,1);
plot(y);
for i=1:3600
y(i)=y(i)+200*sin(2*pi*.5*i/360)+20*sin(2*pi*50*i/360);
end
subplot(4,1,2);
plot(y);
fc=40;
fs=360;
[b,a] = butter(6,fc/(fs/2));
x=filter(b,a,y);
y1=y(21:end);
x1=x(21:end);
subplot(4,1,3);
plot(x1);
[b1,a1] = butter(3,[0.00005 0.01],'stop');
x2=filter(b1,a1,x1);
y3=y1(201:end);
x3=x2(201:end);
% frequency response
[H,W]= freqz(b, a);
magH = abs(H);
angH = angle(H);
% plot mangnitude and phase response
figure(13);
subplot(2,1,1); plot(W/pi, magH); grid
title('Magnitude response');
ylabel('Magnitude');
subplot(2,1,2); plot(W/pi, angH); grid
title('Phase response');
xlabel('Frequency in pi units');
ylabel('Phase in pi Radians');






