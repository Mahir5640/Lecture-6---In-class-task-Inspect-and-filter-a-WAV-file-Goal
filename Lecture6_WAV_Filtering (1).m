%% Lecture 6 - In-class task: Inspect and filter a WAV file
% Requires MATLAB and Signal Processing Toolbox.
% Select a WAV file containing a steady 50 Hz electrical hum.
clear; close all; clc;

[file, folder] = uigetfile('*.wav', 'Select the WAV recording');
if isequal(file,0)
    error('No WAV file selected. Run the script again and select a file.');
end
[x, Fs] = audioread(fullfile(folder,file));
[nSamples, nChannels] = size(x);
duration = nSamples/Fs;
fprintf('File: %s\nSample rate: %g Hz\nDuration: %.2f s\nChannels: %d\n', ...
    file, Fs, duration, nChannels);
if Fs <= 102
    error('Sample rate must exceed 102 Hz for the selected 49-51 Hz stop band.');
end

% Preserve the original channel count for filtering and playback.
% Use the channel average only for diagnostic plots.
xMono = mean(x,2);

%% Inspect the original recording BEFORE filtering
figure('Name','Original WAV inspection','Color','w');
subplot(2,1,1);
plot((0:nSamples-1)'/Fs,xMono,'Color',[0.2 0.3 0.6]);
xlim([0,min(0.5,duration)]); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Original waveform (first 0.5 s)');
subplot(2,1,2);
windowLength = min(1024,nSamples);
if windowLength < 16, error('Recording is too short for spectral analysis.'); end
noverlap = floor(0.75*windowLength);
nfft = max(2048,2^nextpow2(windowLength));
spectrogram(xMono,hamming(windowLength),noverlap,nfft,Fs,'yaxis');
title('Original spectrogram'); colorbar;

%% Design a narrow 50 Hz band-stop filter
% Only use this design if the inspection indicates a 50 Hz hum.
stopBand = [49 51]; % Hz
[b,a] = butter(2,stopBand/(Fs/2),'stop');
y = filtfilt(b,a,x); % Offline zero-phase filtering, channel by channel
filteredMono = mean(y,2);

%% Compare waveforms and spectrograms on consistent scales
segment = (0:min(nSamples-1,max(1,round(0.5*Fs))-1))';
commonLimit = max(abs([xMono(segment+1);filteredMono(segment+1)]));
if commonLimit == 0, commonLimit = 1; end
figure('Name','Original versus filtered WAV','Color','w','Position',[100 100 1100 750]);
subplot(2,2,1);
plot(segment/Fs,xMono(segment+1),'b');
xlim([0,segment(end)/Fs]); ylim([-commonLimit commonLimit]); grid on;
title('Original waveform'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,2,2);
plot(segment/Fs,filteredMono(segment+1),'r');
xlim([0,segment(end)/Fs]); ylim([-commonLimit commonLimit]); grid on;
title('Filtered waveform (49-51 Hz stop band)'); xlabel('Time (s)'); ylabel('Amplitude');
% Compute both spectrograms first, then use the same colour limits.
[S1,F,T] = spectrogram(xMono,hamming(windowLength),noverlap,nfft,Fs);
[S2,~,~] = spectrogram(filteredMono,hamming(windowLength),noverlap,nfft,Fs);
P1 = 20*log10(abs(S1)+eps); P2 = 20*log10(abs(S2)+eps);
cl = [max(-120,min([P1(:);P2(:)])),max([P1(:);P2(:)])];
if cl(1)>=cl(2), cl = cl(2)+[-1 0]; end
subplot(2,2,3);
imagesc(T,F,P1); axis xy; ylim([0,min(500,Fs/2)]); clim(cl); colorbar;
title('Original spectrogram'); xlabel('Time (s)'); ylabel('Frequency (Hz)');
subplot(2,2,4);
imagesc(T,F,P2); axis xy; ylim([0,min(500,Fs/2)]); clim(cl); colorbar;
title('Filtered spectrogram'); xlabel('Time (s)'); ylabel('Frequency (Hz)');
exportgraphics(gcf,'Lecture6_Original_vs_Filtered.png','Resolution',200);

%% Quantify attenuation from filter response (single pass vs filtfilt)
H50 = freqz(b,a,50,Fs);
attenuation50 = 40*log10(max(abs(H50),eps)); % filtfilt applies filter twice
fprintf('Designed stop band: %.0f-%.0f Hz\n',stopBand);
fprintf('Theoretical filtfilt attenuation at 50 Hz: %.2f dB\n',attenuation50);
fprintf('Comparison figure saved: Lecture6_Original_vs_Filtered.png\n');

%% Listen at identical gain (no automatic normalization)
fprintf('Playing original recording, then filtered recording...\n');
soundscDisabled = true; %#ok<NASGU> % Use sound(), not soundsc(), to avoid gain normalization.
sound(x,Fs); pause(duration+0.5);
sound(y,Fs);
fprintf('Listen for reduced hum and any change in useful audio.\n');
