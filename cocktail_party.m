clear, clf, format compact;
clear,format compact;
% Makes a matrix out of wav files in project directory.
% The arreys will all be in the same length as the shortest one.
% Files longer than the shortest one will be truncated.
#disp('Reading .wav files from project directory');
[mic_1,fs1,nb1]=wavread('rightmic.wav'); %Reading file from right microphone.
[mic_2,fs2,nb2]=wavread('leftmic.wav'); %Reading file from left microphone.
size(mic_1)
size(mic_2)
% The below operation makes them 1xN:
fs1
fs2
mic_1=mic_1';
mic_2=mic_2';
if length(mic_1)>length(mic_2)
	mic_1=mic_1(1:length(mic_2));
else
	mic_2=mic_2(1:length(mic_1));
end
size(mic_1)
size(mic_2)
#disp('Playing the recording of the right microphone (closer to music)');
%soundsc(mic_1)
#disp('Playing the recording of the left microphone (closer to me)');
%soundsc(mic_2)

subplot(2,1,1)
plot(mic_1), axis([0 16000 -0.12 0.12]);
title('Right microphone (closer to music)')
xlabel('Sampled points');
subplot(2,1,2)
plot(mic_2), axis([0 16000 -0.12 0.12]);
title('Left microphone (closer to me)')
xlabel('Sampled points');

% I also choose to look at the frequency spectra of the signal:
Fs=fs1; % Sampling frequency
Fs
% Matrix with the frequency spectra from the two microphones:
fftsounds=[real(fft(mic_1,Fs));real(fft(mic_2,Fs))];
f=[1:Fs/2];

figure(2)
subplot(2,1,1)
plot(f,fftsounds(1,f)), axis([0 4000 -15 15]);
title('Frequency spectra of the right microphone')
xlabel('Frequency (Hz)');
subplot(2,1,2)
plot(f,fftsounds(2,f)), axis([0 4000 -15 15]);
title('Frequency spectra of the left microphone')
xlabel('Frequency (Hz)');

% At first I tried the algorithm in the time domain, and it didn't
% manage to separate the two sources very well.
% After that I used the frequency spectra of the two microphone signals
% in the algorithm and it worked out much better.
sounds=[real(fft(mic_1));real(fft(mic_2))];
% N="number of microphones"
% P="number of points"
[N,P]=size(sounds) % P=16000, N=2, in this case.
permute=randperm(P); % Generate a permutation vector.
s=sounds(:,permute); % Time-scrambled inputs for stationarity.
x=s;
mixes=sounds;
% Spheres the data (normalisation).
mx=mean(mixes');
c=cov(mixes');
x=x-mx'*ones(1,P); % Subtract means from mixes.
wz=2*inv(sqrtm(c)); % Get decorrelating matrix.
x=wz*x; % Decorrelate mixes so cov(x')=4*eye(N);
w=pi^2*rand(N); % Initialise unmixing matrix.
M=size(w,2); % M=N usually
sweep=0; oldw=w; olddelta=ones(1,N*N);
Id=eye(M);
% L="learning rate, B="points per block"
% Both are used in sep.m, which goes through the mixed signals
% in batch blocks of size B, adjusting weights, w, at the end
% of each block.
%L=0.01; B=30; sep
%L=0.001; B=30; sep % Annealing will improve solution.
%L=0.0001; B=30; sep % ...and so on
%for multiple sweeps:
L=0.0001; B=30;
for I=1:100
sep; % For details see sep.m
end;
uu=w*wz*mixes; % make unmixed sources
% Plot the two separated vectors in the frequency domain.

figure(3)
subplot(2,1,1)
plot(f,uu(1,f)), axis([0 4000 -22 22]);
title('Frequency spectra of one of the separated signals')
xlabel('Frequency (Hz)');
subplot(2,1,2)
plot(f,uu(2,f)), axis([0 4000 -22 22]);
title('Frequency spectra of the other separated signal')
xlabel('Frequency (Hz)');

% Transform signals back to time domain.
uu(2,:)=real(ifft(uu(2,:)));
uu(1,:)=real(ifft(uu(1,:)));
%disp('Playing the first of the separated vectors');
%soundsc(uu(1,:))
wavwrite(uu(1,:),"out1.wav");
% Plot the vector that is played above.

figure(4);
subplot(2,1,1)
plot(uu(1,:)), axis([0 16000 -0.12 0.12]);
title('Plot of one of the separated vectors (time domain)')
%disp('Playing the second of the separated vectors');

%soundsc(uu(2,:))
wavwrite(uu(2,:),"out2.wav");
% Plot the vector that is played above.

subplot(2,1,2)
plot(uu(2,:)), axis([0 16000 -0.12 0.12]);
title('Plot of the other separated vector (time domain)')
