clc;
clear all;
close all;
% filename = 'M.csv';
% M = csvread(filename);
% csvwrite('M.txt',M)

iter = 10; % number of dataset columns
data_input = zeros(iter, 4*2 + 1); % 4 pilot and alpha in SEFDM
data_output = zeros(iter, 4 + 4); % amplitude and doopler freq
fd = 10;
frame_size = 16;
L = 4;
P = 4;
bitrate = 1e9;
t = (0:1/bitrate: (frame_size+L+P-1)/bitrate).';
alpha = 1;
F = SEFDM_Modulation(alpha,size(t,1));
Tap_Number = 4;
SNR = 100;
for i =1:iter
    Delay = zeros(1,4);
    Doppler_Frequency = rand(4,1);
    Amplitude = raylrnd(1,4,1);
    data_output(i,:) = [Doppler_Frequency;Amplitude];
    S = (randi([0,1],frame_size,1)*2-1) + 1i*(randi([0,1],frame_size,1)*2-1);
    Sx = [S(end-3:end);S(1:6);ones(4,1);S(7:16)];  % 4 cp, 6 signal, 4 pilot, 10 signal
    S_TX = conj(F) * Sx;
    S_RX = channel(t,S_TX,Tap_Number,Amplitude,Doppler_Frequency,Delay);
    Noisy_Signal = awgn(S_RX,SNR,'measured');
    data_input(i,:) = [real(Noisy_Signal(11:14));imag(Noisy_Signal(11:14));alpha];
end

csvwrite('data_output.txt',data_output)
csvwrite('data_input.txt',data_input)

function F = SEFDM_Modulation(alpha,n)
    F = zeros(n,n);
    for i = 0:n-1
        for j = 0:n-1
            F(i+1,j+1) = exp(-1i*i*j*alpha*2*pi/n)/(n^0.5);
        end
    end
end

function y = channel(t,x,Tap_Number,Amplitude,Doppler_Frequency,Delay)
    y = zeros( size( t ));
    sample_rate= t(2) - t(1);
    for i = 1:Tap_Number
         tau = round(Delay(i)*sample_rate); % should change 
         y(1:end-tau) = y(1:end-tau) + Amplitude(i) * exp(1i*2*pi*Doppler_Frequency(i)*t(1+tau:end)).*x(1+tau:end);
    end
end