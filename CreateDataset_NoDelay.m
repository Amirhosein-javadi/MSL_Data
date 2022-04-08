clc;
clear all;
close all;
Fc = 6e9; % 10 Gig
C = 3e8;
fd = 1000;
v = (fd*C/Fc)*3.6;
frame_size = 16;
L = 4;
P = 4;
bitrate = 1e9;
t = (0:1/bitrate: (frame_size+L+P-1)/bitrate).';
alpha_arr = [0.9,0.85,0.8,0.75,0.7];
Tap_Number = 4;
SNR_max = 15;
%%
Train_iter = 1e5; % number of dataset columns
data_input_train = zeros(Train_iter*size(alpha_arr,2), 4*2 + 1); % 4 pilot and alpha in SEFDM
data_output_train = zeros(Train_iter*size(alpha_arr,2), 4 + 4); % amplitude and doopler freq
for j = 1:size(alpha_arr,2)
    alpha = alpha_arr(j);
    F = SEFDM_Modulation(alpha,size(t,1)); % ...... |||| ...... %
    for i =1:Train_iter
        SNR = rand()*SNR_max+15;
        Delay = zeros(1,4);
        Doppler_Frequency = rand(4,1)*fd;
        % rician
        Amplitude = raylrnd(1,4,1);  % parameter ---> tecnical report
        data_output_train((j-1)*Train_iter+i,:) = [Doppler_Frequency;Amplitude];
        S = (randi([0,1],frame_size,1)*2-1) + 1i*(randi([0,1],frame_size,1)*2-1);
        Sx = [S(end-3:end);S(1:6);ones(4,1);S(7:16)];  % 4 cp, 6 signal, 4 pilot, 10 signal
        S_TX = conj(F) * Sx;
        S_RX = channel(t,S_TX,Tap_Number,Amplitude,Doppler_Frequency,Delay);
        Noisy_Signal = awgn(S_RX,SNR,'measured');
        data_input_train((j-1)*Train_iter+i,:) = [real(Noisy_Signal(11:14));imag(Noisy_Signal(11:14));alpha];
    end
end
input_text = ["pilot1_r","pilot2_r","pilot3_r","pilot4_r",...
              "pilot1_I","pilot2_I","pilot3_I","pilot4_I",...
              "alpha"];
output_text = ["Doppler_Freq_1","Doppler_Freq_2","Doppler_Freq_3","Doppler_Freq_4",...
              "Amplitude_1","Amplitude_2","Amplitude_3","Amplitude_4"];   
csvwrite('data_output_train.txt',data_output_train)
csvwrite('data_input_train.txt',data_input_train)
%%
Validation_iter = 1e4;
data_input_validation = zeros(Validation_iter*size(alpha_arr,2), 4*2 + 1);
data_output_validation = zeros(Validation_iter*size(alpha_arr,2), 4 + 4);
for j = 1:size(alpha_arr,2)
    alpha = alpha_arr(j);
    F = SEFDM_Modulation(alpha,size(t,1));
    for i =1:Validation_iter
        SNR = rand()*SNR_max+15;
        Delay = zeros(1,4);
        Doppler_Frequency = rand(4,1)*fd;
        Amplitude = raylrnd(1,4,1);
        data_output_validation((j-1)*Validation_iter+i,:) = [Doppler_Frequency;Amplitude];
        S = (randi([0,1],frame_size,1)*2-1) + 1i*(randi([0,1],frame_size,1)*2-1);
        Sx = [S(end-3:end);S(1:6);ones(4,1);S(7:16)];  % 4 cp, 6 signal, 4 pilot, 10 signal
        S_TX = conj(F) * Sx;
        S_RX = channel(t,S_TX,Tap_Number,Amplitude,Doppler_Frequency,Delay);
        Noisy_Signal = awgn(S_RX,SNR,'measured');
        data_input_validation((j-1)*Validation_iter+i,:) = [real(Noisy_Signal(11:14));imag(Noisy_Signal(11:14));alpha];
    end
end
csvwrite('data_output_validation.txt',data_output_validation)
csvwrite('data_input_validation.txt',data_input_validation)
%%
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