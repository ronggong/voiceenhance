

function [nsd,fs_ns] = implicitWiener(noisyspeech,noise_file,method)
% implementation Enhancement and Bandwidth Compression of Noisy Speech
% input
%       noisyspeech
%       noise_file : noise sample, should be the same format (fs) to noisy
%                    speech, 
%       method : 0 - eq 26b, 1 - MLE eq 14
% output
%       nsd : noisy signal de-noised
if nargin < 3
    method = 0;
end
% 
% if nargin < 4
%     flagVad = 0;
% end
%% stft

[nf,fs_nf] = audioread(noise_file);
[ns,fs_ns] = audioread(noisyspeech);

fftlen = 0.03*fs_nf;
hopsize = 0.015*fs_nf;

[ S_nf,S_freq,time_nf,tau_nf,phase_nf ] = STFT_KI( nf,fs_nf,fftlen,hopsize,0 );
[ S_ns,~,time_ns,tau_ns,phase_ns ] = STFT_KI( ns,fs_ns,fftlen,hopsize,0 );
nFreq = floor(fftlen/2)+1;
S_nf = S_nf(1:nFreq,:);
S_ns = S_ns(1:nFreq,:);
S_A = zeros(size(S_ns,1),tau_ns);


% figure
% imagesc(time_nf,S_freq,20*log10(S_nf));
% set(gca,'YDir','normal')
% xlabel('time (s)')
% ylabel('freq')
% title('noise signal magnitude')
%
figure
imagesc(time_ns,S_freq,20*log10(S_ns));
set(gca,'YDir','normal')
xlabel('time (s)')
ylabel('freq')
title('noisy speech magnitude')

%% procedure
a_nf = mean(S_nf.^2,2);                                               % expectation of noise power spectrum
P_ns = S_ns.^2;                                                       % power spectrum of noisy signal
% N_R = bsxfun(@minus,S_nf,a_nf);                                     % noise residual
% max_NR = max(N_R,[],2);                                             % maximum value of noise residual
% c = 10^(-1.5);
alpha = 2;
for jj = 2:tau_ns-1
    switch method
        case 0
            temp = (P_ns(:,jj) - alpha*a_nf);
            temp(temp<0) = 0;
            S_A(:,jj) = sqrt(temp);                    % eq 26b
        case 1
            temp = P_ns(:,jj) - a_nf;
            temp(temp<0) = 0;
            S_A(:,jj) = 0.5*S_ns(:,jj) + 0.5*sqrt(temp);        % MLE
    end
    
    
%     if flagVad
%         %% additional signal attenuation during non speech activity
%         T = 20*log10(sum(S_A(:,jj)./a_nf)/size(S_A,1));
%         if T < t_th                                                   % frame jj is indicated as noise
%             S_A(:,jj) = c * S_ns(:,jj);                              % attenuation -30dB
%             
%             a_nf = (a_nf*tau_nf + S_ns(:,jj))/(tau_nf+1);            % updata average noise
%             max_NR = max(max_NR, S_ns(:,jj) - a_nf);                 % updata maximum of noise residual
%             tau_nf = tau_nf + 1;
%         end
%     end
    
end

%% reconstruct
nsd = overlapAdd( S_A,phase_ns,fftlen,hopsize );

figure
imagesc(time_ns,S_freq,20*log10(S_A));
set(gca,'YDir','normal')
xlabel('time (s)')
ylabel('freq')
title(['output spectrogram '])

end
