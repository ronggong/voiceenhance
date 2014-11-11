


function [nsd] = martinMini94(ns,fs,ov)

% Spectral Subtraction Based on Minimum Statistics
% input
%      ns: noisy speech
%      fs
%      ov: method of calculating oversubtraction factor, 0-same as in
%           voicebox specsubm.m; 1-Berouti
% output nsd
% Rong GONG, 11/11/2014
if nargin < 3
    ov = 0;
end

%% stft
fftlen = 2^nextpow2(0.032*fs);
hopsize = fftlen/4;
[ S,freq,time_ns,tau,phase ] = STFT_KI( ns,fs,fftlen,hopsize,0 );
nFreq = floor(size(S,1)/2)+1;
S = S(1:nFreq,:);
%% power spectrum
alpha = 0.9231;             % smooth factor
gamma = 0.8187;             % smooth factor
subf = 0.02;                % spectral floor constant

Sp = S.^2;                  % power spectrum
Sps = zeros(nFreq,tau);     % smooth power spectrum used in noise estimation
% Spf = Sps;                  % smooth power spectrum used
Sps(:,1) = Sp(:,1);
Nps = zeros(nFreq,tau);     % noise power spectrum estimate
omin = 1.5;                 % bias factor
SNR = zeros(nFreq,tau);     % SNR
delta = zeros(nFreq,tau);   % oversubtraction
Y =zeros(nFreq,tau);        % output
Q =zeros(nFreq,tau);
osf=4*(1+(0:fftlen/2).'*fs/(fftlen*400)).^(-1);

qL = 1;                     % minimum SNR
qH = 100;                   % maximum SNR
logqL = 10*log10(qL);
logqH = 10*log10(qH);
deltaL = 1;                 % oversubtraction minimum
deltaH = 4;                 % oversubtraction maximum

% noise power buffer intialisation
D = round(0.8*fs/hopsize);      % minimum window
W = 4;                          % number of subwindows

if mod(D,W)
    D = D + mod(D,W);
end
M = D/W;                        % subwindow length
cbi=1 ;                         % index to circular buffer
cbmin=Sps(:,1) ;                % current minimum of the circular buffer
lsmin=Sps(:,1) ;                % current minimum of the last segment
cb=repmat(fftlen/2,nFreq,W) ;   % circular buffer of sps history
lsc= M;                         % length of the last segment

Nps(:,1)=Sps(:,1);
for ii = 2:tau
    % smoothing
    Sps(:,ii) = alpha*Sps(:,ii-1)+(1-alpha)*Sp(:,ii);
    Sp(:,ii) = gamma*Sp(:,ii-1)+(1-gamma)*Sp(:,ii);
    
    % noise power estimation
    if lsc >= M;                % last frame of the subwindow
        cb(:,cbi) = lsmin;
        cbi = mod(cbi,W)+1;
        cbmin = min(cb,[],2);
        lsmin = Sps(:,ii);
        lsc = 1;
    else
        lsmin = min(Sps(:,ii),lsmin);
        lsc = lsc+1;
    end
    Nps(:,ii)=omin*min(lsmin,cbmin);   % new noise estimate
    
    % oversubtraction factor
    if ov
        SNR(:,ii) = 10*log10((Sps(:,ii)-min([Sps(:,ii),Nps(:,ii)],[],2))./Nps(:,ii));  % SNR v2
        %----------------------------------------------------------------------
        indinf5 = (SNR(:,ii)<-5);                                                   % noise greater than speech
        indinf20 = (SNR(:,ii)<=20&SNR(:,ii)>=-5);
        indsup20 = (SNR(:,ii)>20);
        delta(indinf5,ii) = 5;
        delta(indinf20,ii) = 4 - 3*SNR(indinf20,ii)/20;
        delta(indsup20,ii) = 1;                                                     % speech is greater than noise
        %----------------------------------------------------------------------
    else
        delta(:,ii)=0.9*delta(:,ii-1)+(1-0.9)*(1+osf.*Nps(:,ii)./(Nps(:,ii)+Sps(:,ii))); % delta v3
    end
    
    Q(:,ii) = 1-sqrt(delta(:,ii).*Nps(:,ii)./Sp(:,ii));
    for jj = 1:size(S,1);
        Y(jj,ii) = max(S(jj,ii)*Q(jj,ii),subf*sqrt(Nps(jj,ii)));
    end
end

nsd = overlapAdd( Y,phase,fftlen,hopsize );        % reconstruction

% h_nsd = audioplayer(nsd,fs);                       % play sound
% play(h_nsd)
% 
% figure
% plot(time_ns,Sps(87,:),time_ns,Nps(87,:))
% % 
% figure
% imagesc(time_ns,1:size(S,1),20*log10(S));
% colorbar
% set(gca,'YDir','normal')
% xlabel('time (s)')
% ylabel('freq')
% title('ouput speech, voicebox oversubtraction factor')


end


