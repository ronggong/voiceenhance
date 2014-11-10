function [output]=kalmancolor(ns,nf,fs,it)

% OUTPUT=KALMANSIGNALDENOISER(ns,nf,fs,it)
% this purpose of this function is to demonstrate the capability of kalman
% filter for denoising noisy speech (corrupted by colored noise). Kalman
% filtering of noisy speech usually have two steps: 
% 1 . iteratively estimating the AR parameters of noisy segment
% 2 . Filtering the segment
%
% ARGUMENTS
% ns :  noise contaminated speech
% nf :  noise signal, stationary
% fs :  Sampling frequency which should be the same for both signals
% it :  iterative time
%
% Output is the denoised speech signal
%
% Required functions:
% SEGMENT
%Sep-04
%Esfandiar Zavarehei
%edited by Rong GONG 10/10/2014
if nargin < 4
    it = 3;
end

W=fix(.025*fs);                                         %Window length is 25 ms
SP=1;                                   
SpecNs = 13;                                            % noisy speech LPC estimation order
SpecNf = 13;                                            % noise LPC estimation order
SpecP= SpecNs + SpecNf;
Window=ones(W,1);

y=segment(ns,W,SP,Window);                              % noisy frames

H = zeros(1,SpecP); H(1,SpecNs)=1; H(1,end)=1;          % measurement matrix

NsUpper=[zeros(SpecNs-1,1) eye(SpecNs-1)];
NfUpper=[zeros(SpecNf-1,1) eye(SpecNf-1)];
I=eye(SpecP);

[a_ns, Q_ns]=lpc(y,SpecNs);                              % initial noisy speech AR parameters
[a_nf, Q_nf]=lpc(nf,SpecNf);                     % initial noise AR parameters

Q = zeros(SpecP);
P=diag([repmat(Q_ns(1),1,SpecNs) repmat(Q_nf,1,SpecNf)]);   % inital P

o=zeros(1,W*size(y,2));                                     % allocating memory
o(1:SpecNs)=y(1:SpecNs,1)';

hwb = waitbar(0,'Please wait...','Name','Processing');
start=SpecNs+1;

Sp = 0;
t=SpecNs+1;
for n=1:size(y,2)
    waitbar(n/size(y,2),hwb,['Please wait... ' num2str(fix(100*n/size(y,2))) ' %'])
    t_old = t;
    Sp_old = Sp;
    for kk = 1:it
        A_ns=[NsUpper; fliplr(-a_ns(n,2:end))];
        A_nf = [NfUpper; fliplr(-a_nf(1,2:end))];
        A = [A_ns, zeros(SpecNs,SpecNf);zeros(SpecNf,SpecNs),A_nf];
        Q(SpecNs,SpecNs) = Q_ns(n); Q(SpecP,SpecP) = Q_nf;
        for i=start:W
            S_=A*Sp;                                        % one step predicted estimate
            P_=A*P*A'+Q;                                    % error covariance matrix     
            K=(P_*H')/(H*P_*H');                            % kalman gain
            Sp=S_+K*(y(i,n)-H*S_);
            o(t-SpecNs+1 : t)=Sp(1:SpecNs)';                %Notice that the previous SpecP-1 output samples are updated again
            P=(I-K*H)*P_;
            t=t+1;
        end
        start=1;
        if kk < it
            t = t_old;
            Sp = Sp_old;
        end
        [a_ns(n,:), Q_ns(n)] = lpc(o((n-1)*W+1 : n*W),SpecNs);  % update noisy AR parameters
    end
end
close(hwb)
output=o;

[ full_S_nofilter,freqtemp_no,middle_frames_no ] = STFT_KI( ns,fs,W,W,0 );
[ full_S,freqtemp,middle_frames,tau,phase_S ] = STFT_KI( o,fs,W,W,0 );
nFreq = floor(size(full_S,1)/2)+1;
full_S_nofilter = full_S_nofilter(1:nFreq,:);
full_S = full_S(1:nFreq,:);

% for plot
% figure
% imagesc(middle_frames_no,freqtemp_no,20*log10(full_S_nofilter))
% caxis([-20 20])
% colorbar
% set(gca,'YDir','normal')
% xlabel('frame')
% ylabel('frequency Hz')
% title('noisy speech')
% 
% figure
% imagesc(middle_frames,freqtemp,20*log10(full_S))
% caxis([-20 20])
% colorbar
% set(gca,'YDir','normal')
% xlabel('frame')
% ylabel('frequency Hz')
% title(['output of iterative AR coefficients estimation Kalman filter, iteration time = ',num2str(3)])

