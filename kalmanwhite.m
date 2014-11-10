function [output]=kalmanwhite(ns,nf,fs,it)

% OUTPUT=KALMANSIGNALDENOISER(NOISY,CLEAN,FS)
% this purpose of this function is to demonstrate the capability of kalman
% filter for denoising noisy speech (corrupted by white noise). Kalman
% filtering of noisy speech usually have two steps: 
% 1 . iterative estimating the AR parameters of noisy speech segment
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
% edited by Rong GONG 10/11/2014
if nargin<4
    it = 3;
end
W=fix(.025*fs);                         %Window length is 25 ms
SP=1;                                    
SpecP=13;
Window=ones(W,1);
y=segment(ns,W,SP,Window);
R=var(nf);                              % variance white noise

H=[zeros(1,SpecP-1) 1];                 % measurement matrix
FUpper=[zeros(SpecP-1,1) eye(SpecP-1)];
I=eye(SpecP);

a1 = zeros(size(y,2),SpecP+1);
Q1 = zeros(size(y,2),1);
[a, Q]=lpc(y,SpecP);
P=diag(repmat(R(1),1,SpecP));

o=zeros(1,W*size(y,2));                 % allocating memory 
o(1:SpecP)=y(1:SpecP,1)';
hwb = waitbar(0,'Please wait...','Name','Processing');
Sp=y(1:SpecP,1);

start=SpecP+1;
t=SpecP+1;
for n=1:size(y,2)
    waitbar(n/size(y,2),hwb,['Please wait... ' num2str(fix(100*n/size(y,2))) ' %'])
    t_old = t;
    Sp_old = Sp;
    for kk = 1:it
        A=[FUpper; fliplr(-a(n,2:end))];
        for i=start:W
            S_=A*Sp;                                % one step predicted estimate
            P_=A*P*A'+H'*Q(n)*H;                    % error covariance matrix
            K=(P_*H')/(H*P_*H' + R);                % kalman gain
            Sp=S_+K*(y(i,n)-H*S_);
            o(t-SpecP+1 : t)=Sp';                   %Notice that the previous SpecP-1 output samples are updated again
            P=(I-K*H)*P_;
            t=t+1;
        end
        start=1;
        if kk < it
            t = t_old;
            Sp = Sp_old;
        end
        [a(n,:), Q(n)]=lpc(o((n-1)*W+1:n*W),SpecP); % update lpc on filtered signal
    end
end

close(hwb)
output=o;

