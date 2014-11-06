function [ nsd ] = overlapAdd( S_A,phase_ns,fftlen,hopsize )
%OVERLAPADD overlap add method reconstruction temporel signal from
%spectrogram

% input
%       S_A : spectrogram
%       phase_ns : phasegram, two side gram
tau_ns = size(S_A,2);
S_A = [S_A(1:floor(fftlen/2),:);flipud(S_A(1:floor(fftlen/2),:))];

Spec = S_A.*exp(1i*phase_ns);
nsd=zeros((tau_ns-1)*hopsize+fftlen,1);
% weight=sig;
for ii=1:tau_ns
    start=(ii-1)*hopsize+1;
    spec=Spec(:,ii);
    nsd(start:start+fftlen-1)=nsd(start:start+fftlen-1)+real(ifft(spec,fftlen));
end

end

