function [ full_S,freqtemp,middle_frames,tau,phase_S ] = STFT_KI( aud,fs,fftlen,hopsize,obsmeth )
%STFT_KI is a function to calculte the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OBESERVATION INITIALIZATION : AUDIO ANALYSIS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  The matrix xwin is built. 
%               the axis y of xwin is frame axis
%               samples are tiled by the function "buffer"

% output:       xwin, windowed tiled sample matrix

%               tau,  number of frames

% Default ANTESCOFO STAND ALONE analysis parameters
% gamma = -1;
% fftlen = 2048;  % = winsize in Antescofo stand alone
% hopsize = 512;
% nofharm = 10;
% Default MAX PATCH analysis parameters
%  gamma = -1.;
%  fftlen = 4096;
%  hopsize = 1024;
%  nofharm = 10;

if nargin < 5
    obsmeth = 0;
end

%TODO: include parameters pedal & pedalweight

%base = 440; % basis frequency (in Hz) corresponding to 69 MIDI
%ped_weight = 0.;
%obsmeth = 0;    % Observation method, 0=FFT, 1=Constant-Q

% % Load audio
% [aud, fs] = sfread(audiofile);  % "aud" contains audio file
% aud=aud*10.0;
% % crop audio
% if exist('durationAudio','var') && ~isempty(durationAudio)
%     aud = aud(1+floor(fs*durationAudio(1)) : min(end,ceil(fs*durationAudio(end))) );
% end

% Spectral processing inits
%fftlen = 4096;
%hopsize   = 512; %1024;
%hopsize = 1024;
%winsize=fftlen;
winsize = fftlen;
%winsize = 2048;  % WARNING
win = hann(winsize);    % TODO: check analysis window

if obsmeth == 1
    % inits for constant-Q observation
    minFreq = 65.4;maxFreq = 7902;cqthresh=0.0054;cqbins=48;
    sparKernel= sparseKernel(minFreq, maxFreq, cqbins, fs, cqthresh);
    cqQ= 1/(2^(1/cqbins)-1); k=1:ceil(cqbins*log2(maxFreq/minFreq));
    freqtemp = minFreq*2.^((k-1)/cqbins);
    %vdim = length(freqtemp);
elseif obsmeth==0
    freqtemp =  fs/2*linspace(.001,1,fftlen/2+1); % Spectrogram frequencies
    %vdim=fftlen/4;   % in case of FFT only
end

[xbuf,~] = buffer(aud, winsize, winsize-hopsize, 'nodelay'); % return only full frames
% winsize-hopsize is the number of samples overlap
middle_frames = (.5+(1:size(xbuf,2))) * hopsize/fs; % audio frame center time (in s)
% this is not exact center time

winmat = repmat(win, 1, size(xbuf,2));
xwin = xbuf.*winmat;

tau = size(xwin,2); % number of analysis frames

if obsmeth == 1
    for t = 1:tau
        full_S(:,t) = abs(constQ(xwin(:,t)', sparKernel))';
    end
else
    fft_S = fft(xwin,fftlen);
    full_S = abs(fft_S);
    phase_S = atan2(imag(fft_S),real(fft_S));
end

end

