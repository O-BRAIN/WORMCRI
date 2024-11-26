function [Data_filtered] = filter_fir(Data,hp,lp,Fs,fir_order);
%
% Modified klaus.linkenkaer@cncr.vu.nl, 070531.
%
%
%******************************************************************************************************************
% Purpose...
%
% Causal (feedforward) finite impulse response bandpass filter the signal 'F' with a Hamming window.
%
%
%******************************************************************************************************************
% Input parameters...
%
% Data		: data matrix (or vector), time along the 1st dimension!
% hp		: highpass corner (e.g., 8 Hz).
% lp		: lowpass corner (e.g., 1s3 Hz).
% Fs		: sampling frequency.
% fir_order	: filter order in units of seconds (to prevent an influence of sampling frequency!)
%		  Should include at least two cycles of the low-frequency component included in the passband.
%
%******************************************************************************************************************
% Default parameters (can be changed)...

% time window included in the filtering process. 
% Filter orders suitable for alpha and beta oscillations and based on:
% Nikulin. 2005. Neurosci. Long-range temporal correlations in EEG oscillations...

%fir_order = 0.25;	% Use for high time resolution and low frequency resolution of alpha oscillations.
%fir_order = 0.38	% Use for low time resolution and high frequency resolution of alpha oscillations.

%******************************************************************************************************************
% Define fihplter characteristics:

b = fir1(fir_order*Fs,[hp lp]/(Fs/2));
%b = fir1(fir_order*Fs,[lp hp]/(Fs/2));
%******************************************************************************************************************
% Filter the vector 'F' using the filter characteristics from above:

Data_filtered = zeros(size(Data)); 

for i = 1:size(Data,1);
  Data_filtered(i,:) = filter(b,1,Data(i,:));
end

SignalHilb = hilbert(Data_filtered);
phase = unwrap(angle(SignalHilb));
H=phase(1,:);
%plot(H(500:1000))