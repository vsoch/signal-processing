% ------------------------------------------------------------------------
% getFiltCoeffs TUTORIAL HEADER
%-------------------------------------------------------------------------
% getFiltCoeffs is a robustly commented version of getFiltCoeffs
% with the idea of understanding how we get the frequency (Hz values) 
% corresponding with some range of a raw timeseries.  This commentary was
% written by a graduate student with no background in signal processing in
% an effort to make the calculation steps transparent.

% Step 1: What is the data? You have some data that you want to filter, 
% perhaps for highpass, bandpass, or lowpass. This is some power value that 
% is changing over time.  Let's define this timeseries, DATA, as a column 
% vector of T values, each corresponding to one measured timepoint.

% Step 2: How often do we want to sample?  "Sampling" the data means how
% often we reach in and "grab" a measurement.  For example, we might have
% 10 timepoints, each of which is 2 seconds, but decide that we want to sample
% for each 1 second.  This is defined as our "sampling frequency" that we
% we will call the variable Fs.  For example, let's say that our DATA
% timeseries has 10 values (T=10) and each timepoint value is two seconds. 
% (In neuroimaging, this value is called the TR).  If we wanted to sample every
% second, to get the sampling frequency we would calculate:
%
% FS = 1/TR
% where "TR" = the number of seconds represented at each data value.  All
% this is saying is, if I were to look at one timepoint, DATA(T), what
% percentage of that timepoint corresponds to one second?  So if the TR =
% 2, meaning that each timepoint is 2 seconds, the sampling frequency FS =
% 1/TR would = 1/2, meaning that I can break any timepoint into halves to
% break it into seconds. I'm thinking of this "sampling frequency" as a way
% to describe signals (that might have different amounts of times corresponding
% to each DATA(T) ) with a common unit, the second.

% Step 3: What frequency levels am I interested in?  If I am doing any sort
% of filtration, once we are talking about the common units of Hz, I would
% want to keep data within some range (bandpass), filter out signal above
% some high cutoff (highpass) or a low cutoff (lowpass).  Based on reading
% papers, or understanding levels of signal, I need to decide upon a range
% of values that I am interested in.  

% Step 5: What do I do with this range? depends on the filter or method
% that you are applying.  For...
%  bandpass: You want to keep signal between fl and fh.  In this case, you
%  should calculate the center frequency, Fc, which is 0.5*(fh+fl);
%  lowpass: You want to keep everything below fl
%  highpass: You want to keep everything above hp

% Step 5: We now need to specify a filter order, specified by
% floor(Fs * 2/Fc)
% I do not have strong intuition for this calculation.  I read a bit about 
% this, and it is defined as "the maximum delay, in samples used in creating 
% each output sample."  We will call this variable "No"
% Before looking at the calculation, I am thinking that we have some number 
% of actually sampled datapoints the length of DATA, T, and each of those 
% datapoints is ACTUALLY 2 seconds (the TR).  So our sampling frequency, Fs, 
% says that 1/TR of each timepoint is one second.  So maybe what "filter order
% is getting at is the difference between this original number of points,
% and the points we are actually going to use, but normalized by the center 
% of the range of frequencies that we are interested in (Fc).  Thus, it makes
% sense to do Fs * 2 (sampling frequency X TR) to get back to our "base unit",
% 1 second, and then normalize by the center frequency.  For now I wlll
% accept that abstract (likely wrong :P) rationale and keep an ear open for
% what "filter order" means when I watch Signal Processing Lecture.

% Step 5: We now have our data, a range or cutoff that we are interested in, 
% and the sampling frequency that lets us understand our data in a common unit, 
% the second, and some filter order.  The next step

% In Summary:
% 
% DATA: raw timeseries, a column vector of values recorded over time
% T: the length of the data, the number of timepoints
% TR: the number of seconds represented in each timepoint
% Fs: the sampling frequency, how often we dip our spoon into the honey pot
%     1/TR
% fl: Low frequency "cutoff"
% fh: High frequency "cutoff"
% Fc: "center frequency" of range between fl and fh, defined as
%     0.5*(fh+fl)
% No: filter order: max delay in samples used in creating output sample
%     floor(Fs * 2/Fc)

% Step 6: We would call:
%     B = getFiltCoeffs(zeros(1,length(DATA)),Fs,fl,fh,0,No);
% Corresponding to:
%     filtwts = getFiltCoeffs(data,srate,locutoff,hicutoff,epochframes,filtorder,revfilt)
% Continue to code below to see continued walk through!


%--------------------------------------------------------------------------
% GETFILTCOEFFS.m Original Script Header:
%--------------------------------------------------------------------------
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]

% Author: Scott Makeig, Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1997 
% Explanatory commenting added by Vanessa Sochat, Stanford 2012

% eegfilt() -  (high|low|band)-iass filter data using two-way least-squares 
%              FIR filtering. Multiple data channels and epochs supported.
%              Requires the MATLAB Signal Processing Toolbox.
% Usage:
%  >> [smoothdata] = eegfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                             epochframes,filtorder);
%
% See also: firls(), filtfilt()
%

% Copyright (C) 4-22-97 from bandpass.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% As called by example above: B = getFiltCoeffs(zeros(1,length(DATA)),Fs,fl,fh,0,No);
function filtwts = getFiltCoeffs(data,srate,locutoff,hicutoff,epochframes,filtorder, revfilt)

% Variables corresponding to example above
% It is not clear to me why we give an empty vector, but it yields the same
% result.  It just needs to be a row vector, which is salient, because our
% original timeseries I put into a column!
% data = zeros(1,length(DATA));
% srate = Fs;
% locutoff = fl;
% hicutoff = fh;
% epochframes = 0;
% filtorder = No;

% We need minimally the first four args for eegfit
if nargin<4
    fprintf('');
    help eegfilt
    return
end

% We also need the signal processing toolbox
if ~exist('firls')
   error('*** eegfilt() requires the signal processing toolbox. ***');
end

% Frames is the length of the data, chans is the number of channels
% For our timeseries, we have 1 channel, and 165 frames
[chans frames] = size(data);

% Get angry at the user if we don't have a row vector
if chans > 1 & frames == 1,
    help eegfilt
    error('input data should be a row vector.');
end

% The nyquist frequency is "half the sampling frequency of a discrete signal 
% processing system. Aka, folding frequency. [wikipedia]
% Our srate is 0.5 (meaning 1/2 of a timepoint = 1 second) so we are just
% getting at half of that. **NEED TO WRITE ABOUT WHY IMPORTANT HERE!
nyq = srate*0.5;  % Nyquist frequency
%MINFREQ = 0.1/nyq;
MINFREQ = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter 
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

% Obviously the low needs to be less than the high
if locutoff>0 & hicutoff > 0 & locutoff > hicutoff,
    error('locutoff > hicutoff ???\n');
end
% And both values poitive
if locutoff < 0 | hicutoff < 0,
   error('locutoff | hicutoff < 0 ???\n');
end

% I don't think that we are allowed to "break" the resolution of our
% sampling frequency, nyquist? **BETTER UNDERSTANDING HERE
% VANESSA STOPPED HERE
if locutoff>nyq,
    error('Low cutoff frequency cannot be > srate/2');
end

if hicutoff>nyq
    error('High cutoff frequency cannot be > srate/2');
end

if nargin<6
   filtorder = 0;
end
if nargin<7
   revfilt = 0;
end

if isempty(filtorder) | filtorder==0,
   if locutoff>0,
     filtorder = minfac*fix(srate/locutoff);
   elseif hicutoff>0,
     filtorder = minfac*fix(srate/hicutoff);
   end
     
   if filtorder < min_filtorder
        filtorder = min_filtorder;
    end
end

if nargin<5
	epochframes = 0;
end
if epochframes ==0,
    epochframes = frames;    % default
end
epochs = fix(frames/epochframes);
if epochs*epochframes ~= frames,
    error('epochframes does not divide frames.\n');
end

if filtorder*3 > epochframes,   % Matlab filtfilt() restriction
    fprintf('eegfilt(): filter order is %d. ',filtorder);
    error('epochframes must be at least 3 times the filtorder.');
end
if (1+trans)*hicutoff/nyq > 1
	error('high cutoff frequency too close to Nyquist frequency');
end;

if locutoff > 0 & hicutoff > 0,    % bandpass filter
    if revfilt
         fprintf('eegfilt() - performing %d-point notch filtering.\n',filtorder);
    else 
         %fprintf('eegfilt() - performing %d-point bandpass filtering.\n',filtorder);
    end; 
    %fprintf('            If a message, ''Matrix is close to singular or badly scaled,'' appears,\n');
    %fprintf('            then Matlab has failed to design a good filter. As a workaround, \n');
    %fprintf('            for band-pass filtering, first highpass the data, then lowpass it.\n');

    f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
    %fprintf('eegfilt() - low transition band width is %1.1g Hz; high trans. band width, %1.1g Hz.\n',(f(3)-f(2))*srate, (f(5)-f(4))*srate/2);
    m=[0       0                      1            1            0                      0]; 
elseif locutoff > 0                % highpass filter
 if locutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq));
 end
 %fprintf('eegfilt() - performing %d-point highpass filtering.\n',filtorder);
 f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1]; 
 %fprintf('eegfilt() - highpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
 m=[   0             0                   1      1];
elseif hicutoff > 0                %  lowpass filter
 if hicutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - lowpass cutoff freq must be > %g Hz',MINFREQ*nyq));
 end
 %fprintf('eegfilt() - performing %d-point lowpass filtering.\n',filtorder);
 f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1]; 
 %fprintf('eegfilt() - lowpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
 m=[     1           1              0                 0];
else
    error('You must provide a non-0 low or high cut-off frequency');
end
if revfilt
    m = ~m;
end;
filtwts = firls(filtorder,f,m); % get FIR filter coefficients