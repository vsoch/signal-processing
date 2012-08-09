% Understanding Signal Processing!
% Vanessa Sochat (signal processing noob) Stanford University, 08/2012

% If you are new to matlab, this is an example of a script with cells, and
% the idea is that you can run any cell by highlighting it (put cursor in
% it) and then press CNTRL+ENTER.  You can also highlight a single line and
% press F9 to run it, or just copy paste code from here into the terminal.

%% 1. READ IN RAW DATA EXAMPLE
% For full details, read fetFiltCoeffs_details before playing with this!
% The included timeseries data is a functional network, thresholded at
% p=0.05, derived with Independent component analysis (ICA) that is the
% dorsal default mode network! Awesome!  See data_spatial_map.png to see the
% spatial map, and data_timecourse.png for the timeseries (if you plot data
% via matlab, this is what you will see)

% Read in timecourse data for one image        
time_text = dir([ pwd '\data.txt' ]);
time_file = fopen([ pwd '\' time_text(1).name ],'r');
DATA = fscanf(time_file,'%f\t'); % Each timepoint is one TR, 2 seconds
clear time_text time_file;

%% 2. PLOT TIMESERIES
% First visualize raw signal
figure(1);
subplot(2,1,1);
plot(DATA);
title('Raw Timeseries, DATA');

%% 3. This code shows how to get from the timeseries to the 
% frequency data (via FFT)

subplot(2,1,2);
absFFT = 2*abs(fft(DATA)) / length(DATA); 
nFrames = size(DATA,1); % the # of temporal frames: 92
xfreq = 1:nFrames/2; % the different frequencies we'll plot
plot(xfreq,absFFT(xfreq));
xlabel('Frequency, cycles per scan')
ylabel('FFT Units');
title('FFT Calculation, no detrending')

%% 4. Testing calculating "high frequency noise"
% some user specified threshold (25)

HFInitIndex = 25;

% Visualize frequency data 
% PLOT 1: Green line is max, Red line is min
% PLOT 2: Blue line is "high frequency noise threshold"
figure(2); subplot(2,1,1); plot(DATA); title('Raw Timeseries, DATA'); hold on
% Find the max value of the frequency data
[max_value,max_index] = max(DATA); [min_value,min_index] = min(DATA);
line(1:length(DATA),max_value,'color','g');
line(1:length(DATA),min_value,'color','r');
subplot(2,1,2); plot(xfreq,absFFT(xfreq)); title('FFT calculated from raw timeseries');
% Add a line for the user specified threshold, above this threshold is
% considered high frequency noise
% I know the line is crooked and weird, I couldn't get it to show up when I
% did it like this:  line(25,0:1,'color','b');
line(24:25,0:1,'color','b');

% Now calculate high frequency noise as a percentage of total energy.  Too
% much high frequency noise is likely a bad component, but we will use this
% percentage as a feature. Keep track of the following
total_energy = 0; high_freq_noise = 0;
        
% Determine high frequency noise energy
for j = HFInitIndex:length(xfreq)
    high_freq_noise = high_freq_noise + (absFFT(j) * absFFT(j));
end
               
% Determine total energy and energy percent. We will start by summing the
% frequencies up until the threshold, and then add to the high_freq_noise
% to get a total, and then calculate high_freq_noise as a % of total
for j = 1:HFInitIndex
    total_energy = total_energy + (absFFT(j) * absFFT(j));
end
        
% We could check if this % noise is greater than some threshold, 
% and use that as a flag.  If the %-age high frequency
% noise was > 50, this was considered a "bad" component.
total_energy = total_energy + high_freq_noise;

% Percentage of total energy that is "high frequency"
% meaning it is above a user specified threshold (25 = X HZ)
energy_percent = 100 * (high_freq_noise / total_energy);

%% 5. Calculate mean response over time
figure(3);
mean_response = mean(DATA);
plot(DATA); title('Raw Timeseries, T with Mean Response over Time'); hold on
line(1:length(DATA),mean_response,'color','r');
clear mean_response

%% 6. Calculate number peaks, number mins, average amplitude over time
figure(4);
plot(DATA); title('Raw Timeseries (T) with Peaks (Red), Mins (Yel) and Avg Amplitude (at x=0)'); hold on
[max_value,max_index] = max(DATA); [min_value,min_index] = min(DATA);
line(1:length(DATA),max_value,'color','r');
line(1:length(DATA),min_value,'color','r');
% Find the mean response
mean_response = mean(DATA);
line(1:length(DATA),mean_response,'color','b');
% Find the local max and min values of the timeseries
[local_max,maxindex] = findpeaks(DATA);
[local_min,minindex] = findpeaks(-DATA);
% Plot max peaks
for i=1:length(local_max)
    plot(maxindex(i),local_max(i),'*','color','r');
end
for i=1:length(local_min)
    plot(minindex(i),DATA(minindex(i)),'*','color','y');
    local_min(i) = DATA(minindex(i));
end

% number of local max (peaks)
numpeaks = length(local_max);
localmaxpeaks = numpeaks / length(DATA);

% number of local mins
nummins = length(local_min);
localminpeaks = nummins / length(DATA);

% Calculate average min peak and average max peak, take difference to get
% avg amplitude
avg_amplitude = mean(local_max) - mean(local_min);
plot(mean(local_max),'*','color','r');
plot(mean(local_min),'*','color','y');

clear local_max local_min maxindex minindex pks numpeaks nummins avg_amplitude;

%% 7. Calculate distance Between Peaks, biggest jump, average jump
figure(5);
plot(DATA); title('Raw Timeseries (T) with Average min (G) and max (R) distance between peaks');

% Calculate max and min values, then find average distance between them
[max_y,maxindex] = findpeaks(DATA);
[min_y,minindex] = findpeaks(-DATA);
maxindexcopy = [0; maxindex]; minindexcopy = [0; minindex];

for i=1:length(maxindex)
    difference_max(i) =  maxindex(i) - maxindexcopy(i);
end
for i=1:length(minindex)
    difference_min(i) = minindex(i) - minindexcopy(i);
end

% Average distance between maximum peaks
mean_max_dist = mean(difference_max);
% Average distance between minimum peaks
mean_min_dist = mean(difference_min);

% Plot line for mean differences between peaks
line(0:mean_max_dist,.5,'color','b','LineStyle','-');
line(0:mean_min_dist,-.5,'color','b','LineStyle','--');

% If the min value is less than or = to max, it is first in the sorted list
% otherwise, the max value is first in the list, and the remaining values
% fluctuate min,max,min,max,etc.
minmax = [min_y; max_y];
minmaxindex = [minindex; maxindex];
[~,indices_forminmax] = sort(minmaxindex);

for i=1:length(minmax)-1
    contendervalue(i) = abs(minmax(indices_forminmax(i+1)) - minmax(indices_forminmax(i)));
    % contendervalue(i) = abs(T(minmaxindex(i+1)) - T(minmaxindex(i)));
end

% Calculate biggest jump
biggestjump = max(contendervalue);
% Calculate average jump
averagejump = mean(contendervalue);

% Plot biggest jump and average jump
% line(20,0: biggestjump,'color','b');
% line(30,0:averagejump,'color','b');
    
clear contendervalue difference_min difference_max mean_max_dist mean_min_dist maxindexcopy minindexcopy minindex maxindex;

%% 8. Calculate Frequency Bins! (translate a timeseries into buckets of Hz)
% Method 1: Vanessa's convoluted, unnecessarily complicated way:
% Hooray this is where we will use example from getFiltCoeffs!
% .01 to .1 Hz is "range of haemodynamic responses detectable with BOLD fMRI. "low" is 0 to .005 (what is our bandpass at?) (Tohnka)
% Eg., bins used by de Martino et al include:
% Power in band 0-0.008 Hz
% Power in band 0.008-0.02 Hz
% Power in band 0.02-0.05 Hz
% Power in band 0.05-0.1 Hz
% Power in band 0.1-0.25 Hz

TR = 2;             % Set the TR to 2, the seconds / timepoint
Fs = 1/TR;          % Sampling Frequency, what % of a timepoint = 1 second

% To get our "bins" we are going to bandpass filter the data according to
% the ranges, output the filtered data, sum the FFT?
% Define the list of lower and upper values for the filter bins
lowfilts = [ 0     0.1  0.15 ];
upfilts =  [ 0.1   0.15 0.2 ];

% Make a matrix to hold output data, with columns corresponding to bins, 
% and rows the number of bins
X_filtered = zeros(length(lowfilts),length(DATA));

% Calculate for each bin
figure(6)
plot_starts = [1 3 5 ]; 
for i=1:length(lowfilts); 
    subplot(6,1,plot_starts(i));
    % The "bin" is [fl:fh]
    fl = lowfilts(i); % Lower filter threshold
    fh = upfilts(i);  % Upper filter threshold
    Fc = 0.5*(fh + fl);     % Center Frequency
    % I still do not have good intiution for what a filter order is
    No = floor(Fs * 2/Fc); % Filter Order
    % For TR = 2, fl = 0.008, fh = 0.1, No = 18;
    % run getFiltCoeffs, with "data" as an empty row vector to get band B
    B = getFiltCoeffs(zeros(1,length(DATA)),Fs,fl,fh,0,No); %FIR filter Coefficients
    % Plot coefficients
    plot(B); title([ 'Filter Band for Range ' num2str(fl) ' Hz to ' num2str(fh) ' Hz' ]);
    
    % Now use this band to filter the data
    datafilt = filtfilt(B,1,DATA');
    subplot(6,1,plot_starts(i)+1)
    plot(datafilt);
    title([ 'Filtered Data for ' num2str(fl) ' Hz to ' num2str(fh) ' Hz' ])
    X_filtered(i,:) = datafilt;
end

% Now plot each of the filtered series with its FFT
figure(7);
for i=1:length(lowfilts)
    subplot(6,1,plot_starts(i));
    % First plot the new filtered data we just made, for each bin:
    plot(X_filtered(i,:));
    title([ 'Filtered Data for ' num2str(fl) ' Hz to ' num2str(fh) ' Hz' ])
    subplot(6,1,plot_starts(i)+1)
    % Calculate the FFT and sum all the values to get "amount" in the bin
    FFTfiltered = abs(fft(X_filtered(i,:)));
    plot(FFTfiltered);
   
    % If we can take any point along the FFT and multiply by the quantity on X
    % to get the amount of energy at that frequency, then to get the total amount
    % for the bin we should be able to take the entire area under the curve? 
    energy_in_bin(i) = 0;
    for j=1:length(FFTfiltered)
        energy_in_bin(i) = energy_in_bin(i) + FFTfiltered(j);
    end
end
    
%% 9. Method 2 for calculating energy in Hz bins:
% from Kaustubh Supekar, 8/9/2012

bucket_start_freq = 0.01; % starting frequency of your bucket
bucket_end_freq = 0.08;   % end frequency of your bucket
TR = 2; % TR

%ts = rand(180,1); % ts is timeseries data. replace it with actual data.
ntime = length(DATA);
nfft = ntime/2;

bucket_start_index = (nfft*bucket_start_freq)/(1/(2*TR)) ;
bucket_end_index = (nfft*bucket_end_freq)/(1/(2*TR));    
% chosen to remove components having most of the energy in the range f > 0.1 Hz
bucket_start_index = round(bucket_start_index);
bucket_end_index = round(bucket_end_index);

% compute fft
temp = abs(fft(DATA, ntime));
freq_data= temp(1:ntime/2); % just take half the spectrum (because fft is symmetric)

% compute energy in the bucket
bucket_energy = 0;
for column = bucket_start_index:bucket_end_index-1
        bucket_energy = bucket_energy + freq_data(column)^2;
end

% Plot the original timeseries, the fft, and show the bucket
figure(8);
subplot(2,1,1);
plot(DATA); title('Original timeseries');
subplot(2,1,2);
plot(freq_data); title('FFT of Original Timeseries');
% Plot bucket lines, blue is start, red is end
line(bucket_start_index,0:50,'color','b','LineStyle','-');
line(bucket_end_index,0:50,'color','r','LineStyle','-');


%% 10. Method 3 for calculating energy in Hz bins:
% from Rebecca Sawyer, 8/8/2012
 
% We have our signal DATA (in time domain) of length N, sampling frequency fs:
% Again, set the sampling Frequency (FS) based on the TR (2)
TR = 2;             % Set the TR to 2, the seconds / timepoint
fs = 1/TR;          % Sampling Frequency, what % of a timepoint = 1 second
N = length(DATA);

X = fft(DATA,4*N); % This will give you better resolution in the frequency domain.
% Now calculate the f (in Hz)
f = (1:4*N-1)/(4*N-1)*fs;

% To find the appropriate indices for your bins, you can either do it manually 
% by just looking at the values of f or you can try to determine it automatically 
% doing something like:
f1 = 0;   % Hz
ind1 = 1; % because f1=0
f2 = 0.1; % Hz
diff = abs(f-f2);
[m ind2] = min(diff);

% Then to get the area:
area = trapz(X(ind1:ind2));

% From Rebecca:
%  Frequency at each point of FFT should be 
% (0:N-1)*360/(N-1) 
% (-(N-1)/2:(N-1/2))*180/((N-1)/2) for fftshift data (-pi to pi)
% For this data since we just have 0 to pi, we would want:
% (0:N-1)*360/(N-1)
% N is the number of samples (length of FFT output)


%% Experimenting with FFT and IFFT
% Katie says that IFFT should reverse FFT
% View timeseries data, "DATA"
figure(10); 
subplot(5,1,1);
plot(T); 
title('Raw timeseries data');

% Calculate FFT, without cutting in half
FFT = abs(fft(DATA,length(DATA)));
subplot(5,1,2);
plot(FFT)
title('FFT of T, with FFT units');

% Test inverse FFT function
subplot(5,1,3);
IFFTy = ifft(FFT);
plot(IFFTy);
title('Inverse of FFT with IFFTY - not the original T!')

% Whoops, try taking FFT again without absolute value
subplot(5,1,4);
FFT2 = fft(DATA,length(DATA));
plot(FFT2);
title('FFT without taking absolute value - zomg what is happening!')
IFFTy2 = ifft(FFT2);
subplot(5,1,5);
plot(IFFTy2);
title('Inverse of FFT without taking absolute value - the original T!');