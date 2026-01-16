clc;
clear
close all;
%% www.onlinebme.com
% filter bank CCA (FBCCA) in SSVEP frequency detection (demo code)
% by    Mohammad Norizadeh Cherloo,
%       Homa Kashefi Amiri,
%       Amir Mohammad Mijani,
%       Liang Zhan,
%       Mohammad Reza Daliri

%% define prameters (Fs,data length, channels and Nh,...)
Fs=250;% sampling rate

% -------------------------------------------------------------------------
% Electrode index mapping (from your table)
% -------------------------------------------------------------------------
idx_Pz  = 48;
idx_PO5 = 54;
idx_PO3 = 55;
idx_POz = 56;
idx_PO4 = 57;
idx_PO6 = 58;
idx_O1  = 61;
idx_Oz  = 62;
idx_O2  = 63;

% 9 channels: Chen 2015 set
chn = [idx_Pz idx_PO5 idx_PO3 idx_POz idx_PO4 idx_PO6 idx_O1 idx_Oz idx_O2];

% number of harmonics
Nh=5;

% data lenght in seconds (0.5,1,1.5,2,2.5 and 3 were considered in our study)
duration=2.5;
time= linspace(0,6,1500);
position= find(time>=0.5 & time<=0.5+duration); % index of EEG signal

% design a band-pass butterworth filter
[b,a]= butter(3,[8 90]/(Fs/2), 'bandpass');

% load frequency-phase information of stimuli
load('dataset\SSVEP-BCI-Data\Freq_Phase.mat')
fstim= freqs;
% build label for each stimulus which will be used for evaluatoin
freqs= repmat(freqs,1,6);
y_true= repmat(1:40,1,6);

% define frequency bands and their weights
freq_bands= [(1:9)*8; ones(1,9)*90];
freq_bands(1,1:end)=freq_bands(1,1:end)-2;

sigma=1;
Nsb=[1:size(freq_bands,2)]';
w= exp(-Nsb/ (2*(sigma^2)) );

%% Filter bank schemes over 8–88 Hz

sigma = 1;   % for exponential weights

% ----- M1: equal 8 Hz bands from 8 to 88 Hz -----
edges = 8:8:88;                % 8,16,...,88
% Bands: [8–16], [16–24], ..., [80–88]
freq_bands_M1 = [edges(1:end-1); edges(2:end)];
Nsb1 = size(freq_bands_M1, 2);
n1 = (1:Nsb1)';
w_M1 = exp(-n1 / (2*sigma^2)); % weights for M1

% ----- M2: bands aligned to harmonics -----
% Use min and max stimulus frequencies and harmonics 1..Nh.
fmin = min(fstim);
fmax = max(fstim);
margin = 2;  % small margin around harmonic ranges

freq_bands_M2 = zeros(2, Nh);
for h = 1:Nh
    low  = max(8,  h * fmin - margin);
    high = min(88, h * fmax + margin);
    freq_bands_M2(:, h) = [low; high];
end
Nsb2 = size(freq_bands_M2, 2);
n2 = (1:Nsb2)';
w_M2 = exp(-n2 / (2*sigma^2)); % weights for M2

% ----- M3: cumulative bands (Chen-style) -----
% Increasing low cut, fixed high cut at 88 Hz
lows  = 8:8:80;         % e.g., 8,16,...,80
highs = 88 * ones(size(lows));
freq_bands_M3 = [lows; highs];
Nsb3 = size(freq_bands_M3, 2);
n3 = (1:Nsb3)';
w_M3 = exp(-n3 / (2*sigma^2)); % weights for M3

% ----- Pack into a struct for easy looping -----
fbSchemes(1).name = 'M1 equal 8 Hz';
fbSchemes(1).freq_bands = freq_bands_M1;
fbSchemes(1).w = w_M1;

fbSchemes(2).name = 'M2 harmonic-aligned';
fbSchemes(2).freq_bands = freq_bands_M2;
fbSchemes(2).w = w_M2;

fbSchemes(3).name = 'M3 cumulative';
fbSchemes(3).freq_bands = freq_bands_M3;
fbSchemes(3).w = w_M3;

nSchemes = numel(fbSchemes);


%% Construct sine-cosine reference signal for each stimulus according to equation 2
Xref = mySinCosReference(fstim,duration,Nh,Fs);


%%
%% SSVEP frequency detection using filter bank CCA (FBCCA)
Accuracy = zeros(35, nSchemes);   % 35 subjects x 3 filter-bank schemes

for sbj = 1:35
    load(['dataset/SSVEP-BCI-Data/S', num2str(sbj), '.mat']);
    EEGdata= cat(3, data(:,:,:,1), data(:,:,:,2), data(:,:,:,3), ...
                     data(:,:,:,4), data(:,:,:,5), data(:,:,:,6));
    clear data

    nTrials = size(EEGdata, 3);

    % Loop over filter-bank schemes
    for m = 1:nSchemes
        freq_bands = fbSchemes(m).freq_bands;
        w          = fbSchemes(m).w;

        y_pred = zeros(1, nTrials);

        for i = 1:nTrials
            X = EEGdata(:, position, i)';   % [time x channels]
            % apply designed band-pass filter [8–90 Hz]
            X = filtfilt(b,a,X);

            % compute canonical correlation in each sub-band
            nSubBands = size(freq_bands, 2);
            Rho_sb = zeros(nSubBands, size(Xref,3)); % subbands x freqs

            for sb = 1:nSubBands
                % filter data using one of the frequency bands
                [b2,a2] = butter(3, [freq_bands(:,sb)]/(Fs/2), 'bandpass');
                X_sb = filtfilt(b2,a2, X(:, chn));   % [time x nCh]

                for k = 1:size(Xref,3)
                    [~,~,r] = myCCA(X_sb', Xref(:,:,k));
                    Rho_sb(sb,k) = r(1);    % first canonical correlation
                end
            end

            % combine correlation coefficients according to Eq. (6)
            W = repmat(w, 1, size(Rho_sb,2));   % same weights for all freqs
            Rho = sum(W .* (Rho_sb.^2), 1);     % 1 x nFreqs

            % determine stimulus frequency
            [~, ind] = max(Rho);
            y_pred(i) = ind;
        end

        % Performance evaluation for this subject + scheme
        C = confusionmat(y_true, y_pred);
        acc = sum(diag(C)) / sum(C(:));  % 0–1

        Accuracy(sbj, m) = acc * 100;
        fprintf('Subj %2d, %s: Acc = %.2f %%\n', ...
            sbj, fbSchemes(m).name, Accuracy(sbj,m));
    end
end

% Average across subjects
plusminu = char(177);
meanAcc = mean(Accuracy, 1); 
stdAcc  = std(Accuracy, 0, 1) / sqrt(35);

for m = 1:nSchemes
    fprintf('Average accuracy %s: %.2f %s %.2f %%\n', ...
        fbSchemes(m).name, meanAcc(m), plusminu, stdAcc(m));
end

%% Optional: simple bar plot
x = 1:nSchemes;

figure;
errorbar(x, meanAcc, stdAcc, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlim([0.5 nSchemes+0.5]);
xticks(x);
xticklabels({fbSchemes.name});
xtickangle(20);
ylabel('Accuracy (%)');
title('FBCCA accuracy for different filter-bank schemes');
grid on;










%% SSVEP frequency detection using filter bank CCA (FBCCA)
for sbj=1:35
    load(['dataset/SSVEP-BCI-Data/S',num2str(sbj),'.mat/','S',num2str(sbj),'.mat'])
    EEGdata= cat(3,data(:,:,:,1),data(:,:,:,2),data(:,:,:,3),data(:,:,:,4),...
        data(:,:,:,5),data(:,:,:,6));
    clear data
    % frequency recognition
    y_pred= zeros(1,size(EEGdata,3));
    for i=1:size(EEGdata,3)
        X= EEGdata(:,position,i)';
        % apply designed band-pass filter[8-90Hz]
        X= filtfilt(b,a,X);
        % calculate cannonical correlation between the EEG sub-bands and each of the reference signals(Xref)
        for sb=1:size(freq_bands,2)
            % filter data using one of the frequency bands
            [b2,a2]= butter(3,[freq_bands(:,sb)]/(Fs/2), 'bandpass');
            X_sb= filtfilt(b2,a2,X(:,chn));
            for k= 1:size(Xref,3)
                [~,~,temp(:,k)] = myCCA(X_sb',Xref(:,:,k));
            end
            Rho_sb(sb,:)= max(temp);
            temp=[];
        end
        % combine correlation coefficients according to equation(6)
        W= repmat(w,1,size(Rho_sb,2));
        Rho= sum(W.* (Rho_sb.^2));
        % determine the the stimulus frequency of EEG signal(X)
        [mx,ind]= max(Rho);
        y_pred(i)= ind;
        Rho_sb=[];
    end
    %% Performance evaluation
    C= confusionmat(y_true,y_pred); %cunfusion matrix
    Accuracy(sbj)= sum(diag(C)) / sum(C(:)) *100; % accuracy
    disp(['Accuracy(',num2str(sbj),'): ', num2str(Accuracy(sbj)),' %'])
end
plusminu=char(177);
stderror= std( Accuracy ) / sqrt( length( Accuracy ));
tderror= std( Accuracy ) / sqrt( length( Accuracy ));
Ave_Acc_across_sbjs= mean(Accuracy );
disp(['Average accuracy: ',num2str(mean(Accuracy))...
    ,' ',plusminu,' ',num2str(stderror),' %'])


