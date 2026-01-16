clc;
clear
close all;

%% define prameters 
Fs=250;% sampling rate

% Nine channels that are used for analysis
% [O2, Oz, O1, PO6, PO4, POZ, PO3, PO7, and P8]
%chn =[52 53 55 56 57 58 61 62 63];

% data lenght in seconds 
duration = 1.5;
time = linspace(0,6,1500);
position = find(time>=0.5 & time<=0.5+duration); % index of EEG signal

% design a band-pass butterworth filter
[b,a] = butter(3,[8 90]/(Fs/2), 'bandpass');


% load frequency-phase information of stimuli
load('dataset\SSVEP-BCI-Data\Freq_Phase.mat')
fstim= freqs;
% build label for each stimulus which will be used for evaluatoin
% freqs= repmat(freqs,1,6);
%y_true= repmat(1:40,1,6);

%% Construct sine-cosine reference signal for each stimulus according to equation 2
% number of harmonics
Nh=5;
Xref = mySinCosReference(fstim,duration,Nh,Fs);


%%
% Test for the number of electrodes

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

% -------------------------------------------------------------------------
% Electrode sets for factor "number of electrodes"
% -------------------------------------------------------------------------

% 1 channel: Oz
set_1ch  = [idx_Oz];

% 3 channels: O1, Oz, O2
set_3ch  = [idx_O1 idx_Oz idx_O2];

% 5 channels: Pz, POz, O1, Oz, O2
set_5ch = [idx_Pz idx_POz idx_O1 idx_Oz idx_O2];

% 7 channels: Pz, PO3, POz, PO4, O1, Oz, O2
set_7ch = [idx_Pz idx_PO3 idx_POz idx_PO4 idx_O1 idx_Oz idx_O2];

% 9 channels: Chen 2015 set
set_9ch = [idx_Pz idx_PO5 idx_PO3 idx_POz idx_PO4 idx_PO6 idx_O1 idx_Oz idx_O2];

electrodeSets = {set_1ch, set_3ch, set_5ch, set_7ch, set_9ch};
setLabels     = {'1ch (Oz)', '3ch (O1,Oz,O2)', '5ch', '7ch', '9ch (Chen)'};
nSets         = numel(electrodeSets);


y_true = repmat(1:40,1,6);

Accuracy = zeros(35, nSets); % 35 subjects x 4 electrode sets
ITR      = zeros(35, nSets); % same size

trialLengthSec = duration;   % 2.5 s, adjust if you use a different effective trial time
breakSec       = 0;          % set non-zero if you want to include inter-trial interval
M = numel(fstim);            % number of targets (40)

for sbj = 1:35
    load(['dataset/SSVEP-BCI-Data/S',num2str(sbj),'.mat'])
    
    % concatenate all trials to construct a 3D matrix
    EEGdata = cat(3, data(:,:,:,1), data(:,:,:,2), data(:,:,:,3), ...
                     data(:,:,:,4), data(:,:,:,5), data(:,:,:,6));
    clear data

    nTrials = size(EEGdata, 3);

    % loop over electrode sets
    for s = 1:nSets
        chn = electrodeSets{s};
        y_pred = zeros(1, nTrials);

        for i = 1:nTrials
            X = EEGdata(:, position, i)';  % [time x channels]
            % band-pass filter [8-90 Hz]
            X = filtfilt(b, a, X);

            % compute canonical correlation for each stimulus frequency
            for j = 1:size(Xref,3)
                [~, ~, temp(:,j)] = myCCA(X(:, chn)', Xref(:,:,j));
            end
            Rho = max(temp);            % first canonical correlation per freq
            [~, ind] = max(Rho);        % pick the max
            y_pred(i) = ind;
        end

        % Performance evaluation for this subject + electrode set
        C = confusionmat(y_true, y_pred);
        acc = sum(diag(C)) / sum(C(:));     % 0..1
        Accuracy(sbj, s) = acc * 100;       % store as percentage

        % ITR (bits/min)
        ITR(sbj, s) = computeITR(M, acc, trialLengthSec, breakSec);

        fprintf('Subj %2d, %s: Acc = %.2f %% , ITR = %.2f bits/min\n', ...
                sbj, setLabels{s}, Accuracy(sbj,s), ITR(sbj,s));
    end
end

% Average across subjects
plusminu = char(177);

meanAcc = mean(Accuracy, 1);               % 1 x nSets
stdAcc  = std(Accuracy, 0, 1) / sqrt(35);  % standard error

meanITR = mean(ITR, 1);
stdITR  = std(ITR, 0, 1) / sqrt(35);

for s = 1:nSets
    fprintf('Average accuracy %s: %.2f %s %.2f %%\n', ...
        setLabels{s}, meanAcc(s), plusminu, stdAcc(s));
    fprintf('Average ITR %s: %.2f %s %.2f bits/min\n', ...
        setLabels{s}, meanITR(s), plusminu, stdITR(s));
end

figure;
subplot(1,2,1);
bar(meanAcc);
set(gca, 'XTickLabel', setLabels, 'XTickLabelRotation', 20);
ylabel('Accuracy (%)');
title('Accuracy vs number of electrodes');

subplot(1,2,2);
bar(meanITR);
set(gca, 'XTickLabel', setLabels, 'XTickLabelRotation', 20);
ylabel('ITR (bits/min)');
title('ITR vs number of electrodes');
