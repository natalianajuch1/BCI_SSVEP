clc;
clear
close all;
%% Standard CCA in SSVEP frequency detection (demo code)
% by    Mohammad Norizadeh Cherloo,
%       Homa Kashefi Amiri,
%       Amir Mohammad Mijani,
%       Liang Zhan,
%       Mohammad Reza Daliri

%% define prameters (Fs,data length, channels and Nh,...)
Fs=250;% sampling rate

% Nine channels that are used for analysis
% [O2, Oz, O1, PO6, PO4, POZ, PO3, PO7, and P8]
chn=[52 53 55 56 57 58 61 62 63];
% chn=1:64;
% data lenght in seconds (0.5,1,1.5,2,2.5 and 3 were considered in our study)
duration = 2.5;
time= linspace(0,6,1500);
position= find(time>=0.5 & time<=0.5+duration); % index of EEG signal

% design a band-pass butterworth filter
[b,a]= butter(3,[8 90]/(Fs/2), 'bandpass');

% load frequency-phase information of stimuli
load('C:\Users\marty\Desktop\BCI\Project\SSVEP-BCI-Data\Freq_Phase.mat')
fstim= freqs;
% build label for each stimulus which will be used for evaluatoin
% freqs= repmat(freqs,1,6);
y_true= repmat(1:40,1,6);

all_correct_rho   = [];
all_incorrect_rho = [];

%% Construct sine-cosine reference signal for each stimulus according to equation 2
% number of harmonics
Nh=5;
Xref = mySinCosReference(fstim,duration,Nh,Fs);

%% SSVEP frequency detection using Standard CCA
for sbj=1:35
    load(['C:\Users\marty\Desktop\BCI\Project\SSVEP-BCI-Data\S\S',num2str(sbj),'.mat'])
    % concatenate all trials to costruct a 3 dimension matrix
    EEGdata= cat(3,data(:,:,:,1),data(:,:,:,2),data(:,:,:,3),data(:,:,:,4),...
        data(:,:,:,5),data(:,:,:,6));
    clear data

    y_pred= zeros(1,size(EEGdata,3));
    % frequency recognition
    for i=1:size(EEGdata,3)
        X= EEGdata(:,position,i)'; % Extracts post-stimulus EEG signal
        % apply designed band-pass filter[8-90Hz] (Removes low-frequency drift, preserving SSVEP harmonics)
        X= filtfilt(b,a,X);
        % calculate cannonical correlation between the EEG signal(X) and each of the reference signals(Xref)
        for j= 1:size(Xref,3)
            [Wx,~,r] = myCCA(X(:,chn)', Xref(:,:,j));
            temp(:,j) = r;
            if sbj == 1 && i == 1
                Wx_all{j} = Wx(:,1);   % store first canonical component
            end
        end
        % calculates the maximum correlation between the EEG signal(X) and each of the reference signals(Xref) 
        % For each trial, Rho = [r1 r2 ... r40] is the feature vector
        Rho= max(temp);
        % determine the stimulus frequency of EEG signal(X)
        [mx,ind]= max(Rho);
        y_pred(i)= ind;

        % Store Rho values (maximum correlation value)
        if y_pred(i) == y_true(i)
            all_correct_rho   = [all_correct_rho, mx];
        else
            all_incorrect_rho = [all_incorrect_rho, mx];
        end
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

%% Plot confusion matrix
figure;
imagesc(C);
colorbar;
xlabel('Predicted Frequency');
ylabel('True Frequency');
title(['Confusion Matrix – Subject ', num2str(sbj)]);

%% Plot accuracy across subjects
figure;
bar(Accuracy);
xlabel('Subject');
ylabel('Accuracy (%)');
title('CCA Accuracy Across Subjects');
grid on;

% Add mean accuracy line
hold on;
yline(mean(Accuracy),'r--','LineWidth',2);

% Print mean value on the plot
text(1, mean(Accuracy)+1, ...
    ['Mean = ', num2str(mean(Accuracy),'%.2f'), '%'], ...
    'Color','r','FontWeight','bold');

legend('Subject Accuracy','Mean Accuracy');

%% Correlation score distribution (Correct vs Incorrect)
% This shows why errors occur and the separation power of CCA features
data = [all_correct_rho, all_incorrect_rho];
group = [ones(1,length(all_correct_rho)), 2*ones(1,length(all_incorrect_rho))];

figure;
boxplot(data, group, 'Labels', {'Correct','Incorrect'});
ylabel('Max CCA Correlation');
title('CCA Correlation Distribution');
grid on;

%% Plot correlation scores
figure;
plot(Rho);
xlabel('Stimulus index');
ylabel('Canonical correlation');
title('CCA correlation scores');

% Load electrode positions
chanlocs = readlocs('C:\Users\marty\Downloads\eeglab2025.1.0\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp'); % EEGLAB comes with this file

%% Plot spatial filter
figure;
topoplot(Wx_all{1}, chanlocs(chn));
title('CCA Spatial Filter (1st Component)');
colorbar;
