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
duration=2.5;
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
%% Construct sine-cosine reference signal for each stimulus according to equation 2
% number of harmonics
Nh_list = 1:10;
Acc_vs_Nh = zeros(length(Nh_list), 35);  % harmonics × subjects

for h = 1:length(Nh_list)
    Nh = Nh_list(h);
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
            X= EEGdata(:,position,i)'; % EEG signal
            % apply designed band-pass filter[8-90Hz]
            X= filtfilt(b,a,X);
            % calculate cannonical correlation between the EEG signal(X) and each of the reference signals(Xref)
            for j= 1:size(Xref,3)
                [~,~,temp(:,j)] = myCCA(X(:,chn)',Xref(:,:,j));
            end
            % calculates the maximum correlation between the EEG signal(X) and each of the reference signals(Xref)
            Rho= max(temp);
            % determine the the stimulus frequency of EEG signal(X)
            [mx,ind]= max(Rho);
            y_pred(i)= ind;
        end
        %% Performance evaluation
        C= confusionmat(y_true,y_pred); %cunfusion matrix
        Accuracy(sbj)= sum(diag(C)) / sum(C(:)) *100; % accuracy
        Acc_vs_Nh(h, sbj) = Accuracy(sbj);
        disp(['Accuracy(',num2str(sbj),'): ', num2str(Accuracy(sbj)),' %'])
    end
end
plusminu=char(177);
stderror= std( Accuracy ) / sqrt( length( Accuracy ));
tderror= std( Accuracy ) / sqrt( length( Accuracy ));
Ave_Acc_across_sbjs= mean(Accuracy );
disp(['Average accuracy: ',num2str(mean(Accuracy))...
    ,' ',plusminu,' ',num2str(stderror),' %'])

meanAcc = mean(Acc_vs_Nh, 2);
stdErr  = std(Acc_vs_Nh, [], 2) / sqrt(35);

errorbar(Nh_list, meanAcc, stdErr, '-o','LineWidth',2);
xlabel('Number of harmonics N_h');
ylabel('Accuracy (%)');
grid on;
title('Effect of Number of Harmonics on CCA Performance');
