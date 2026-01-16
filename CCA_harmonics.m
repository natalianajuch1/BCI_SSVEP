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
time= linspace(0,6,1500);

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
%Nh=5;
%Xref = mySinCosReference(fstim,duration,Nh,Fs);
time_list = [0.5 1 1.5 2 2.5 3];
Nh_list = 1:10;
Acc_time = zeros(length(Nh_list), length(time_list), 35);

%% SSVEP frequency detection using Standard CCA
% For each Nh, rebuild reference signals
for h = 1:length(Nh_list)
    Nh = Nh_list(h); % update Nh for the current iteration
    disp(['Processing Nh = ', num2str(Nh)])

    for t = 1:length(time_list)
        duration = time_list(t);
        disp(['  Time window = ', num2str(duration), ' s'])
        position= find(time>=0.5 & time<=0.5+duration); % index of EEG signal

        Xref = mySinCosReference(fstim, duration, Nh, Fs); % rebuild reference signals
    
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

                temp = zeros(1, size(Xref,3));

                % calculate cannonical correlation between the EEG signal(X) and each of the reference signals(Xref)
                for j= 1:size(Xref,3)
                    [~,~,r] = myCCA(X(:,chn)',Xref(:,:,j));
                    temp(j) = max(r);
                end
                [~, y_pred(i)] = max(temp);
            end
            %% Performance evaluation
            C= confusionmat(y_true,y_pred); %cunfusion matrix
            Acc_vs_Nh(h,t,sbj) = sum(diag(C)) / sum(C(:)) * 100;
            disp(['Accuracy = ', num2str(Acc_vs_Nh(h,t,sbj)), ' %'])
        end
    end
end

meanAcc = mean(Acc_vs_Nh, 3); 

figure; hold on;
for h = 1:length(Nh_list)
    plot(time_list, meanAcc(h,:), '-o', 'LineWidth', 2, ...
        'DisplayName', ['N_h = ', num2str(Nh_list(h))]);
end

xlabel('Time window length (s)');
ylabel('Accuracy (%)');
title('Accuracy vs Time for Different Numbers of Harmonics (CCA)');
legend('Location','southeast');
grid on;
