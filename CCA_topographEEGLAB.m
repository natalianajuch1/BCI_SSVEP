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
% Load EEGLAB channel locations
chanlocs = readlocs('Standard-10-20-Cap81.ced');

% SSVEP channels
target_labels = {'O2','Oz','O1','PO8','PO4','POz','PO3','PO7'};

% Find indices in the layout file (for plotting)
chn_layout = find(ismember({chanlocs.labels}, target_labels));

chn_dataset=[52 53 55 56 57 58 61 62];
% Select only used channels
%x = x_all(chn);
%y = y_all(chn);

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

    TopoData = zeros(length(chn_dataset), size(EEGdata,3));
    % frequency recognition
    for i=1:size(EEGdata,3)
        X= EEGdata(:,position,i)'; % EEG signal
        % apply designed band-pass filter[8-90Hz]
        X= filtfilt(b,a,X);
        % calculate cannonical correlation between the EEG signal(X) and each of the reference signals(Xref)
        for j= 1:size(Xref,3)
            for c = 1:length(chn_dataset)
                [~,~,r] = myCCA(X(:,chn_dataset(c))',Xref(:,:,j));
                temp_chan(c,j) = max(r);
            end
        end
        % calculates the maximum correlation between the EEG signal(X) and each of the reference signals(Xref)
        Rho_chan= max(temp_chan,[],2);
        TopoData(:,i) = Rho_chan;
        % determine the the stimulus frequency of EEG signal(X)
        Rho = max(temp_chan,[],1);
        [~,ind]= max(Rho);
        y_pred(i)= ind;
    end
    TopoMean(:,sbj) = mean(TopoData,2);

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

%% Plot Topological Map
maxTopo = max(TopoMean(:));

nRows = 5;
nCols = 7;

for sbj = 1:35
    % Create a vector with all electrodes initialized to zero
    vals_full = zeros(length(chanlocs),1);
    
    % Assign measured occipital/parietal channel values
    vals_full(chn_layout) = TopoMean(:,sbj);
    
    % Plot topoplot with fixed color scale
    figure('Name',['Subject ', num2str(sbj)],'NumberTitle','off');
    topoplot(vals_full, chanlocs, 'electrodes','labels','style','both');
    colorbar;
    caxis([0 maxTopo]); % fixed color scale across all subjects
    title(['Subject ', num2str(sbj), ' – CCA Topography'])
end

figure('Name','All Subjects – CCA Topography','NumberTitle','off');
for sbj = 1:35
    % Assign measured occipital/parietal channel values to full electrode vector
    vals_full = zeros(length(chanlocs),1);
    vals_full(chn_layout) = TopoMean(:,sbj);
    
    % Subplot for each subject
    subplot(nRows, nCols, sbj);
    topoplot(vals_full, chanlocs, 'electrodes','on','style','both');
    caxis([0 max(TopoMean(:))]);           % fixed color scale across all subjects
    title(['S', num2str(sbj)], 'FontSize', 8);
end

% Add a single colorbar for the entire figure
h = colorbar('Position',[0.92 0.1 0.02 0.8]); % adjust position
ylabel(h,'CCA Contribution','FontSize',10);

%% Plot accuracy across subjects
figure;
bar(Accuracy);
xlabel('Subject');
ylabel('Accuracy (%)');
title('CCA Accuracy Across Subjects (Single channel)');
grid on;

% Add mean accuracy line
hold on;
yline(mean(Accuracy),'r--','LineWidth',2);

% Print mean value on the plot
text(1, mean(Accuracy)+1, ...
    ['Mean = ', num2str(mean(Accuracy),'%.2f'), '%'], ...
    'Color','r','FontWeight','bold');

legend('Subject Accuracy','Mean Accuracy');