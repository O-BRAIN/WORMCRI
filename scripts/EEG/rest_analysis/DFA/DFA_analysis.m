clear all
close all
%%
scriptdir = '/data/p_02191/Analysis/Nadine/scripts/matlab/DFA_code'
addpath(scriptdir)
addpath('/data/u_naherzog_software/eeglab/eeglab2022.0/')   
eeglab %eeglab needs to be in the workspace to do topoplot
%define variables

Fs = 500;                   % sampling frequency.
DFA_SmallTime = 0.5; 		% Smallest time window (in seconds) to be computed in DFA.
DFA_LargeTime = 180;        % Largest time window (in seconds) to be computed in DFA.
DFA_SmallTimeFit = 2; 		% Smallest time window (in seconds) to be include in the DFA fit.
DFA_LargeTimeFit = 25;		% Largest time window (in seconds) to be include in the DFA fit.
DFA_Overlap = 0.5;		    % Overlap between windows in DFA.
DFA_Plot = 1

%Alpha CONTROLS
DFAalpha = zeros(62,78);  %electrodes = 61+1, N subjects = 78
test = num2cell(DFAalpha)
hp = 8;   %set frequency pass band, 2-4 for delta, 5-8 for theta, 8-14 for alpha
lp = 14;
%fir_order = 58    	    % Use 58th-order for high time resolution and low frequency resolution.
fir_order = 2 / hp    	%   fir_order - Filterorder in seconds, usually good with 2/hp


%% read in subject list 
datdir = '/path/to/where/clean/data/is/'
cd(datdir)
sublist = dir();
sublist = sublist(~[sublist(:).isdir] & contains({sublist.name}, '.set'));%select only files with .set ending and non directories
%% DFA loop
%subject x electrode loop

for j = 1:length(sublist)
    subid = sublist(j).name
try
    EEG = pop_loadset('filename',subid, 'filepath', datdir);
    Data = EEG.data;
    cd(scriptdir) %go into script directory to find functions below... can also be done by addpath()?
    for i = 1:size(Data,1)
        data = Data(i,:);%data needs to be a vector ... calculate DFA for each electrode
        data_filt = filter_fir(data,hp,lp,Fs,fir_order);
        fprintf('Computing DFA for electrode %i\n',i);
        [DFA_x,DFA_y,DFA_exp] = Scaling_DFA(abs(hilbert(data_filt)),Fs,DFA_SmallTime,DFA_LargeTime,DFA_SmallTimeFit,DFA_LargeTimeFit,DFA_Overlap,DFA_Plot);
        DFAalpha(i+1,j-2) = DFA_exp; %store DFA_exp per electrode per subject in DFAalpha matrix
    end
cd(datdir)%go back to data directory to load new subject file
end
end

%save table
% DFAalpha = cell2table(DFAalpha)
writetable(DFAalpha,"DFAtheta.xlsx")
%% plotting

% t = 1: 104192
% x = Data_filt
% env = abs(hilbert(Data_filt))
% figure 
% plot(t,x)
% hold on
% plot(t,[-1;1]*env,'r','LineWidth',2)
% grid on
% zoom on
% axis([log10(min(DFA_x/Fs))-0.1 log10(max(DFA_x/Fs))+0.1 log10(min(DFA_y(3:end)))-0.1 log10(max(DFA_y))+0.1])
% xlabel('log_{10}(time), [Seconds]','Fontsize',12)


% DFAs
chanlocs = EEG.chanlocs
DFA = zeros(61,1);
pValues = zeros(61,1);
DFA2 = zeros(64,13);
DFAmedPRE = zeros(64,6);
DFAmedPost = zeros (64,6);
DFActrlPRE = zeros (64,8);
DFActrlPost = zeros (64,8);
DFA1 = DFAmedPRE
DFA2 = DFAmedPost
DFA1 = DFActrlPRE
DFA2 = DFActrlPost
 
% DFA_Ctrl = zeros (64,17);


meanDiff = zeros(64,1);
meanDiff = nanmean(DFA2, 2)-nanmean(DFA1, 2);%for topoplot compute difference for each electrode for DFAalphaPre and DFAalphaPost...
chanlocs = EEG.chanlocs  
 pValues = nan(64, 1);
% for channelIdx = 1 : 64
%     % ttest is paired, ttest2 is unpaired
%     
%     [~, pValues(channelIdx)] = ttest(DFA1(channelIdx,:), DFA2(channelIdx,:));
%     tValue(channelIdx) = tinv(pValues(channelIdx),31);
% end
    
figure
    topoplot(avgDFA, chanlocs,'electrodes','on', 'headrad', 'rim', 'maplimits',[-3 3], 'style', 'map', 'numcontour', 2, 'circgrid', 100, 'gridscale', 300, 'shading', 'flat','emarker2', {find(pValues<=0.05),'o','w',1,0.1});
    colorbar('fontsize',18);
     % caxis([-0.01 0.07])
  caxis([min(meanDiff(:)), max(meanDiff(:))])


