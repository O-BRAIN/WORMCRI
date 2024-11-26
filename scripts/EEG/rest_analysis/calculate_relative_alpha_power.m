%calculate relative alpha power per subject
clear all
addpath('/data/u_naherzog_software/eeglab/eeglab2022.0/')   
eeglab

datapath = '/data/p_02191/Analysis/Nadine/EEG/rest/preprocessing/data/'
cd(datapath);
sublist = dir();
sublist = sublist([sublist(:).isdir]);
substart = 3;
subend = 81 ;
sublist(substart).name  %check substart = S002?
sublist(subend).name    %check subend = S098?

vals = {'SubID', 'relAlphaPow'};
count = 2

for sub = substart:subend%length(sublist);
    subid = sublist(sub).name

try
    filepath = [datapath,subid]
    EEG = pop_loadset('filename',[subid, 'post_ICA.set'], 'filepath', [filepath]);
    [PSD, freqs] = pwelch(EEG.data',[],[],[],EEG.srate);
    totalPower = mean(bandpower(PSD, freqs, 'psd'));
    alphaPower = mean(bandpower(PSD, freqs, [8 13],'psd'));

    relativeAlphaPower = alphaPower/totalPower;

    vals{count,1} = subid;
    vals{count,2} = relativeAlphaPower
    count = count + 1;
end
end

vals = cell2table(vals)
writetable(vals,'relativeAlphaPower.xlsx')

    

    
   