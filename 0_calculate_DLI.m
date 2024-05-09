clear
clc

datapath = '/share/inspurStorage/home1/zhangruohan/Documents/LanguageComprehension';
language_data = dir(fullfile(datapath,'sorted_modified_HCPex_data'));
language_data(1:2) = [];
load(fullfile(datapath,'data_supply','TimingInfo.mat'));

NumNode = 426; % brain regions
WinLen = 20; % set length of window
interval = 1; % sliding step 1 TR

for i = 1:length(language_data)
    disp(num2str(i))
    load(fullfile(language_data(i).folder,language_data(i).name));
    NumSub = size(atlasdata,2);
    for j = 1:NumSub
        atlasdata{j} = atlasdata{j}(:,StartTimePoints:TimingInfo(i).new_EndTimePoints);
    end
    tc = size(atlasdata{1},2);
    NumWin = (tc - WinLen)/interval + 1; % number of window
    dLI = zeros(NumNode,NumWin,NumSub);
    for s = 1:NumSub
        for t = 1:NumWin
            GS_L = mean(atlasdata{s}([1:180 361:393],t:t+WinLen-1));
            GS_R = mean(atlasdata{s}([181:360 394:426],t:t+WinLen-1));
            dLI(:,t,s) = fisherz(corr(GS_L',atlasdata{s}(:,t:t+WinLen-1)'))' - fisherz(corr(GS_R',atlasdata{s}(:,t:t+WinLen-1)'))';    
        end     
    end
    
    dLI_results(i).dLI = dLI;
    dLI_results(i).MLI = squeeze(mean(dLI,2));   % mean laterality index
    dLI_results(i).LF = squeeze(std(dLI,[],2));   % laterality fluctuations
    dLI_results(i).NumNode = NumNode;
    dLI_results(i).tc = tc;
    dLI_results(i).WinLen = WinLen;
    dLI_results(i).interval = interval;
    dLI_results(i).NumWin = NumWin;
end

savefile = fullfile(datapath,'results','dLI_values.mat');
save(savefile,'dLI_results');

