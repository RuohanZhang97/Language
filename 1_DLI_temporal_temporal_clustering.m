clear
clc

datapath = '/share/inspurStorage/home1/zhangruohan/Documents/LanguageComprehension';
load(fullfile(datapath,'results','dLI_values.mat'));
NumSub = 20; % number of subjects
NumNode = 426; % brain regions
WinLen = 20; % set length of window
interval = 1; % sliding step 1 TR

NumWin_AllArticle = [];
ind_NumWin_AllArticle = [];
for i = 1:length(dLI_results)
    NumWin_AllArticle = [NumWin_AllArticle dLI_results(i).NumWin];
    ind_NumWin_AllArticle = [ind_NumWin_AllArticle repelem(i,dLI_results(i).NumWin)];
end

%% Step 1: temporal clustering at individual level
for i = 1:NumSub
    disp(num2str(i))
    dLI_state_results(i).NoSub = i;
    sub_dLI = [];
    for j = 1:length(dLI_results)
        sub_dLI = [sub_dLI dLI_results(j).dLI(:,:,i)];
    end
    dLI_state_results(i).sub_dLI = sub_dLI;
    W = corr(sub_dLI);
    W(logical(eye(size(W,1)))) = 0;
    n  = size(W,1);             % number of nodes
    M  = 1:n;                   % initial community affiliations
    Q0 = -1; 
    Q1 = 0;            % initialize modularity values
    tic;
    while Q1-Q0>1e-100           % while modularity increases
        Q0 = Q1;                % perform community detection
        [M, Q1] = community_louvain(W,[],M,'negative_sym');
    end
    toc;
    
    sub_dLI_state_centroid = struct();
    for k = 1:max(M)
        eval(['sub_dLI_state_centroid.state' num2str(k) ' = mean(sub_dLI(:,M==k),2);'])
    end
    
    dLI_state_results(i).SubStateSeries = M;
    dLI_state_results(i).SubStateCentroid = sub_dLI_state_centroid;
    dLI_state_results(i).NumNode = NumNode;
    dLI_state_results(i).WinLen = WinLen;
    dLI_state_results(i).interval = interval;
end

%% Step 2: temporal clustering at group level
NumSubStateCentroid = 0;
ind_NoSubStateCentroid = [];
AllSubStateCentroid = [];
for i = 1:length(dLI_state_results)
    state_names = fieldnames(dLI_state_results(i).SubStateCentroid);
    NumSubStateCentroid = NumSubStateCentroid + length(state_names);
    ind_NoSubStateCentroid = [ind_NoSubStateCentroid; i * ones(length(state_names),1)];
    for j = 1:length(state_names)
        M = eval(['dLI_state_results(i).SubStateCentroid.' state_names{j}]);
        AllSubStateCentroid = [AllSubStateCentroid M];
    end
end

W = corr(AllSubStateCentroid);
W(logical(eye(size(W,1)))) = 0;
n  = size(W,1);             % number of nodes
M  = 1:n;                   % initial community affiliations
Q0 = -1; 
Q1 = 0;            % initialize modularity values
tic;
while Q1-Q0>1e-100           % while modularity increases
    Q0 = Q1;                % perform community detection
    [M, Q1] = community_louvain(W,[],M,'negative_sym');
end
toc;
GroupStateCentroid.ind_GroupState = M;
GroupStateCentroid.StateCentroid = [ind_NoSubStateCentroid GroupStateCentroid.ind_GroupState];

AllSubNewStateSeries = [];
for i = 1:length(dLI_state_results)
    relabel = GroupStateCentroid.StateCentroid(GroupStateCentroid.StateCentroid(:,1) == i, 2);
    sub_new_state_series = zeros(size(dLI_state_results(i).SubStateSeries));
    for j = 1:length(relabel)
        sub_new_state_series(dLI_state_results(i).SubStateSeries == j) = relabel(j);
    end
    AllSubNewStateSeries = [AllSubNewStateSeries; sub_new_state_series'];
end

RelabelStateSeries = struct();
for i = 1:length(NumWin_AllArticle)
    RelabelStateSeries(i).NoArticle = i;
    RelabelStateSeries(i).NewStateSeries = AllSubNewStateSeries(:,ind_NumWin_AllArticle==i)'; 
    RelabelStateSeries(i).NumWin = NumWin_AllArticle(1,i);
    RelabelStateSeries(i).WinLen = WinLen;
    RelabelStateSeries(i).interval  = interval ;
    RelabelStateSeries(i).NumNode = NumNode;
end
    
savefile = fullfile(datapath,'results','dLI_temporal_clustering.mat');
save(savefile,'RelabelStateSeries','GroupStateCentroid');

