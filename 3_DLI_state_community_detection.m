clear
clc

datapath = '/share/inspurStorage/home1/zhangruohan/Documents/LanguageComprehension';
load(fullfile(datapath,'results','dLI_values.mat'));
load(fullfile(datapath,'results','dLI_temporal_clustering.mat'));

% spatial correlation of dynamic laterality index
K = max(GroupStateCentroid.ind_GroupState);
dLI_allsub = [];
dLI_stateseries_allsub = [];
for i = 1:length(dLI_results)
    dLI_allsub = [dLI_allsub reshape(dLI_results(i).dLI, size(dLI_results(i).dLI,1),[],1)];
    dLI_stateseries_allsub = [dLI_stateseries_allsub; reshape(RelabelStateSeries(i).NewStateSeries,[],1)];
end

mean_dLI_spatialcorr_allsub.allstate = corr(dLI_allsub');
for k = 1:K
    str = ['state' num2str(k)];
    mean_dLI_spatialcorr_allsub.(str) = corr(dLI_allsub(:,dLI_stateseries_allsub==k)');
end

% spatial communitiy detection
str = {'allstate','state1','state2'};
prop_thre = 0.3;
for i = 1:length(str)
    dLI_SpatialCorr(i).state = str{i};
    W = mean_dLI_spatialcorr_allsub.(str{i});
    dLI_SpatialCorr(i).spatial_corr = W;
    W(logical(eye(size(W,1)))) = 0;
    
    n  = size(W,1);             % number of nodes
    M  = 1:n;                   % initial community affiliations
    Q0 = -1; 
    Q1 = 0;            % initialize modularity values
    tic;
    while Q1 - Q0 > 1e-100           % while modularity increases
        Q0 = Q1;                % perform community detection
        [M, Q1] = community_louvain(W .* (W > 0),[],M);
    end
    toc;
    
    dLI_SpatialCorr(i).community_ind = M;
    dLI_SpatialCorr(i).community_num = max(M);

end

savefile = fullfile(datapath,'results','dLI_StateSpatialCommunity.mat');
save(savefile,'dLI_SpatialCorr');



