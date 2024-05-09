clear
clc

datapath = '/home1/zhangruohan/Documents/LanguageComprehension/';
load(fullfile(datapath,'results','dLI_temporal_clustering.mat'));

StateSeries = struct();
OccurrenceRateResults = struct();
MeanLifetimeResults = struct();
TransitionNumberResults = struct();
for i = 1:length(RelabelStateSeries)
    for j = 1:size(RelabelStateSeries(i).NewStateSeries,2)
        [StateSeries(i).state{j,1},StateSeries(i).state{j,2}] = StateOccurrence(RelabelStateSeries(i).NewStateSeries(:,j));
        
        sub_OR_result(j).NoSub = j;
        sub_OR_result(j).occurrence_rate = OccurrenceRate(StateSeries(i).state{j,1},StateSeries(i).state{j,2});        
        OccurrenceRateResults(i).Result = sub_OR_result;
        OccurrenceRateResults(i).NumWin = RelabelStateSeries(i).NumWin;
        
        sub_ML_result(j).NoSub = j;
        sub_ML_result(j).mean_lifetime = MeanLifetime(StateSeries(i).state{j,1},StateSeries(i).state{j,2});        
        MeanLifetimeResults(i).Result = sub_ML_result;
        MeanLifetimeResults(i).NumWin = RelabelStateSeries(i).NumWin;
        
        sub_TN_result(j).NoSub = j;
        sub_TN_result(j).transition_number = TransitionNumber(StateSeries(i).state{j,1});        
        TransitionNumberResults(i).Result = sub_TN_result;
        TransitionNumberResults(i).NumWin = RelabelStateSeries(i).NumWin;
    end
end

savefile = fullfile(datapath,'results','dLI_temporal_properties.mat');
save(savefile,'StateSeries','OccurrenceRateResults','MeanLifetimeResults','TransitionNumberResults');

