clear
clc

datapath = '/home1/zhangruohan/Documents/LanguageComprehension/';
load(fullfile(datapath,'results','dLI_temporal_properties.mat'));

%% sort OccurrenceRateResults for each state
for i = 1:length(OccurrenceRateResults)
    for j = 1:length(OccurrenceRateResults(i).Result)
        State1_OccurrenceRateResults(i).Result(j).NoSub = OccurrenceRateResults(i).Result(j).NoSub;
        ind_state = find([OccurrenceRateResults(i).Result(j).occurrence_rate.state] == 1);
        if ~isempty(ind_state)
            occurrence_rate_value = [OccurrenceRateResults(i).Result(j).occurrence_rate.value];
            State1_OccurrenceRateResults(i).Result(j).occurrence_rate = occurrence_rate_value(ind_state);
        else
            State1_OccurrenceRateResults(i).Result(j).occurrence_rate = 0;
        end
    end
    State1_OccurrenceRateResults(i).NumWin = OccurrenceRateResults(i).NumWin;
end

for i = 1:length(OccurrenceRateResults)
    for j = 1:length(OccurrenceRateResults(i).Result)
        State2_OccurrenceRateResults(i).Result(j).NoSub = OccurrenceRateResults(i).Result(j).NoSub;
        ind_state = find([OccurrenceRateResults(i).Result(j).occurrence_rate.state] == 2);
        if ~isempty(ind_state)
            occurrence_rate_value = [OccurrenceRateResults(i).Result(j).occurrence_rate.value];
            State2_OccurrenceRateResults(i).Result(j).occurrence_rate = occurrence_rate_value(ind_state);
        else
            State2_OccurrenceRateResults(i).Result(j).occurrence_rate = 0;
        end
    end
    State2_OccurrenceRateResults(i).NumWin = OccurrenceRateResults(i).NumWin;
end

%% sort MeanLifetimeResults for each state
for i = 1:length(MeanLifetimeResults)
    for j = 1:length(MeanLifetimeResults(i).Result)
        State1_MeanLifetimeResults(i).Result(j).NoSub = MeanLifetimeResults(i).Result(j).NoSub;
        ind_state = find([MeanLifetimeResults(i).Result(j).mean_lifetime.state] == 1);
        if ~isempty(ind_state)
            mean_lifetime_value = [MeanLifetimeResults(i).Result(j).mean_lifetime.value];
            State1_MeanLifetimeResults(i).Result(j).mean_lifetime = mean_lifetime_value(ind_state);
        else
            State1_MeanLifetimeResults(i).Result(j).mean_lifetime = 0;
        end
    end
    State1_MeanLifetimeResults(i).NumWin = MeanLifetimeResults(i).NumWin;
end

for i = 1:length(MeanLifetimeResults)
    for j = 1:length(MeanLifetimeResults(i).Result)
        State2_MeanLifetimeResults(i).Result(j).NoSub = MeanLifetimeResults(i).Result(j).NoSub;
        ind_state = find([MeanLifetimeResults(i).Result(j).mean_lifetime.state] == 2);
        if ~isempty(ind_state)
            mean_lifetime_value = [MeanLifetimeResults(i).Result(j).mean_lifetime.value];
            State2_MeanLifetimeResults(i).Result(j).mean_lifetime = mean_lifetime_value(ind_state);
        else
            State2_MeanLifetimeResults(i).Result(j).mean_lifetime = 0;
        end
    end
    State2_MeanLifetimeResults(i).NumWin = MeanLifetimeResults(i).NumWin;
end

savefile = fullfile(datapath,'results','sorted_dLI_temporal_properties.mat');
save(savefile,'StateSeries','State1_OccurrenceRateResults','State2_OccurrenceRateResults', ...
    'State1_MeanLifetimeResults','State2_MeanLifetimeResults', 'TransitionNumberResults');

