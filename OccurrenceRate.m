function result = OccurrenceRate(state_order,occurrence)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    state = unique(state_order);
    result = struct();
    for k = 1:length(state)
        result(k).state = state(k);
        result(k).value = sum(occurrence(state_order == state(k)))/sum(occurrence);
    end
end

