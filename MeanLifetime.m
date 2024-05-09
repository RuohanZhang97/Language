function result = MeanLifetime(state_order,occurrence)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    state = unique(state_order);
    result = struct();
    for k = 1:length(state)
        result(k).state = state(k);
        result(k).value = mean(occurrence(state_order == state(k)));
    end
end

