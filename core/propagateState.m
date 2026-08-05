function nextState = propagateState(currentState, worldVelocity, params)
    if nargin < 3
        params = defaultConfig();
    end

    nextState = currentState;
    nextState(1:2) = currentState(1:2) + worldVelocity(1:2) * params.dt;
    nextState(3) = currentState(3) + worldVelocity(3) * params.dt;
end
