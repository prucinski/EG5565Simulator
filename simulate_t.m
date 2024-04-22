%Function for runnign the simulation every 15 seconds
function [] = simulate_t(app)
    %every 15 seconds run the "simulate" function
    disp("simulate_t getting cheeky")
    t = timer('ExecutionMode', 'fixedRate', 'Period', 15, 'TimerFcn', @(~, ~)simulate(app));
    t.TasksToExecute = 50;      %do it for no longer than 50 times (obviously in a real scenario this can be set to inf)
    start(t)                    %Start the task and do it 10 times

    %can be ended prematurely with command "stop(timerfindall)"
end