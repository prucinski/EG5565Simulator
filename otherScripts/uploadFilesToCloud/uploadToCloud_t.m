function [] = uploadToCloud_t(key, secretID)
    % Extract username and password from command-line arguments
    % Display the extracted username and password
    disp(['Key ID provided: ', key]);
    disp(['secretID provided: ', secretID]);   

    arg = strcat('uploadToCloud("', key,'", "', secretID,'")');
    t = timer('ExecutionMode', 'fixedRate', 'Period', 15, 'TimerFcn', arg); %every 15 seconds
    t.TasksToExecute = 10;      %do it for no longer than 10 timies
    start(t)                    %Start the task and do it 10 times

    %can be ended prematurely with command "stop(timerfindall)"
end