function [] = maintainCloudConnection_t(key, secretID, app)
    % start and maintain file updates from the cloud
    %arg = strcat('connectToCloud("', key,'", "', secretID,'", "app")');
    t = timer('ExecutionMode', 'fixedRate', 'Period', 15, 'TimerFcn', {@(~, ~)connectToCloud(key, secretID, app)}); %every 15 seconds
    
    t.TasksToExecute = 50;      %do it for no longer than 50 times (obviously in a real scenario this can be set to inf)
    t.StartDelay = 15;
    start(t)                    %Start the task and do it 10 times

    %can be ended prematurely with command "stop(timerfindall)"
end

