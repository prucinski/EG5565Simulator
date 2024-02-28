function any = setAppToStatus(status, app)
    %unsafe
    if(status == 0)
        app.StatusTextArea.Value = 'Container is NOT ready for transportation.';
        app.Lamp.Color = 'red';  
    %safe
    elseif(status == 1)
        app.StatusTextArea.Value = 'Container is ready for transportation.';
        app.Lamp.Color = 'green';   
    %unknown
    elseif(status == -1)
        app.StatusTextArea.Value = 'Status unknown.';
        app.Lamp.Color = 'yellow';  
    end
    any = 0;
end