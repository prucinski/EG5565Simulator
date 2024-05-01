function any = setAppToStatus(status, app)
    %unsafe
    stringBuilder = "";
    if(status == 0)
        stringBuilder = 'Container is NOT ready for transportation.';
        app.Lamp.Color = 'red';  
    %safe
    elseif(status == 2)
        stringBuilder = 'Container is ready for transportation.';
        app.Lamp.Color = 'green';   
    %unknown
    elseif(status == -1)
        stringBuilder = 'Status unknown.';
        app.Lamp.Color = 'white';  
    elseif(status == 1)
        stringBuilder = 'Warning: the container parameters are not fully within safe limits.';
        app.Lamp.Color = 'yellow';  
    end
    app.StatusTextArea.Value = char(sprintf("% s \nSpectrum created on: %s", stringBuilder, app.dateToDisplayInBox));
    any = 0;
end