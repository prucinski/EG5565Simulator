function [] = openFBGFile(strainFile, strainPath, refFile, refPath, app)
    if(refPath ~= "")       %for the cloud solution, we're providing all as one path
        strainFull = strainPath + "\" + strainFile;
    else
        strainFull = strainFile;
    end
    if(strainPath ~= "")
        refFull = refPath + "\" + refFile;
    else
        refFull = refFile;
    end
    fbgRespAsMatrix = readmatrix(strainFull);         %depends on file type - right now anything in csv format works
    refFbgRespAsMatrix = readmatrix(refFull);
    app.dateToDisplayInBox = extractDate(strainFile); %keep the date from file for later
    %fileAsMatrix(1,2)
    fileLen = length(fbgRespAsMatrix);
    %Extracting wavelengths
    x = fbgRespAsMatrix(:, 1);
    y = fbgRespAsMatrix(:, 2);
    refX = refFbgRespAsMatrix(:, 1);
    refY = refFbgRespAsMatrix(:, 2);
    %extract peak wavelength
    app.peakWavelength = getPeakFromMatrix(fbgRespAsMatrix);
    app.peakRefWavelength = getPeakFromMatrix(refFbgRespAsMatrix);
    % TODO: I was pondering what would happen if the x axis was different, and
    % I found that the app would break (although currently manufacturing files
    % like this is impossible in the simulator)
    cla(app.UIAxes);
    hold(app.UIAxes, 'on');
    plot(app.UIAxes, x, y);
    plot(app.UIAxes, x, refY);
    xlim(app.UIAxes, [fbgRespAsMatrix(1,1) fbgRespAsMatrix(end,1)]);
    ylim(app.UIAxes, "auto");
    legend(app.UIAxes, {'str + temp FBG', 'ref. FBG'},'Location','southwest')
end
