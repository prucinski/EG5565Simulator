fileAsMatrix = readmatrix(app.fbgPath + "\" + app.fbgFile);         %depends on file type - right now anything in csv format works
%fileAsMatrix(1,2)
fileLen = length(fileAsMatrix);
%Extracting wavelengths
x = fileAsMatrix(:, 1);
y = fileAsMatrix(:, 2);

plot(app.UIAxes, x, y);
xlim(app.UIAxes, [fileAsMatrix(1,1) fileAsMatrix(end,1)]);
ylim(app.UIAxes, "auto");