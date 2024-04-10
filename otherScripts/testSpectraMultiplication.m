path = "C:\EG5565\simulatedResponses\strain_temp\strain_20240410_134921.csv";
fbgRespAsMatrix = readmatrix(path);        


x = fbgRespAsMatrix(:, 1);
y = fbgRespAsMatrix(:, 2);


laserY = fbgRespAsMatrix(:, 2);
%simulate laser response
for i = 1:numel(x)
    if i > 800 && i < 1200
        laserY(i) = 0.8;
    else
        laserY(i) = 0;
    end
end

figure;

subplot(1, 3, 1); % 1 row, 3 columns, plot 1
plot(x, y);
xlim([fbgRespAsMatrix(1,1) fbgRespAsMatrix(end,1)]);
ylabel("Relative intensity")
ylim("auto");
title('FBG spectral response');

% Plot 2
subplot(1, 3, 2); % 1 row, 3 columns, plot 2
plot(x, laserY);
title('Laser spectrum');
ylim([0 1]);

% Plot 3
subplot(1, 3, 3); % 1 row, 3 columns, plot 3
plot(x, y.*laserY);
title('Resulting spectrum');

