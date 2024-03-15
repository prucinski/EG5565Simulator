%repeat (or start) calculation, status unknown
setAppToStatus(-1, app);

%% Tank dimensions, as per mechanical team, in meters
%default tank
if(app.NewTankCheckBox == false)  
    height = 2.4;
    width = 2.2;
    lengthOfElipsoidOnEnd = 0.4;
    lengthOfCylindrical = 5.858 - 2*lengthOfElipsoidOnEnd;
%custom tank
else
    height = app.HeightEditField.Value;
    width = app.WidthEditField.Value;
    lengthOfElipsoidOnEnd = app.EllipsoidEndRadiusEditField.Value;
    lengthOfCylindrical = app.TotallLengthEditField.Value - 2*lengthOfElipsoidOnEnd;
end
%% Sensor design dimensions, as per electrical team - TODO
spacing = 7.75 * 1e6; % mm * nm
n = 1e-4;
braggWavelength = getBraggWavelength(spacing, n);
%braggWavelength = 1550;
%add wavelength shift to the graph
wavelengthShift = braggWavelength - app.peakWavelength*1e9;         %TODO: consistent across app
%display in the graph
textPosition = [0.95, 0.9]; % Adjust text position as needed
text(app.UIAxes, textPosition(1), textPosition(2), sprintf('Wavelength shift: %.2f nm', wavelengthShift), 'Units', 'Normalized', 'HorizontalAlignment', 'right');

%% determine total volume of the tank first
%ellipse surface area: minor axis * major axis * pi
volumeOfCylindricalSection = height/2 * width/2 * pi * lengthOfCylindrical;
%On ends, there are two halfs of an elipsoid that together actually
%make a full elipsoid
volumeOfElipsoidSection = 4/3* pi *height/2 * width/2 * lengthOfElipsoidOnEnd;
totalVolume = volumeOfCylindricalSection + volumeOfElipsoidSection;
totalVolume = totalVolume * 1000; %conversion from m^3 to liters

%% Calculate the height of liquid based off of a height/volume input

%TODO: Height is meant to be acquired from spectral data

%Option 1 - liquid height
if(strcmp(app.ChooseInputDropDown.Value, app.ChooseInputDropDown.Items(1)))
    liquidHeight = app.ValueEditField.Value/1000;
    disp(liquidHeight);
    %fantastic source - https://www.had2know.org/academics/ellipse-segment-tank-volume-calculator.html
    firstTerm = acos(1-2*liquidHeight/height);
    secondTerm = (1-2*liquidHeight/height)*sqrt(4*liquidHeight/height-4*(liquidHeight^2)/(height^2));
    filledCylindricalVol = height/2*width/2*lengthOfCylindrical*(firstTerm - secondTerm);
    %ellipsoid at tend: https://math.stackexchange.com/questions/1380958/
    ellipsoidMultiplier = pi*lengthOfElipsoidOnEnd*width/2*liquidHeight^2/(3*(height/2)^2);
    filledEllipsoidVol = ellipsoidMultiplier*(3*height/2-liquidHeight); 
    newVolume = (filledCylindricalVol + filledEllipsoidVol) *1000; %times 1000 for conversion into liters
    
%Option 2 - liquid volume
elseif(strcmp(app.ChooseInputDropDown.Value, app.ChooseInputDropDown.Items(2)))
    newVolume = app.ValueEditField.Value;
end
disp(['total volume of Tank:', num2str(totalVolume), ...
    ', current volume of Tank', num2str(newVolume)]); %apparently it's around 23186.7 liters
%paste the new volume into the application
app.VolumeGauge.Value = newVolume/totalVolume * 100;
app.VolumelitersEditField.Value = newVolume;
if(app.VolumeGauge.Value > 100)
   setAppToStatus(0, app); 
   errordlg('Volume levels above 100%','Volume Error');
end

%% Calculate the mass of the liquid
density = app.SpecificGravityEditField.Value;
mass = density * newVolume;
app.TotalMasskgEditField.Value = mass;










