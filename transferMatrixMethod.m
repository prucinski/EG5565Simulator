%% Parameters
%with app
try
    cla(app.FBGStrainedGraph);                         %clear the figure
    core_refractive = app.RefIndexEditField.Value;     %just for reference of variables
    core_cladding = app.CladdingIndexEditField.Value;                               
    changeInRefractiveIndex = app.DRefIndexEditField.Value; 
    gratingPeriod = app.GratingPeriodEditField.Value;                 %m
    gratingLength = app.GratingLengthEditField.Value;                 %m 
    photoelasticCoefficient = app.PeCoeffEditField.Value;
    thermalExpansionCoefficient = app.ThermalCoeffEditField.Value;    %1/K
    strainInSection = app.StrainEditField.Value;
    dTemp = app.DTempEditField.Value;

%for using without the app
catch   
    core_refractive = 1.47;
    core_cladding = 1.457;                             
    changeInRefractiveIndex = 1e-4;
    gratingPeriod = 5.27212e-7;                   %m
    gratingLength = 10*1e-3;                      %m 
    photoelasticCoefficient = 0.215;
    %the photoealastic coefficient can also be calculated from Eq 3 from
    %here: https://link.springer.com/article/10.1007/BF02323100
    thermalExpansionCoefficient = 2.6e-6;       %source - google made it up TODO
    strainInSection = 0.0005;
    dTemp = -5;

end
n_1 = core_refractive;        %just for reference of variables
n_2 = core_cladding;  
dNeff = changeInRefractiveIndex;
gratingP = gratingPeriod;
L= gratingLength;  
braggL = getBraggWavelength(gratingP, core_refractive);         %m  
pe = photoelasticCoefficient;
alpha = thermalExpansionCoefficient;
dndT = 19*alpha*n_1; %in in silica, the change of refractive index because of temperature
% accounts for 95% of the thermal response... Hence this goofy workaround
% (paper 3). "trust me, Im an engineer"

%maximum N
M = 2*n_1*gratingLength/braggL;

%paper1: https://ijcsi.org/papers/IJCSI-9-1-2-368-374.pdf
%paper2:https://www.sciencedirect.com/science/article/pii/S235271101630022X
%paper3: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=618377

currentLarray = linspace(1.549*1e-6,1.551*1e-6,1000);
y_result = zeros(1,200);
j = 1;
%% reference spectrum
for currentL = currentLarray
    %% Equation for fiber propagation constant B and constants k and y (eq 12, paper1)
    B = 2*pi*core_refractive/currentL;
    dB = B - pi/gratingP;
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
    %% Prepare for the transfer matrix method
    N = 10;
    sectionL = L/N;
    
    T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
    T_12 = -1i*k/y*sinh(y*sectionL);
    T_21 = 1i*k/y*sinh(y*sectionL);
    T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
    
    T_i = [T_11, T_12;T_21, T_22];
    
    T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
    for i = 2:N
        T= T*T_i;
    end
    %% extract reflected & transmitted amplitude from boundary conditions
    % 4.8.3 & 4.8.4 from Fiber Bragg Gratings
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    shouldBeOne = real(transmissivity)^2+ imag(transmissivity)^2 + real(reflectivity)^2+ imag(reflectivity)^2;  %all OK, distance from origin is 1
    R = k^2*sinh(y*L)^2/(dB^2*sinh(y*L)^2 + k^2*cosh(y*L)^2);
    %y_result(end+1) = R;    this has weird arifacts for some reason
    y_result(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;


end
%% Plot the reference spectrum
try
    plot(app.FBGRefGraph, currentLarray, y_result);
    xlim(app.UIAxes, [currentLarray(1) currentLarray(end)]);
    ylim(app.UIAxes, "auto");
    textPosition = [0.95, 0.9]; % Adjust text position as needed
    text(app.FBGRefGraph, textPosition(1), textPosition(2), sprintf('Bragg wavelength: %.2f nm', braggL*1e9), 'Units', 'Normalized', 'HorizontalAlignment', 'right');
catch
    plot(currentLarray, y_result);
end

%% Now, the strained spectrum 
y_result_strained = zeros(1,200);
j = 1;
for currentL = currentLarray 
    %% - change induced by strain & temperature
    %new grating period (eq 2, paper 2 + eq3, paper3)
    newGratingP = gratingP*(1 + (1 - pe)*strainInSection + (alpha +dndT/n_1)*dTemp); 
    newN = n_1; %note - there should be a change of the effective refractive index
    newBraggL = getBraggWavelength(newGratingP, newN);%note - this is only useful for uniform stress
    B = 2*pi*core_refractive/currentL; 
    dB = B - pi/newGratingP;        
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
    N = 10;
    sectionL = (2/3)*L/N;
    
    T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
    T_12 = -1i*k/y*sinh(y*sectionL);
    T_21 = 1i*k/y*sinh(y*sectionL);
    T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
    
    T_i = [T_11, T_12;T_21, T_22];
    
    T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
    for i = 2:N
        T= T*T_i;
    end
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    y_result_strained(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;
end
newBraggL_strain = newBraggL;   %for plotting later
y_result_temp = zeros(1,200);
j = 1;
for currentL = currentLarray 
    %% Change induced by only the temperature on the same fibre (no mechanical strain)
    newGratingP = gratingP*( 1 + (alpha +dndT/n_1)*dTemp); 
    newN = n_1; %note - there should be a change of the effective refractive index
    newBraggL = getBraggWavelength(newGratingP, newN);%note - this is only useful for uniform stress
    B = 2*pi*core_refractive/currentL; 
    dB = B - pi/newGratingP;        
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
    N = 10;
    sectionL = (1/3)*L/N;         %to distinguish between spectra, temperature spectrum will have half the length
    T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
    T_12 = -1i*k/y*sinh(y*sectionL);
    T_21 = 1i*k/y*sinh(y*sectionL);
    T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
    
    T_i = [T_11, T_12;T_21, T_22];
    T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
    for i = 1:N
        T= T*T_i;
    end
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    y_result_temp(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;
end
%hello git
%% Plot the strained spectrum
try
    hold(app.FBGStrainedGraph, 'on');
    plot(app.FBGStrainedGraph, currentLarray, y_result_strained);
    plot(app.FBGStrainedGraph, currentLarray, y_result_temp);
    textPosition = [0.95, 0.9]; % Adjust text position as needed
    text(app.FBGStrainedGraph, textPosition(1), textPosition(2), sprintf('Bragg wavelength (temp/temp+strain): %.2f/%.2f nm', newBraggL*1e9, newBraggL_strain*1e9), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'right');
    xlim(app.UIAxes, [currentLarray(1) currentLarray(end)]);
    ylim(app.UIAxes, "auto");
    hold(app.FBGStrainedGraph, 'off');
    app.simulatedXResp = currentLarray;
    app.simulatedYResp_strain = y_result_strained;
    app.simulatedYResp_temp = y_result_temp;
catch
    figure;
    hold on;
    plot(currentLarray, y_result_strained);
    plot(currentLarray, y_result_temp);
end





