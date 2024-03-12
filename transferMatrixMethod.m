%% Parameters
%with app
try
    core_refractive = app.RefIndexEditField.Value;     %just for reference of variables
    core_cladding = app.CladdingIndexEditField.Value;                               
    changeInRefractiveIndex = app.DRefIndexEditField.Value; 
    gratingPeriod = app.GratingPeriodEditField.Value;                 %m
    gratingLength = app.GratingLengthEditField.Value;                 %m 
    photoelasticCoefficient = app.PeCoeffEditField.Value;
    strainInSection = app.StrainEditField.Value;
    dTemp = app.DTempEditField.Value;

%for using without the app
catch   
    core_refractive = 1.47;
    core_cladding = 1.457;                             
    changeInRefractiveIndex = 1e-4;
    gratingPeriod = 5.27214e-7;                   %m
    gratingLength = 10*1e-3;                      %m 
    photoelasticCoefficient = 0.215;
    strainInSection = 0.001;
    dTemp = 0;

end
n_1 = core_refractive;        %just for reference of variables
n_2 = core_cladding;  
dNeff = changeInRefractiveIndex;
gratingP = gratingPeriod;
L= gratingLength;  
braggL = getBraggWavelength(gratingP, core_refractive);         %m  
pe = photoelasticCoefficient;

%maximum N
M = 2*n_1*gratingLength/braggL;

%paper1: https://ijcsi.org/papers/IJCSI-9-1-2-368-374.pdf
%paper2:https://www.sciencedirect.com/science/article/pii/S235271101630022X

currentLarray = linspace(1.5494*1e-6,1.5506*1e-6,200);
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
catch
plot(currentLarray, y_result);
%figure;
%plot(currentLarray, currentLarray);
end

%% Now, the strained spectrum
y_result_strained = zeros(1,200);
j = 1;
for currentL = currentLarray
    %% Change induced by strain
    %new grating period (eq 2, paper 2)
    newGratingP = gratingP*(1 + (1 - pe)*strainInSection); 


    B = 2*pi*core_refractive/currentL; %I'll come back to this later
    dB = B - pi/newGratingP;        %TODO: This is different with strained
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
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
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    y_result_strained(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;
end

%% Plot the strained spectrum
try
plot(app.FBGStrainedGraph, currentLarray, y_result_strained);
xlim(app.UIAxes, [currentLarray(1) currentLarray(end)]);
ylim(app.UIAxes, "auto");
catch
figure;
plot(currentLarray, y_result_strained);
end





