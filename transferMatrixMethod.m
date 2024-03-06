%% Parameters
core_refractive = 1.47; n_1 = core_refractive;        %just for reference of variables
core_cladding = 1.457; n_2 = core_cladding;                              
changeInRefractiveIndex = 1e-4; dNeff = changeInRefractiveIndex;
gratingPeriod = 5.27214e-7; gratingP = gratingPeriod;
braggL = getBraggWavelength(gratingP, core_refractive);         %m  
gratingLength = 15*1e-3; L= gratingLength;                       %m       

%maximum N
M = 2*n_1*gratingLength/braggL;

currentL = braggL;  %This will get simulated based on the input force and temperature - for now, constant

%paper: https://ijcsi.org/papers/IJCSI-9-1-2-368-374.pdf

currentLarray = linspace(1.5494*1e-6,1.5506*1e-6,200);
y_result = [];
for currentL = currentLarray
    %% Equation for fiber propagation constant B and constants k and y (eq 12)
    B = 2*pi*core_refractive/currentL;
    dB = B - pi/gratingPeriod;
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
    y_result(end+1) = real(reflectivity)^2+ imag(reflectivity)^2;

    %TODO: figure out a way to make the spectrum from the reflectivity &
    %transmissivity values for each wavelength

end

plot(currentLarray, y_result);




