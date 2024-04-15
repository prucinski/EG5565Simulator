%h - host structure, f - fiber, p - protective coating, a - adhesive
% I have no clue what the units for sectionL, fiberL are meant to be so I
% will assume m
function [e_t, e_m] = getThermalAndMechanical(dT, strain, fiberL, sectionL) 
    r_p = 125;          %micrometers, radius of coating
    r_f = 62.5;         %micrometers, radius of fiber
    h =   0.1e6;        %convert to micrometers, thickness of host material (tank by def 0.1m)
    E_h = 72;           %GPa, Young's modulus of host structure
    E_f = 72;           %GPa, Young's modulus of optical fiber
    b_rp = 0.2;         %1, gap between adhesive and fiber 
    G_a = 0.714e-3;     %MPa convert to GPa, shear modulus of adhesive
    G_p = 2.2e-3;       %MPa covert to GPa, shear modulus of coating
    L_f = 0.5e6*fiberL; %convert to micrometers, half the bonding length (length of fibre)
    S_f = sectionL*1e6;   %convert to micrometers, section length
    a_h = 23;           %microe/K, coeff of thermal expansion of host
    a_f = 0.5;          %microe/K, coeff of thermal expansion of fiber

    %% get the lambda term (paper 4, eq 4)
    term0 = (pi*r_f^2/(2*h*r_p*E_h))+1/E_f;
    term1 = 2*r_p/(pi*r_f^2)*term0;
    integralFunction = @(theta)  1./(r_p.*(1-sin(theta))./G_a+r_p./G_p.*log(r_p./r_f));
    term2 = integral(integralFunction, 0, acos(b_rp));
    lambdaTerm = sqrt(term1*term2);
    disp(lambdaTerm);
    %% get thermal strain
    %This paper too: https://sensors.myu-group.co.jp/sm_pdf/SM1255.pdf
    %God knows the units they use... Nobody cares about reproducibility
    %anymore
    term_e1 = (a_h - a_f)*dT/term0;
    term_e2 = (1 - (cosh(lambdaTerm*S_f))/(cosh(lambdaTerm*L_f)));
    e_t = term_e1*term_e2 + a_f*dT;

    %% get mechanical strain
    e_m = strain/term0*term_e2;
    disp("Strain transferred " + e_m);
end
