%h - host structure, f - fiber, p - protective coating, a - adhesive
% I have no clue what the units for sectionL, fiberL are meant to be so I
% will assume m
function [e_t, e_m] = getThermalAndMechanical(dT, strain, fiberL, sectionL, term0, lambdaTerm, E_f, a_h, a_f) 

    L_f = 0.5*fiberL;   %ms, half the bonding length (length of bonded fibre)
    x = sectionL;   %convert to micrometers, section length
    %x = 0;

    %convert to consistent units
    E_f = E_f*1e9;  %GPa to Pa
    a_h = a_h*1e-6; %ue/K to e/K
    a_f = a_f*1e-6; %ue/K to e/K
        
    %good explaination here: https://www.mdpi.com/1424-8220/11/7/6926
    
    %term0 = (pi*r_f^2/(2*h*r_p*E_h))+1/E_f;
    %lambda Term and term0 are calculated only once for efficiency
    %% get thermal strain
    %This paper too: https://sensors.myu-group.co.jp/sm_pdf/SM1255.pdf
    term_e1 = (a_h - a_f)*dT/(E_f*term0);
    term_e2 = (1 - (cosh(lambdaTerm*x))/(cosh(lambdaTerm*L_f)));
    e_t = term_e1*term_e2 + a_f*dT;

    %% get mechanical strain
    e_m = strain/(E_f*term0)*term_e2;
    %disp("Strain transferred " + e_m);
end
