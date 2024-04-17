function [E_f, E_h, r_f, r_p, b_rp, h, a_f, a_h, G_p, G_a] = getTermsFromTable(app)
    %best I could do for error detection
    try
      E_f = str2double(app.UITable.Data{1,2});
        if(isnan(E_f)), app.UITable.Data{1,2} = {'72'}; error('E_f is not a number');   end
      E_h = str2double(app.UITable.Data{1,5});
        if(isnan(E_h)), app.UITable.Data{1,5} = {'72'}; error('E_h is not a number');   end
    
     r_f = str2double(app.UITable.Data{4,2});
        if(isnan(r_f)), app.UITable.Data{4,2} = {'62.5'}; error('r_f is not a number');   end     
     r_p = str2double(app.UITable.Data{4,3});
        if(isnan(r_p)), app.UITable.Data{4,3} = {'125'}; error('r_p is not a number');   end
     b_rp = str2double(app.UITable.Data{4,4});
        if(isnan(b_rp)), app.UITable.Data{4,4} = {'0.2'}; error('b_rp is not a number');   end
     h = str2double(app.UITable.Data{4,5});
        if(isnan(h)), app.UITable.Data{4,5} = {'0.1'}; error('h is not a number');   end

     a_f = str2double(app.UITable.Data{5,2});
        if(isnan(a_f)), app.UITable.Data{5,2} = {'0.5'}; error('a_f is not a number');   end
     a_h = str2double(app.UITable.Data{5,5});
        if(isnan(a_h)), app.UITable.Data{5,5} = {'23'}; error('a_h is not a number');   end
     
     G_p = str2double(app.UITable.Data{3,3});
        if(isnan(G_p)), app.UITable.Data{3,3} = {'2.02'}; error('G_p is not a number');   end
     G_a = str2double(app.UITable.Data{3,4});
        if(isnan(G_a)), app.UITable.Data{3,4} = {'0.714'}; error('G_a is not a number');   end
          
    catch ME 
        errordlg("Problem with a variable in the advanced settings tab:  " + ...
                    ME.message + ".Assuming default for variable ", "Data error");
        E_h = 72;           %GPa, Young's modulus of host structure
        E_f = 72;           %GPa, Young's modulus of optical fiber
        G_a = 0.714;        %MPa, shear modulus of adhesive
        G_p = 2.02;         %MPa, shear modulus of coating
        r_p = 125;          %micrometers, radius of coating
        r_f = 62.5;         %micrometers, radius of fiber
        h =   0.1;          %ometers, thickness of host material (tank by def 0.1m)
        b_rp = 0.2;         %1, gap between adhesive and fiber 
        a_h = 23;           %microe/K, coeff of thermal expansion of host
        a_f = 0.5;          %microe/K, coeff of thermal expansion of fiber
    end
end