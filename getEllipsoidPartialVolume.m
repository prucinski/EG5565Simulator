function newVolume = getEllipsoidPartialVolume(liquidHeight, height, width, lengthOfCylindrical, lengthOfElipsoidOnEnd)
    %fantastic source - https://www.had2know.org/academics/ellipse-segment-tank-volume-calculator.html
    firstTerm = acos(1-2*liquidHeight/height);
    secondTerm = (1-2*liquidHeight/height)*sqrt(4*liquidHeight/height-4*(liquidHeight^2)/(height^2));
    filledCylindricalVol = height/2*width/2*lengthOfCylindrical*(firstTerm - secondTerm);
    %ellipsoid at tend: https://math.stackexchange.com/questions/1380958/
    ellipsoidMultiplier = pi*lengthOfElipsoidOnEnd*width/2*liquidHeight^2/(3*(height/2)^2);
    filledEllipsoidVol = ellipsoidMultiplier*(3*height/2-liquidHeight); 
    newVolume = (filledCylindricalVol + filledEllipsoidVol) *1000; %times 1000 for conversion into liters

end