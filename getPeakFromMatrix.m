%go through a structure and find the point which has the highest arb. value
function peakWavelength = getPeakFromMatrix(matrix)
    peakWavelength = -1;
    peakValue = -1;
    [noOfRows, noOfColumns] = size(matrix);
    for i = 1:noOfRows
        %if value is higher than one previously found, update
        if matrix(i, 2) > peakValue
            peakWavelength = matrix(i,1);
            peakValue = matrix(i, 2);
        end
    end
    
    