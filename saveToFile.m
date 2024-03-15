%save to a file. Extension should be provided in fileName
function success = saveToFile(x, y, path, fileName)
    try
        matrixToSave = transpose([x; y]);
        pathAndFile = append(path,fileName);
        disp(pathAndFile);
        writematrix(matrixToSave, pathAndFile);
        success = 1;
    catch   %should probably do better error catching with this
        success = 0;
    end
end
    