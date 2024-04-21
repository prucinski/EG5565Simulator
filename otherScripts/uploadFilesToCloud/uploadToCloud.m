function [] = uploadToCloud(key, secretID)
    accessKeyID = key;
    sAccessKey = secretID; %also needs access key ID (which is answer{2}, only for testing purposes
    setenv("AWS_ACCESS_KEY_ID",accessKeyID);
    setenv("AWS_SECRET_ACCESS_KEY",sAccessKey);
    setenv("AWS_DEFAULT_REGION","eu-north-1");
    ds_ref = tabularTextDatastore("s3://eg5565data/temp");
    ds_strain = tabularTextDatastore("s3://eg5565data/strain_temp");
    lengthOfDatastores = length(ds_ref.Files);
    if(lengthOfDatastores > 800)
        disp("The limit of the free tier for the PUT operation is close to be reached. Please use this more sparingly" + ...
            "for the rest of the month.");
    elseif(lengthOfDatastores > 970)
        disp("The limit of the free tier for the PUT operation has been reached. The program will stop.");
        return;
    end
    oldestFileRef = char(ds_ref.Files(lengthOfDatastores - 4));  %get a file that's oldish, cycle every 5 files
    oldestFileStrain = char(ds_strain.Files(lengthOfDatastores - 4));
    oldestFileRefTable = readmatrix(oldestFileRef);              %file retrieved
    oldestFileStrainTable = readmatrix(oldestFileStrain);              %file retrieved
    warning('off','all')
    currentDateTime = datetime('now', 'Format', 'yyyymmdd_HHMMSS');
    dateTimeAsStr = datestr(currentDateTime, 'yyyymmdd_HHMMSS');
    warning('on','all')
    newNameRef = "temp_" + dateTimeAsStr + ".csv";
    newNameStrain = "strain_" + dateTimeAsStr +".csv";
    disp("Oldest temp file renamed to " + newNameRef + " and uploaded to the cloud");
    disp("Oldest strain file renamed to " + newNameStrain + " and uploaded to the cloud")
    %writematrix(oldestFileRefTable, "s3://eg5565data/temp/" + newNameRef);
     %writematrix(oldestFileRefTable, "s3://eg5565data/strain/" + newNameStrain);
    %oldestFile_table =readtable(oldestFile);


end