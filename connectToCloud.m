%returns 1 if connection successful, 0 if unsuccessful and updates the
%file
%
function status = connectToCloud(key, secretID, app)

    setenv("AWS_ACCESS_KEY_ID",key);
    setenv("AWS_SECRET_ACCESS_KEY",secretID);
    setenv("AWS_DEFAULT_REGION","eu-north-1");
    try
        ds_ref = tabularTextDatastore("s3://eg5565data/temp");
        ds_strain = tabularTextDatastore("s3://eg5565data/strain_temp");
    catch
        errordlg("Error: wrong credentials provided", "Login error");
        if(~isempty(timerfindall))
            stop(timerfindall);
        end
        status = 0;
        return;
    end
    pathAndFilenameRef = char(ds_ref.Files(end));  %get the newest (last)
    pathAndFilenameStrain = char(ds_strain.Files(end));
    openFBGFile(pathAndFilenameStrain,"" , pathAndFilenameRef,"", app);
    status = 1;

end