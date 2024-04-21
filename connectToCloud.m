function [] = connectToCloud(key, secretID)


    setenv("AWS_ACCESS_KEY_ID",accessKeyID);
    setenv("AWS_SECRET_ACCESS_KEY",sAccessKey);
    setenv("AWS_DEFAULT_REGION","eu-north-1");
    ds_ref = tabularTextDatastore("s3://eg5565data/temp");
    ds_strain = tabularTextDatastore("s3://eg5565data/strain_temp");
    a = char(ds_ref.Files(end));
    newestStrain =readtable(a);
    %newestTemp = readtable(ds_strain.Files(1));
    plot(newestStrain);

end