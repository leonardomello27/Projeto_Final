appFolder = fullfile(fileparts(mfilename("fullpath")), '..');
mlappFile = fullfile(appFolder, 'ACC_SLSimApp');
run(mlappFile);
