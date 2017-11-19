function download_dependencies()
    install_dir = strsplit(mfilename('fullpath'),filesep);
    install_dir = strjoin(install_dir(1:end-1), filesep);
    cd(install_dir);

    %% download bbci toolbox and unzip    
    disp('downloading bbci-toolbox')
    bbci_url = 'https://github.com/bbci/bbci_public/archive/master.zip';
    filename = fullfile('external','bbci_tmp.zip');
    outfilename = websave(filename,bbci_url);
    unzip(outfilename,'external');
    delete(outfilename);
    
    %% download head model    
    disp('downloading head-model')
    head_url = 'https://www.parralab.org/nyhead/sa_nyhead.mat';
    filename = fullfile('data','sa_nyhead.mat');
    % outfilename = websave(filename,head_url);
    urlwrite(head_url, filename);
    
    %% configure bbci toolbox
    addpath(fullfile('external','bbci_public-master'))
    startup_bbci_toolbox()
end