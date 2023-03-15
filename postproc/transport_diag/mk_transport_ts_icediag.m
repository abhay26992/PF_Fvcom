function [u,v]=mk_transport_ts_icediag(ind,my_case)

% Load and crop the 2016 model outputs

%qsa and qsb = spatial distribution experiments
%qsc and qsd = temporal distribution experiments

ncvars =  {'u','v'};

if strcmp(my_case,'std')
    
    projectdir = 'D:/Abhay/Model_Out/PF_std_yr_2016'; 
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        u{K} = ncread(this_file, ncvars{1});
        u{K} = u{K}(ind,:,:);
        v{K} = ncread(this_file, ncvars{2});
        v{K} = v{K}(ind,:,:);
    end

elseif strcmp(my_case,'qsa')
    
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_a'; 
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        u{K} = ncread(this_file, ncvars{1});
        u{K} = u{K}(ind,:,:);
        v{K} = ncread(this_file, ncvars{2});
        v{K} = v{K}(ind,:,:);
    end
    
elseif strcmp(my_case,'qsb')
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_b'; 
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        u{K} = ncread(this_file, ncvars{1});
        u{K} = u{K}(ind,:,:);
        v{K} = ncread(this_file, ncvars{2});
        v{K} = v{K}(ind,:,:);
    end
    
elseif strcmp(my_case,'qsc')
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_c';  
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        u{K} = ncread(this_file, ncvars{1});
        u{K} = u{K}(ind,:,:);
        v{K} = ncread(this_file, ncvars{2});
        v{K} = v{K}(ind,:,:);
    end    
    
elseif strcmp(my_case,'qsd')
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_d';  
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );
    
    for K = 1 : num_files
        this_file  = filenames{K};
        u{K} = ncread(this_file, ncvars{1});
        u{K} = u{K}(ind,:,:);
        v{K} = ncread(this_file, ncvars{2});
        v{K} = v{K}(ind,:,:);
    end    
      
end

u=u{1,1};v=v{1,1};
       