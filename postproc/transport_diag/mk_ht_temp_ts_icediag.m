function temp=mk_ht_temp_ts_icediag(my_case)

ncvars =  {'temp'};

if strcmp(my_case,'qsa')
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_a';
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        temp{K} = ncread(this_file, ncvars{1});
        temp{K} = temp{K}(:,:,:);        
    end
    
elseif strcmp(my_case,'qsb') 
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_b';
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        temp{K} = ncread(this_file, ncvars{1});
        temp{K} = temp{K}(:,:,:);
    end
    
elseif strcmp(my_case,'qsc') 
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_c';
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        temp{K} = ncread(this_file, ncvars{1});
        temp{K} = temp{K}(:,:,:);
    end   
    
elseif strcmp(my_case,'qsd')
    projectdir = 'D:/Abhay/Model_Out/PF_Q_sg.exp_d';
    dinfo = dir( fullfile(projectdir, '*.nc') );
    num_files = length(dinfo);
    filenames = fullfile( projectdir, {dinfo.name} );

    for K = 1 : num_files
        this_file  = filenames{K};
        temp{K} = ncread(this_file, ncvars{1});
        temp{K} = temp{K}(:,:,:);
    end   

      
end

temp=temp{1,1};
       