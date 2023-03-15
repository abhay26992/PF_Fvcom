% Create monthly climatological forcing files from the 2007-2009 A4
% boundary data that were extracted and stored in separate .mat files

% Split in 3 parts. This is part (1) of (3). See further:

% --> A4_merge_clima.m (Part 2)
% --> A4_write_clima.m (Part 3)

% Part 1: A4 from 2007-2009 were not linked with dates and were therefore
% extracted clubbed together in one file with ~365 time steps. This needs
% to be split month-wise in order to prepare the mean-monthly climatology.

% Written by Abhay Prakash (abhay.prakash@natgeo.su.se)

clear all;close all;clc;

%% Splitting and saving as mean monthly climatology (for that year):

months=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];


load ./A4_2007_2017/Boundary/2007-09/pf_o_1_nest_2008__all_days.mat %File (Change as needed)
t_id=datestr(nclm.time(1:end-1)); %We exclude Jan-01-YYYY+1 for the YYYY climatology


%There are many ways to do this, I do it this way, quick and dirty but
%does the job

for ii=1:size(months,1)
    kk=0;
    for jj=1:length(t_id)
        if contains(t_id(jj,:),months(ii,:),'IgnoreCase',1)==1
            kk=kk+1;
            [temp{ii}(:,:,kk)]=nclm.temp(:,:,jj);
            [salt{ii}(:,:,kk)]=nclm.salt(:,:,jj);
            [u{ii}(:,:,kk)]=nclm.u(:,:,jj);
            [v{ii}(:,:,kk)]=nclm.v(:,:,jj);
            [ubar{ii}(:,kk)]=nclm.ubar(:,jj);
            [vbar{ii}(:,kk)]=nclm.vbar(:,jj);
            [zeta{ii}(:,kk)]=nclm.zeta(:,jj);
        end
    end
    
    m_temp(:,:,ii)=mean(temp{1,ii},3);
    m_salt(:,:,ii)=mean(salt{1,ii},3);
    m_u(:,:,ii)=mean(u{1,ii},3);
    m_v(:,:,ii)=mean(v{1,ii},3);
    m_ubar(:,ii)=mean(ubar{1,ii},2);
    m_vbar(:,ii)=mean(vbar{1,ii},2);
    m_zeta(:,ii)=mean(zeta{1,ii},2);
    
end

%% Save as struct:

clima.temp=m_temp;
clima.salt=m_salt;
clima.u=m_u;
clima.v=m_v;
clima.ubar=m_ubar;
clima.vbar=m_vbar;
clima.zeta=m_zeta;

save ./A4_2007_2017/Boundary/2007-09/clima_08 clima; 

%% Quick check:

t=mean(m_temp,1);
t_p=squeeze(t(:,1,:)); %Seasonality check : Layer 1, all time-steps

plot(t_p);grid on;         
            
    