fabm=true;
year=['2013';'2014';'2015'];
month=['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'];
eom=[31 28 31 30 31 30 31 31 30 31 30 31];

for i=1
  disp(['year=' year(i,:)])
     
  for j=5 %:length(month) 
	  disp(['month=' month(j,:)])
  if j>1
    tstart=datenum(str2num(year(i,:)),str2num(month(j-1,:)),eom(j-1),1,0,0);
    tend=datenum(str2num(year(i,:)),j,eom(j),24,0,0);
  elseif j==1&i>1
    tstart=datenum(str2num(year(i-1,:)),str2num(month(12,:)),eom(12),1,0,0);
    tend=datenum(str2num(year(i,:)),j,eom(j),24,0,0);
  elseif j==1&i==1
    tstart=datenum(str2num(year(i,:)),str2num(month(j,:)),1,1,0,0);
    tend=datenum(str2num(year(i,:)),j,eom(j),24,0,0);
  end

  time=tstart:1/24:tend;
  casename=get_cstr;
 
filename=[casename '_nest_',year(i,:),'_' month(j,:)];
 
disp(['Loading ' 'Nesting/',filename,'.mat'])
struct=load(['Nesting/',filename,'.mat']);

nclm=struct.nclm;

load ngrd


end

fileprefix=['input/',filename];

disp('writing nest climatology to netcdf file')

 if fabm
   load fabm_nest
   write_FVCOM_nest(ngrd, fileprefix, nclm,time,fclm);
 else
   write_FVCOM_nest(ngrd, fileprefix, nclm,time);
 end
end
