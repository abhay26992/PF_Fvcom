function kb=get_kb

%Usage: kb=get_kb
%Reads kb from the input/casename_sigma.dat file

%casestr=get_cstr;
%fid=fopen(['input/' casestr '_sigma.dat']);

fid=fopen('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/input/sigma.dat','r');

test=1;
while test
 tline=fgetl(fid);
 if length(tline)>22
  if tline(1:22)=='NUMBER OF SIGMA LEVELS', test=0; end
 end
end

fclose(fid);
i=find(tline=='=');
kb=str2num(tline(i+1:end));
