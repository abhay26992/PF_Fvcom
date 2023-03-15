function fout=node2cell(fin,nv)
%Usage fout=node2cell(fin,nv)
t1=fin(nv(:,1),:);
t2=fin(nv(:,2),:);
t3=fin(nv(:,3),:);
fout=(t1+t2+t3)/3;
