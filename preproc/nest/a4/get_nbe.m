function [nbe,isonb,isbce]=get_nbe(Mobj)

if ~exist('./nbe.mat','file')
%[nbe,isonb,isbce]=get_nbe(Mobj);
   nc=Mobj.nElems;
   nn=Mobj.nVerts;
   nv=Mobj.tri;
   nbet=zeros(nc,3);
   cells=zeros(nn,50);
   cellcnt=zeros(nn,1);
   nbe=zeros(nc,3); 
   for I=1:nc
     N1 = nv(I,1) ; cellcnt(N1) = cellcnt(N1)+1;
     N2 = nv(I,2) ; cellcnt(N2) = cellcnt(N2)+1;
     N3 = nv(I,3) ; cellcnt(N3) = cellcnt(N3)+1;
     cells(nv(I,1),cellcnt(N1)) = I;
     cells(nv(I,2),cellcnt(N2)) = I;
     cells(nv(I,3),cellcnt(N3)) = I;
   end
   %if(maxval(cellcnt) > 50)write(ipt,*)'bad',maxval(cellcnt)
   for I=1:nc
     N1 = nv(I,1);
     N2 = nv(I,2);
     N3 = nv(I,3);
     for J1 = 1:cellcnt(N1) 
       for J2 = 1:cellcnt(N2) 
         if ((cells(N1,J1) == cells(N2,J2)) & cells(N1,J1) ~= I), nbe(I,3) = cells(N1,J1); end
       end
     end
     for J2 = 1:cellcnt(N2) 
       for J3 = 1:cellcnt(N3) 
	 if ((cells(N2,J2) == cells(N3,J3)) & cells(N2,J2) ~= I), nbe(I,1) = cells(N2,J2); end
       end
     end
     for J1 = 1:cellcnt(N1) 
       for J3 = 1:cellcnt(N3) 
         if ((cells(N1,J1) == cells(N3,J3)) & cells(N1,J1) ~= I), nbe(I,2) = cells(N3,J3); end
       end
     end
   end

%!
%!----IF ELEMENT ON BOUNDARY SET isbce(I)=1 AND isonb(J)=1 FOR BOUNDARY NODES J
%
    isbce=zeros(nc,1);
    isonb=zeros(nn,1);
    for I=1:nc 
     if (min([nbe(I,1),nbe(I,2),nbe(I,3)])==0)  %ELEMENT ON BOUNDARY
       isbce(I) = 1;
       if (nbe(I,1) == 0) 
         isonb(nv(I,2)) = 1; isonb(nv(I,3)) = 1;
       end
       if (nbe(I,2) ==0) 
         isonb(nv(I,1)) = 1; isonb(nv(I,3)) = 1;
       end
       if (nbe(I,3) ==0)
         isonb(nv(I,1)) = 1; isonb(nv(I,2)) = 1;
       end
     end
    end
   %%%%%%%%open boundary 
   
   for n=1:length(Mobj.read_obc_nodes)
     I_OBC_N=Mobj.read_obc_nodes{n};
     IOBCN=length(I_OBC_N);
     for I=1:IOBCN
       isonb(I_OBC_N(I))=2;
     end
   end
   
   IBCETMP=0;
   for I=1:nc
     ITMP1=isonb(nv(I,1));
     ITMP2=isonb(nv(I,2));
     ITMP3=isonb(nv(I,3));

     if (sum(isonb(nv(I,1:3))) == 4)
       isbce(I)=2;
       IBCETMP =IBCETMP+1;
     elseif (sum(isonb(nv(I,1:3))) > 4)
%#  if defined (MULTIPROCESSOR)
%       PRINT*,'SORRY, THE BOUNDARY CELL',EGID(I),'IS NOT GOOD FOR MODEL.'
%#  else
       disp(['SORRY, THE BOUNDARY CELL ' num2str(I) ' IS NOT GOOD FOR MODEL.'])
%#  endif       
       disp('IT HAS EITHER TWO SIDES OF OPEN BOUNDARY OR ONE OPEN BOUNDARY')
       disp('AND ONE SOLID BOUNDARY. PLEASE CHECK AND MODIFIED IT.')
       %disp('THIS MESSAGE IS IN SUBROUTINE TRIANGLE_GRID_EDGE (TGE.F)')
       %PRINT*,'STOP RUNNING...'
%#  if %defined (PLBC)
%       CALL PSTOP
%#  endif

     end
   end

   for I=1:nc
     if ((nbe(I,1)+nbe(I,2)+nbe(I,3) == 0)&(isbce(I) ~= 2)), isbce(I)=3; end
     if ((nbe(I,1)+nbe(I,2) == 0)&(isbce(I) ~= 2)), isbce(I)=3; end
     if ((nbe(I,2)+nbe(I,3) == 0)&(isbce(I) ~= 2)), isbce(I)=3; end
     if ((nbe(I,1)+nbe(I,3) == 0)&(isbce(I) ~= 2)), isbce(I)=3; end
   end

  save nbe nbe isonb isbce
else
    
 load nbe 
end
   
