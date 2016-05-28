function [DXB,DXM,DXD] =DERIV2(NT,NF,XLIQ,YVAP,X,Y,XBB,YBB,XDD,XFF,VB,BOT,R,D,F,MB,MD)
% DERIV2.M
%
%PERSAMAAN DERIVATIF HOLDUP KOMPONEN LIQUID PADA TIAP TRAY
NC=2;
for J=1:NC;
   DXB(J)=(XLIQ(1)*X(J,1)-VB*YBB(J)-BOT*XBB(J))/MB;
   N=1;
   DXM(J,1)=XLIQ(N+1)*X(J,N+1)+VB*YBB(J)-XLIQ(N)*X(J,N)-YVAP(N)*Y(J,N);
   for N=2:(NT-1);
      DXM(J,N)=XLIQ(N+1)*X(J,N+1)+YVAP(N-1)*Y(J,N-1)-XLIQ(N)*X(J,N)-YVAP(N)*Y(J,N);
      if N==NF;
         DXM(J,N)=XLIQ(N+1)*X(J,N+1)+YVAP(N-1)*Y(J,N-1)-XLIQ(N)*X(J,N)-YVAP(N)*Y(J,N)+F*XFF(J);
      end
   end
   N=NT;
   DXM(J,N)=R*XDD(J)+YVAP(N-1)*Y(J,N-1)-XLIQ(N)*X(J,N)-YVAP(N)*Y(J,N);
   DXD(J)=(YVAP(NT)*Y(J,NT)-(R+D)*XDD(J))/MD;
end
return
