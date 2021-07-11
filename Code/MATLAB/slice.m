function [SLICE,xaxis,yaxis]=slice(x,y,z,wdir,fname)

 mesh = readFile(strcat(wdir,'/stationCoord.txt'),0,3);
 field = readFile(strcat(wdir,'/',fname),0,3);
 field = field(:,3);

 Z1 = find( z(1) >= mesh(:,3),1,'last'); 
 Z2 = find( z(end) <= mesh(:,3),1);
 if (length(Z1) == 0)
  Z1 = 1;
 end
 if (length(Z2) == 0)
  Z2 = length(mesh(:,3));
 end 

 if ( length(z) == 1 )
   cstFound = 1;
   if ( abs(mesh(Z1,3)-z(1)) < abs(mesh(Z2,3)-z(end)) )
     Z2 = Z1;
   else
     Z1 = Z2;
   end
 else
   cstFound = 0;
 end

 SLICE = []; yaxis = []; xaxis = [];
 i = Z1;
 while i <= Z2
   I = find( mesh(i,3) == mesh(:,3) );
   subMesh = mesh(I,1:2);
   subField = field(I);
   X1 = find( x(1) >= subMesh(:,1),1,'last');
   X2 = find( x(end) <= subMesh(:,1),1);
   if (length(X1) == 0)
     X1 = 1;
   end
   if (length(X2) == 0)
     X2 = length(subMesh(:,1));
   end
   if ( length(x) == 1 && abs(subMesh(X1,1)-x(1)) < abs(subMesh(X2,1)-x(end)) )
     X2 = X1;
   elseif ( length(x) == 1 )
     X1 = X2;
   end 
   if ( ~( length(x) == 1 && ~cstFound ) )
     xaxis = [];
   end
   SLICE_z = [];
   j = X1;
   while j <= X2
    J = find( subMesh(j,1) == subMesh(:,1) );
    subSubMesh = subMesh(J,2);
    subSubField = subField(J);
    Y1 = find( y(1) >= subSubMesh,1,'last');
    Y2 = find( y(end) <= subSubMesh,1);
    if (length(Y1) == 0)
     Y1 = 1;
    end
    if (length(Y2) == 0)
     Y2 = length(subSubMesh);
    end
    if ( length(y) == 1 && abs(subSubMesh(Y1)-y(1)) < abs(subSubMesh(Y2)-y(end)) )
      Y2 = Y1;
    elseif ( length(y) == 1 )
      Y1 = Y2;
    end
    if ( length(x) == 1 && ~cstFound )
      SLICE = [SLICE,subSubField((Y1:Y2)')];
      yaxis = subSubMesh((Y1:Y2)');
    elseif ( cstFound )
      SLICE = [SLICE,subSubField((Y1:Y2)')];
      yaxis = subSubMesh((Y1:Y2)');
      xaxis = [xaxis;subMesh(j,1)];
    else
      SLICE_z = [SLICE_z,subSubField(Y1)];
      xaxis = [xaxis;subMesh(j,1)];
    end
    j = J(end)+1;
   end
   if ( length(SLICE_z) > 0 )
     SLICE = [SLICE;SLICE_z];
     yaxis = [yaxis;mesh(i,3)];
   elseif ( length(x) == 1 && ~cstFound )
     xaxis = [xaxis;mesh(i,3)];
   end
   i = I(end)+1;
 end
end
