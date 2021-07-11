name='stationCoord.txt';

[x,y,z]=makeMesh(1.05,1,1.5);

fid=fopen(name,'w');
for i = 1:length(z)
 for j = 1:length(x)
   for k=1:length(y)
    fprintf(fid,'\t%6.1f\t%6.1f\t%6.2f\n',x(j),y(k),z(i));
   end
 end
end

fclose(fid);

