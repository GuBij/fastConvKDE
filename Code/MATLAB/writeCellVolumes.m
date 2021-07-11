name='cellVolumes.txt';

[x,y,z,Vx,Vy,Vz]=makeMesh(1.05,1,1.5);

fid=fopen(name,'w');
for i = 1:length(z)
 for j = 1:length(x)
   for k=1:length(y)
    fprintf(fid,'\t%10.4f\n',Vx(j)*Vy(k)*Vz(i)); %,0.0,1.5);
   end
 end
end

fclose(fid);

