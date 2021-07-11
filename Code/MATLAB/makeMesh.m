function [x,y,z,Vx,Vy,Vz]=makeMesh(r,l0,z0)

maxValAt = [997,0,8];
domRange = [[778,-104,0];[1380,104,73]];
L = [[maxValAt-domRange(1,:)];[domRange(2,:)-maxValAt]];
n = ceil(log((L/l0+1/(r-1)+0.5)/(1/(r-1)+0.5))/log(r))-1;

n(1,3) = min(floor(log(z0/l0)/log(r))-1,n(1,3));

x = [fliplr(maxValAt(1) - l0*(r.^(1:(n(1,1)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)),maxValAt(1),maxValAt(1) + l0*(r.^(1:(n(2,1)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)];
y = [fliplr(maxValAt(2) - l0*(r.^(1:(n(1,2)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)),maxValAt(2),maxValAt(2) + l0*(r.^(1:(n(2,2)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)];
if (n(1,3) > 0 )
 z = fliplr(maxValAt(3) - l0*(r.^(1:(n(1,3)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5));
 z(z < 0 ) = [];
 zz = z(1);
 z = [(0.5*z0):z0:zz,z,maxValAt(3),maxValAt(3) + l0*(r.^(1:(n(2,3)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)];
 Vz = [z0*ones(1,ceil((zz-0.5*z0)/z0)),fliplr(l0*r.^(1:(n(1,3)+1))),l0,l0*r.^(1:(n(2,3)+1))];
else
 z = [maxValAt(3),maxValAt(3) + l0*(r.^(1:(n(2,3)+1))*(1/(r-1)+0.5)-1/(r-1)-0.5)];
 Vz = [l0,l0*r.^(1:(n(2,3)+1))];
end
Vx = [fliplr(l0*r.^(1:(n(1,1)+1))),l0,l0*r.^(1:(n(2,1)+1))];
Vy = [fliplr(l0*r.^(1:(n(1,2)+1))),l0,l0*r.^(1:(n(2,2)+1))];
end
