name='stationCoord_Table4.txt';
arcs=[50,100,200,400,800];
angle = [270,360,0,90];
spacing = [2,2,2,2,1];

coord = [];
for i = 1:length(arcs)
  phi = [angle(1):spacing(i):angle(2),spacing(i):spacing(i):angle(4)];
  coord = [coord;[arcs(i)*cos(phi'*pi/180),arcs(i)*sin(phi'*pi/180),1.5*ones(length(phi),1)]];
end

plot(coord(:,1),coord(:,2),'ok'); hold on;

towers = [325,339,353,7,21,35];
heights = [0.5,1.0,1.5,2.5,4.5,7.5,10.5,13.5,17.5];

for i = 1:length(towers)
  coord = [coord;[ones(length(heights),1)*[100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180)],heights']];
  plot(100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180),'or'); hold on;
end
hold off;

fid=fopen(name,'w');
fprintf(fid,'\t%6.2f\t%6.2f\t%6.2f\n',coord');
fclose(fid);

name='stationCoord_Fig3_exp22_exp40.txt';
arcs=[50,100,200,400,800];
angle = [330,360,0,30];;
spacing = 0.6*ones(1,length(arcs));

coord = [];
for i = 1:length(arcs)
  phi = [angle(1):spacing(i):angle(2),spacing(i):spacing(i):angle(4)];
  coord = [coord;[arcs(i)*cos(phi'*pi/180),arcs(i)*sin(phi'*pi/180),1.5*ones(length(phi),1)]];
end

plot(coord(:,1),coord(:,2),'ok'); hold on;

towers = [325,339,353,7,21,35];
heights = 0.25:0.25:17.5;

for i = 1:length(towers)
  coord = [coord;[ones(length(heights),1)*[100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180)],heights']];
  plot(100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180),'or'); hold on;
end
hold off;

fid=fopen(name,'w');
fprintf(fid,'\t%6.2f\t%6.2f\t%6.2f\n',coord');
fclose(fid);

name='stationCoord_Fig3_exp29.txt';
arcs=[50,100,200,400,800];
angle = [10,70];
spacing = 0.6*ones(1,length(arcs));

coord = [];
for i = 1:length(arcs)
  phi = angle(1):spacing(i):angle(2);
  coord = [coord;[arcs(i)*cos(phi'*pi/180),arcs(i)*sin(phi'*pi/180),1.5*ones(length(phi),1)]];
end

plot(coord(:,1),coord(:,2),'ok'); hold on;

towers = [325,339,353,7,21,35];
heights = 0.25:0.25:17.5;

for i = 1:length(towers)
  coord = [coord;[ones(length(heights),1)*[100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180)],heights']];
  plot(100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180),'or'); hold on;
end
hold off;

fid=fopen(name,'w');
fprintf(fid,'\t%6.2f\t%6.2f\t%6.2f\n',coord');
fclose(fid);

name='stationCoord_Fig4.txt';
arcs=[50,100,200,400,800];
angle = [270,360,0,90];
spacing = [2,2,2,2,1];

coord = [];
for i = 1:length(arcs)
  phi = [angle(1):spacing(i):angle(2),spacing(i):spacing(i):angle(4)];
  coord = [coord;[arcs(i)*cos(phi'*pi/180),arcs(i)*sin(phi'*pi/180),1.5*ones(length(phi),1)]];
end

plot(coord(:,1),coord(:,2),'ok'); hold on;

towers = [325,339,353,7,21,35];
heights = 0.25:0.25:17.5;

for i = 1:length(towers)
  coord = [coord;[ones(length(heights),1)*[100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180)],heights']];
  plot(100*cos(towers(i)*pi/180),100*sin(towers(i)*pi/180),'or'); hold on;
end
hold off;

fid=fopen(name,'w');
fprintf(fid,'\t%6.2f\t%6.2f\t%6.2f\n',coord');
fclose(fid);
