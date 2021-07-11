
[x,y,z]=makeMesh(1.05,1,1.5); %(1.08,1.0,2.0);

Cnew = readFile('../PPG/exp24/homTurb/20s/CONCV_20s.txt',0,3);
Cnew = Cnew(:,3);
%Cnew = readFile('../PPG/exp45/newPMV3/118s/CONCV_0.002.txt',0,3);
%Cnew = Cnew(:,3);
Cold = readFile('../PPG/exp24/homTurb_slowPM/20s/CONCV_20s.txt',0,3);
Cold = Cold(:,3);
[SLICE,xaxis,yaxis]=slice(191,0,z,'../PPG/exp24/homTurb/20s','CONC_PUFF0.txt');
%[SLICE0,xaxis0,yaxis0]=slice(993,0,z,'../PPG/exp24/homTurb_slowPM/104s','CONC_PUFF0.txt');

plot(SLICE,xaxis,'-k','lineWidth',1.1); hold on;
%plot(SLICE0,xaxis0,'--k','lineWidth',1.1); hold on;
plot(Cnew,z,'-.k','lineWidth',1.1); hold on;
plot(Cold,z,'--k','lineWidth',1.1); hold off;
ylabel('height [m]','FontSize',18,'interpreter','latex');
xlabel('c [mg/m$^3$]','FontSize',18,'interpreter','latex');

title('$x = 191$ m, $y = 0$ m','FontSize',18,'interpreter','latex');
set(gcf,'renderer','painters','Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,'conc_NHT_20s_2','-dpdf');
%saveas(gcf,'N_120s.png');

