%%
close all
clear

%%
z_bar=linspace(0,1,250); %z_bar=z/l
t_bar=linspace(0,1.5,500); %t_bar=t/tau
[Z_bar,T_bar]=meshgrid(z_bar,t_bar);
D=.05; %dimensionless dispersion number
X=(Z_bar-T_bar)./sqrt(4.*D.*T_bar);

%%
%figure
%plot3(Z_bar(:),T_bar(:),X(:),'.k');
%surf(Z_bar,T_bar,X,X);

%%
C=1/2*(1-erf(X));

%%
h=figure;
h.Position=[100 0 360 1080];
%plot3(Z_bar(:),T_bar(:),C(:),'.k');

z_disc=(0:15)/15;
[Z_disc,T_disc]=meshgrid(z_disc,t_bar);
X_disc=(Z_disc-T_disc)./sqrt(4.*D.*T_disc);
C_disc=1/2*(1-erf(X_disc));

scram=[1,2,15,16,4,3,14,13,5,8,9,12,6,7,10,11];
z_scram=z_disc(scram);
[Z_scram,T_scram]=meshgrid(z_scram,t_bar);
X_scram=(Z_scram-T_scram)./sqrt(4.*D.*T_scram);
C_scram=1/2*(1-erf(X_scram));


%subplot(3,1,1)
subplot('Position',[0.15 .726 .75 .246])
% scatter3(Z_scram(:),T_disc(:),C_disc(:),10,C_disc(:),'filled');
scatter3(Z_disc(:)*15+1,T_scram(:),C_scram(:),10,C_scram(:),'filled');
title('Apparent Response');
xlabel({'Apparent','Sensor Index'})
xlim([1,16])
xticks([1,6,11,16])
ylabel('Time')
zlabel('Tracer Concentration')

% subplot(3,1,2)% subplot(3,1,2)
subplot('Position',[0.15 .41 .75 .246])
scatter3(Z_disc(:)*16+1,T_disc(:),C_disc(:),10,C_disc(:),'filled');
title('Data-driven Rectification')
xlabel({'Revealed','Sensor Index'})
xlim([1,16])
xticks([1,6,11,16])
xticklabels({'I','VI','XI','XVI'})
ylabel('Time')
zlabel('Tracer Concentration')

% subplot(3,1,3)% subplot(3,1,2)
subplot('Position',[0.15 .01 .75 .33])
z_scale=10;
t_scale=20;
r1=t_scale:t_scale:length(t_bar);
r2=z_scale:z_scale:length(z_bar);
surf(Z_bar(r1,r2),T_bar(r1,r2),C(r1,r2),C(r1,r2));
title('True Space-Time Response')
xlabel({'Tube Length'})
ylabel('Time')
zlabel('Tracer Concentration')

c=colorbar; c.Label.String='Tracer Concentration';
c.Location='southoutside';
%% horizontal version
h=figure;
h.Position=[2000 200 1080 360];
%plot3(Z_bar(:),T_bar(:),C(:),'.k');

z_disc=(0:15)/15;
[Z_disc,T_disc]=meshgrid(z_disc,t_bar);
X_disc=(Z_disc-T_disc)./sqrt(4.*D.*T_disc);
C_disc=1/2*(1-erf(X_disc));

scram=[1,2,15,16,4,3,14,13,5,8,9,12,6,7,10,11];
z_scram=z_disc(scram);
[Z_scram,T_scram]=meshgrid(z_scram,t_bar);
X_scram=(Z_scram-T_scram)./sqrt(4.*D.*T_scram);
C_scram=1/2*(1-erf(X_scram));


% subplot(1,3,1)
subplot('Position',[0.045 0.12 .225 .8])
% scatter3(Z_scram(:),T_disc(:),C_disc(:),10,C_disc(:),'filled');
scatter3(Z_disc(:)*15+1,T_scram(:),C_scram(:),10,C_scram(:),'filled');
title('Apparent Response');
xlabel({'Apparent','Sensor Index'})
xlim([1,16])
xticks([1,6,11,16])
ylabel('Time')
zlabel('Tracer Concentration')

% subplot(1,3,2)
subplot('Position',[.36 0.12 .225 .8])
scatter3(Z_disc(:)*16+1,T_disc(:),C_disc(:),10,C_disc(:),'filled');
title('Data-driven Rectification')
xlabel({'Revealed','Sensor Index'})
xlim([1,16])
xticks([1,6,11,16])
xticklabels({'I','VI','XI','XVI'})
ylabel('Time')
zlabel('Tracer Concentration')

% subplot(1,3,3)
% subplot('Position',[0.15 .726 .75 .246]) <----column example
subplot('Position',[.675 0.12 .303 .8])
z_scale=10;
t_scale=20;
r1=t_scale:t_scale:length(t_bar);
r2=z_scale:z_scale:length(z_bar);
surf(Z_bar(r1,r2),T_bar(r1,r2),C(r1,r2),C(r1,r2));
title('True Space-Time Response')
xlabel({'Tube Length'})
ylabel('Time')
zlabel('Tracer Concentration')

c=colorbar; c.Label.String='Tracer Concentration';

%%
cmap=parula(256);
t_ind=round(length(t_bar)/2);
c_snap=C_disc(t_ind,:);
c_ind=max(1,ceil(c_snap*256));

blc=[0 0;... % location of bottom left corner
     1 0;...
     1 1;...
     0 1;...
     0 2;...
     0 3;...
     1 3;...
     1 2;...
     2 2;...
     2 3;...
     3 3;...
     3 2;...
     3 1;...
     2 1;...
     2 0;...
     3 0];

h=figure;
h.Position=[200 200 1000 300];
for i=1:16
    subplot(1,4,1)
    rectangle('Position',[blc(i,1),blc(i,2),1,1],'FaceColor',cmap(c_ind(i),:),'EdgeColor','k',...
              'LineWidth',3)
    subplot(1,4,2:4)
    rectangle('Position',[i-1,0,1,1],'FaceColor',cmap(c_ind(i),:),'EdgeColor','k',...
              'LineWidth',3)
end
subplot(1,4,1)
xlim([0 4])
ylim([0 4])
yticks([0,2,4])
axis square
subplot(1,4,2:4)
xlim([0 16])
ylim([0 1])
axis equal
ax1=gca; ax1.XAxisLocation='origin';
set(gca,'color','none');
ax1.YAxis.Visible = 'off'; % remove y-axis


%%
xyz=hilbert_curve(z_bar);
figure
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),20,z_bar,'filled')

%%
h=figure;
h.Position=([100 100 1080 540]);
for i=1:length(t_bar)
    t=t_bar(i);
    x=(z_bar-t)./sqrt(4.*D.*t);
    c=1/2*(1-erf(x));
    subplot(1,2,1)
    plot(z_bar,c)
    xlabel('z (dimensionless length along tube')
    ylabel('c (concentration of tracer')
    title(['Dispersed Plug Flow, t/\tau = ',sprintf('%0.2f',t)])
    subplot(1,2,2)
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),20,c,'filled')
    c=colorbar; c.Label.String='c (concentration of tracer';
    title(['Hilbert Curve, t/\tau = ',sprintf('%0.2f',t)])
    drawnow
    F(i)=getframe(gcf);
end

%%
writerObj=VideoWriter('rawlings_plug_flow_hilbert.avi');
writerObj.FrameRate=60;
open(writerObj);
for i=1:length(F)
    frame=F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);

