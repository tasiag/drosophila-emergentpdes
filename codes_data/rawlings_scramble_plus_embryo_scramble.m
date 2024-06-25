%% Load Embryo Data
close all
clear
load './fitted midline figures/fitted_midline.mat' data
load ./embryo_embedding.mat phi1 phi2
phi_dists=squareform(pdist([phi1,phi2]));
[max_dists,max_rows]=max(phi_dists);
[max_max_dist,max_col]=max(max_dists);
max_row=max_rows(max_col);

%% Reorganize and interpolate data
spacetime=squeeze(data(max_row,1:320,1:60))';
C=parula(64);
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(spacetime(:)),max(spacetime(:)),L),1:L,spacetime));
image = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
image=image(end:-1:1,:,:);
% imshow(image, [], 'XData', [0 1], 'YData', [0 1],'InitialMagnification','fit');

grid_line_width=0.1;

%% Plot embryo data and scrambled data
ny=30; nym1=1/ny;
ly=2; lym1=nym1/ly;
sy=1:ly:(ny*ly);
ey=ly:ly:(ny*ly);
my=(sy+ey)/2;
nx=40; nxm1=1/nx;
lx=8;  lxm1=nxm1/lx;
sx=1:lx:(nx*lx);
ex=lx:lx:(nx*lx);
mx=(sx+ex)/2;
% dx=.05;
% grid=0:dx:1-dx;
% [S,P]=meshgrid(grid,grid); %S varies with j, P varies with i

h=figure;
h.Position=1*[0 0 1600 333];
h.Color=[1 1 1];
subplot(1,4,3)
% subplot('Position',[0.5 0 .25 1])
xrange=[0 320]; yrange=[0 60];
imagesc(image)
title('(C){\it Drosophila} data','FontSize',6)
set(gca,'XTick',[.5,size(image,2)+.5],'XTickLabel',xrange)
set(gca,'YTick',[.5,size(image,1)+.5],'YTickLabel',yrange(end:-1:1))
% custom_ticks_image_linear(xrange,yrange,size(image))
custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5],'LineWidth',grid_line_width)
custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5],'LineWidth',grid_line_width)
xlabel('s')
ylabel('t')

image_observed=image;
rng(0)
keep_i=rand(ny,1)>0.25;
keep_j=rand(nx,1)>0.25;
shade=0;
shade_tile=repmat( reshape(hsv2rgb([1/12,1,1]),1,1,3) , ly , lx , 1 );
% subplot(2,2,3)
for i=1:ny
    for j=1:nx
        if keep_i(i) && keep_j(j)
        else
            image_observed(sy(i):ey(i),sx(j):ex(j),:) = ...
                shade * shade_tile + ...
                (1-shade) * image_observed(sy(i):ey(i),sx(j):ex(j),:);
        end
%         text(S(i,j)+dx/2,P(i,j)+dx/2,sprintf('%d,%d',i,j))
    end
end
% imagesc(image_observed)
% title('Partially observed data (pre-shredded)')
% set(gca,'XTick',[.5,size(image,2)+.5],'XTickLabel',xrange)
% set(gca,'YTick',[.5,size(image,1)+.5],'YTickLabel',yrange(end:-1:1))
% % custom_ticks_image_linear(xrange,yrange,size(image))
% custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5])
% custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5])
% xlabel('s')
% ylabel('t')
% % grid on
% % xticks(0:lx:4000)
% % yticks(0:ly:3000)
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% % xlabel('s')
% % ylabel('t')
% 
image_scrambled=NaN(size(image_observed));
subplot(1,4,4)
% subplot('Position',[0.75 0 .25 1])
[~,i_ind]=sort(rand(ny,1));
[~,j_ind]=sort(rand(nx,1));
for i=1:ny
    for j=1:nx
%         if keep_i(i_ind(i)) && keep_j(j_ind(j))
            image_scrambled(sy(i):ey(i),sx(j):ex(j),:)=image_observed(sy(i_ind(i)):ey(i_ind(i)),sx(j_ind(j)):ex(j_ind(j)),:);
%         else
%             plot(S(i,j)+dx/2,P(i,j)+dx/2,'.k')
%         end
%         text(S(i,j)+dx/2,P(i,j)+dx/2,sprintf('%d,%d',i_ind(i),j_ind(j)))
    end
end
imagesc(image_scrambled)
title('(D) Disorganized{\it Drosophila} data','FontSize',6)
which_tick_x=round([1,nx/5:nx/5:nx-1,nx]);
set(gca,'XTick',mx(which_tick_x),'XTickLabel',which_tick_x)
which_tick_y=round([1,ny/5:ny/5:ny-1,ny]);
set(gca,'YTick',my(which_tick_y),'YTickLabel',which_tick_y)
custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5],'LineWidth',grid_line_width)
custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5],'LineWidth',grid_line_width)
xlabel('s index')
ylabel('t index')
% grid on
% xticks(0:lx:4000)
% yticks(0:ly:3000)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% xlabel('s index')
% ylabel('t index')

% image_matrix=NaN(sum(keep_i)*ly,sum(keep_j)*lx,3);
% subplot(2,2,2)
% for i_true=1:ny
%     for j_true=1:nx
%         i=sum(keep_i(i_ind(1:i_true)));
%         j=sum(keep_j(j_ind(1:j_true)));
%         if keep_i(i_ind(i_true)) && keep_j(j_ind(j_true))
%             image_matrix(sy(i):ey(i),sx(j):ex(j),:)=image_scrambled(sy(i_true):ey(i_true),sx(j_true):ex(j_true),:);
%         end
%     end
% end
% imagesc(image_matrix)
% title('Observed data (scrambled, in matrix form)')
% which_tick_x=round([1,nx/5:nx/5:sum(keep_j)-1,sum(keep_j)]);
% set(gca,'XTick',mx(which_tick_x),'XTickLabel',which_tick_x)
% which_tick_y=round([1,ny/5:ny/5:sum(keep_i)-1,sum(keep_i)]);
% set(gca,'YTick',my(which_tick_y),'YTickLabel',which_tick_y)
% custom_grid_image(sx(2:sum(keep_j))-.5,[],'Color',1*[1.0 0.5 0.5])
% custom_grid_image([],sy(2:sum(keep_i))-.5,'Color',1*[0.5 1.0 0.5])
% xlabel('s index')
% ylabel('t index')
% % grid on
% % xticks(0:lx:(nx*sum(keep_j)))
% % yticks(0:ly:(ny*sum(keep_i)))
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% % xlabel('s index')
% % ylabel('t index')
% 

%% Generate rawlings data
z_bar=linspace(0,1,16); %z_bar=z/l
t_bar=linspace(0,1.5,25); %t_bar=t/tau
[Z_bar,T_bar]=meshgrid(z_bar,t_bar);
D=.05; %dimensionless dispersion number
X=(Z_bar-T_bar)./sqrt(4.*D.*T_bar); %placeholder
C_rawl=1/2*(1-erf(X)); %concentration
C_rawl(1,1)=0.5; % limit as t_bar goes to 0+ with z_bar=0;

%%
spacetime=C_rawl;
t=t_bar;
C=parula(64);
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(spacetime(:)),max(spacetime(:)),L),1:L,spacetime));
image = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
image = image(end:-1:1,:,:);

%%
ny=25; nym1=1/ny;
ly=1; lym1=nym1/ly;
sy=1:ly:(ny*ly);
ey=ly:ly:(ny*ly);
my=(sy+ey)/2;
nx=16; nxm1=1/nx;
lx=1;  lxm1=nxm1/lx;
sx=1:lx:(nx*lx);
ex=lx:lx:(nx*lx);
mx=(sx+ex)/2;
% dx=.05;
% grid=0:dx:1-dx;
% [S,P]=meshgrid(grid,grid); %S varies with j, P varies with i

% h=figure;
% h.Position=1*[0 0 1200 1000];
subplot(1,4,1)
xrange={0,1}; yrange=[min(t(:)) max(t(:))];
imagesc(image)
title('(A) Rawlings data','FontSize',6)
set(gca,'XTick',[.5,size(image,2)+.5],'XTickLabel',xrange)
set(gca,'YTick',[.5,size(image,1)+.5],'YTickLabel',yrange(end:-1:1))
% custom_ticks_image_linear(xrange,yrange,size(image))
custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5],'LineWidth',grid_line_width)
custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5],'LineWidth',grid_line_width)
xlabel('s')
ylabel('t')

image_observed=image;
rng(0)
keep_i=rand(ny,1)>0.25;
keep_j=rand(nx,1)>0.25;
shade=.7;
shade_tile=repmat( reshape(hsv2rgb([1/12,1,1]),1,1,3) , ly , lx , 1 );
% subplot(2,2,3)
% for i=1:ny
%     for j=1:nx
%         if keep_i(i) && keep_j(j)
%         else
%             image_observed(sy(i):ey(i),sx(j):ex(j),:) = ...
%                 shade * shade_tile + ...
%                 (1-shade) * image_observed(sy(i):ey(i),sx(j):ex(j),:);
%         end
% %         text(S(i,j)+dx/2,P(i,j)+dx/2,sprintf('%d,%d',i,j))
%     end
% end
% imagesc(image_observed)
% title('Partially observed data (pre-shredded)')
% set(gca,'XTick',[.5,size(image,2)+.5],'XTickLabel',xrange)
% set(gca,'YTick',[.5,size(image,1)+.5],'YTickLabel',yrange(end:-1:1))
% % custom_ticks_image_linear(xrange,yrange,size(image))
% custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5])
% custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5])
% xlabel('s')
% ylabel('t')
% % grid on
% % xticks(0:lx:4000)
% % yticks(0:ly:3000)
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% % xlabel('s')
% % ylabel('t')

image_scrambled=NaN(size(image_observed));
subplot(1,4,2)
[~,i_ind]=sort(rand(ny,1));
[~,j_ind]=sort(rand(nx,1));
for i=1:ny
    for j=1:nx
%         if keep_i(i_ind(i)) && keep_j(j_ind(j))
            image_scrambled(sy(i):ey(i),sx(j):ex(j),:)=image_observed(sy(i_ind(i)):ey(i_ind(i)),sx(j_ind(j)):ex(j_ind(j)),:);
%         else
%             plot(S(i,j)+dx/2,P(i,j)+dx/2,'.k')
%         end
%         text(S(i,j)+dx/2,P(i,j)+dx/2,sprintf('%d,%d',i_ind(i),j_ind(j)))
    end
end
imagesc(image_scrambled)
title('(B) Disorganized Rawlings data','FontSize',6)
which_tick_x=1:3:16;%round([1,nx/4:nx/5:nx-1,nx]);
set(gca,'XTick',mx(which_tick_x),'XTickLabel',which_tick_x)
which_tick_y=round([1,ny/5:ny/5:ny-1,ny]);
set(gca,'YTick',my(which_tick_y),'YTickLabel',which_tick_y)
custom_grid_image(sx(2:end)-.5,[],'Color',1*[1.0 0.5 0.5],'LineWidth',grid_line_width)
custom_grid_image([],sy(2:end)-.5,'Color',1*[0.5 1.0 0.5],'LineWidth',grid_line_width)
xlabel('s index')
ylabel('t index')
% grid on
% xticks(0:lx:4000)
% yticks(0:ly:3000)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% xlabel('s index')
% ylabel('t index')

% image_matrix=NaN(sum(keep_i)*ly,sum(keep_j)*lx,3);
% subplot(2,2,2)
% for i_true=1:ny
%     for j_true=1:nx
%         i=sum(keep_i(i_ind(1:i_true)));
%         j=sum(keep_j(j_ind(1:j_true)));
%         if keep_i(i_ind(i_true)) && keep_j(j_ind(j_true))
%             image_matrix(sy(i):ey(i),sx(j):ex(j),:)=image_scrambled(sy(i_true):ey(i_true),sx(j_true):ex(j_true),:);
%         end
%     end
% end
% imagesc(image_matrix)
% title('Observed data (scrambled, in matrix form)')
% which_tick_x=round([1,nx/5:nx/5:sum(keep_j)-1,sum(keep_j)]);
% set(gca,'XTick',mx(which_tick_x),'XTickLabel',which_tick_x)
% which_tick_y=round([1,ny/5:ny/5:sum(keep_i)-1,sum(keep_i)]);
% set(gca,'YTick',my(which_tick_y),'YTickLabel',which_tick_y)
% custom_grid_image(sx(2:sum(keep_j))-.5,[],'Color',1*[1.0 0.5 0.5],'LineWidth',1)
% custom_grid_image([],sy(2:sum(keep_i))-.5,'Color',1*[0.5 1.0 0.5],'LineWidth',1)
% xlabel('s index')
% ylabel('t index')
% % grid on
% % xticks(0:lx:(nx*sum(keep_j)))
% % yticks(0:ly:(ny*sum(keep_i)))
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% % xlabel('s index')
% % ylabel('t index')
% 
% 
% 
