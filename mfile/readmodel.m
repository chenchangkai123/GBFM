clc;
clear;

% fname = "sea_data_munk_b2000.vp";
fname = "sea_data_munk_b2000.vp";
file = fopen(fname,"r");

data=fread(file,"float");
data_real=reshape(data,[300,200,200]);

data_real_x=permute(data_real(:,:,:),[3 2 1]);

n=200;
dx=10;
nx=20;
bottom_n=100;

[X,Y,Z]=meshgrid(1:dx:dx*n,1:dx:dx*n,1:dx:dx*(n+bottom_n));
% [X,Y,Z]=meshgrid(1:dx*2:dx*n,1:dx*2:dx*n,1:dx*2:dx*(n+bottom_n));
S1=1000;
S2=1000;
S3=1000;
figure(2)
slice(X,Y,Z,data_real_x,S1,S2,S3);
shading flat;
xlabel('x/m','FontSize',15,'FontAngle','italic','FontWeight','bold');
ylabel('y/m','FontSize',15,'FontAngle','italic','FontWeight','bold');
zlabel('z/m','FontSize',15,'FontAngle','italic','FontWeight','bold');
c = colorbar('southoutside');
c.Label.String = 'AMP';
c.FontSize=15;
caxis([1500 1700]);
c.FontAngle="italic";
c.FontWeight="bold";
set(gca,'ZDir','reverse','FontSize',10,'FontAngle','italic','FontWeight','bold');