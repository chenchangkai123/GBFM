clc;
clear;

fname = "./data/sea_real/su/data_real_2000_p.su.shot1";

file = fopen(fname,"r");

data=fread(file,"float");
data_real=reshape(data,[6060,419]);.

data_real_x=data_real(61:4060,:);

figure(2)
imagesc(data_real_x);

shading flat;
xlabel('x/m','FontSize',15,'FontAngle','italic','FontWeight','bold');
ylabel('z/m','FontSize',15,'FontAngle','italic','FontWeight','bold');
c = colorbar('southoutside');

caxis([-0.0000002 0.0000002]);
c.Label.String = 'AMP';
colormap('gray');
c.FontSize=15;
c.FontAngle="italic";
c.FontWeight="bold";
set(gca,'ZDir','reverse','FontSize',10,'FontAngle','italic','FontWeight','bold');