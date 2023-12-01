clc;
clear;

grid_size = [400, 600]; % 400x400 grid
dz = 10; % 10m spacing between depths
c0 = 1500; % Reference sound speed in m/s
epsilon = 7.37e-3; % Strength of sound speed variation
z0 = 1300; % Depth of minimum sound speed in meters
B = 1300; % Scale depth in meters

% velocity_distribution = MakeSeawaterModel(grid_size, dz, c0, epsilon, z0, B);
velocity_distribution = MakeModel(400,400,1);
% for j=1:400
%     for i=401:600
%         s=sin((j)*0.1)*20;
%         if s<i-400
%             velocity_distribution(i,j)=2500;
%         end
% 
%     end
% end

imagesc(velocity_distribution);
% 
% writefile('MunkModel_sin.vp',velocity_distribution);
