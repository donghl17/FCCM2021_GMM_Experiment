clear all; clc;
filename=['PointCloud6.csv'];
cloud = importdata(filename);
result=cloud.data(:,2:4);
cmap = colormap;
d = length(cmap);
zLow = min(result(:,3));
zDelta = (max(result(:,3))-min(result(:,3)))/d;
figure(4);
hold on
for ppp = 1:length(cmap)
    filter = (result(:,3) > zLow & result(:,3) <= zLow+zDelta);
    plot3(result(filter,1), result(filter,2), result(filter,3), '.', 'Color', cmap(ppp,:))
    zLow = zLow+zDelta;
end
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal
grid on
hold off
