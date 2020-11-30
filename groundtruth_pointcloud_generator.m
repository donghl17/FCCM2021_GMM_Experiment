data_total=[];
for k=0:1:35
    filename=['PointCloud',int2str(k),'.csv'];
    cloud = importdata(filename);
    datat=cloud.data(:,2:4);
%     S=find(datat(:,3)<2);
%     datat=datat(S,:);
    data_total=[data_total;datat];
end
load('data_total.mat');