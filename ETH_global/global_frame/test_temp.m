clear all;close all;clc;
%%calculate data transmission in map registration
sum=0;
for k=3:2:17
    if k~=9
filename=['PointCloud',int2str(k),'.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
load(['target2_frame',int2str(k),'_xyz.mat']);
choose=(data(:,1)>xt).*(data(:,1)<(xt+6)).*(data(:,2)>yt).*(data(:,2)<(yt+6)).*(data(:,3)>zt).*(data(:,3)<(zt+6));
data=data(choose==1,:);
[hang,lie]=size(data);
sum=sum+hang;
    end
end
% sum=sum/7;

%%calculating pose error
load groundtruth_pose.mat
groundtruth=rotationlist;
load total_pose.mat
total_error=mean(sum(abs(groundtruth-rotationlist)));
load semantic_pose.mat
semantic_error=mean(sum(abs(groundtruth-rotationlist)));


rotationlist=zeros(3,7);
for k=3:2:17
    if k~=9
%     filename=['regis_semantic_',int2str(k),'.mat'];
%     load(filename);
%     R=RT(1:3,1:3)';
filename='icpList.csv';
pose = importdata(filename);
R=[pose.data(k+1,1:3);pose.data(k+1,5:7);pose.data(k+1,9:11)];
T=[pose.data(k+1,4);pose.data(k+1,8);pose.data(k+1,12)];

    thetax=atan(R(3,2)/R(3,3));
    thetay=atan(-1*R(3,1)/sqrt(R(3,2).^2+R(3,3).^2));
    thetaz=atan(R(2,1)/R(1,1));
    if k<=7
    rotationlist(1,(k-1)/2)=thetax;
    rotationlist(2,(k-1)/2)=thetay;
    rotationlist(3,(k-1)/2)=thetaz;
    else
    rotationlist(1,(k-3)/2)=thetax;
    rotationlist(2,(k-3)/2)=thetay;
    rotationlist(3,(k-3)/2)=thetaz;
    end
    end
end

% load target2.mat
% ptCloud1=pointCloud([x,y,z]);

% filename='icpList.csv';
% pose = importdata(filename);
% RT=zeros(4,4);
% RT(1,:)=pose.data(1,1:4);
% RT(2,:)=pose.data(1,5:8);
% RT(3,:)=pose.data(1,9:12);
% RT(4,:)=pose.data(1,13:16);
% RT=RT';
% tform=affine3d(RT);
% ptCloud1=pctransform(ptCloud1,tform);

for k=7:2:17
load(['target2_frame',int2str(k),'_xyz.mat'])
    
    
filename=['PointCloud0.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
choose=(data(:,1)<xt+6);%.*(data(:,1)<(xt+12)).*(data(:,2)>yt-6).*(data(:,2)<(yt+12)).*(data(:,3)>zt-12).*(data(:,3)<(zt+6));
data_ori=data(choose==1,:);
ptCloud1 = pointCloud(data);

filename=['PointCloud',int2str(k),'.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
% ptCloud2 = pointCloud(data);
choose=(data(:,1)<xt+6);%.*(data(:,1)<(xt+12)).*(data(:,2)>yt-6).*(data(:,2)<(yt+12)).*(data(:,3)>zt-12).*(data(:,3)<(zt+6));
data=data(choose==1,:);
ptCloud2 = pointCloud(data);



% filename='icpList.csv';
% pose = importdata(filename);
% RT=zeros(4,4);
% RT(1,:)=pose.data(k+1,1:4);
% RT(2,:)=pose.data(k+1,5:8);
% RT(3,:)=pose.data(k+1,9:12);
% RT(4,:)=pose.data(k+1,13:16);
% RT=RT';
% tform=affine3d(RT);
% ptCloud3=pctransform(ptCloud2,tform);


%tform = pcregistericp(ptCloud1,ptCloud2 ,'Extrapolate',true);
tform = pcregisterndt(ptCloud2,ptCloud1,0.5 );
tform = invert(tform);
RT=tform.T;
save( ['regis_semantic_',int2str(k),'.mat'],'RT');
end
aaa=1;







%%
load target2_frame1_true.mat;
%mengtecalo sampling
num=10000;
pointnum=round(cp.*num/sum(cp));
result1=[];
for p=1:1:length(cp)
    mu_t = mu(p,:);
    sigma_t = sigma(:,:,p);
    R = chol(sigma_t);
    temp = repmat(mu_t,pointnum(1,p),1) + randn(pointnum(1,p),3)*R;
    result1=[result1;temp];
end
for i=3:2:7
filename=['target2_frame',int2str(i),'.mat']; 
load(filename) ;
%mengtecalo sampling
num=10000;
pointnum=round(cp.*num/sum(cp));
result2=[];
for p=1:1:length(cp)
    mu_t = mu(p,:);
    sigma_t = sigma(:,:,p);
    R = chol(sigma_t);
    temp = repmat(mu_t,pointnum(1,p),1) + randn(pointnum(1,p),3)*R;
    result2=[result2;temp];
end
ptCloud1 = pointCloud(result1);
ptCloud2 = pointCloud(result2);
%tform = pcregistericp(ptCloud1,ptCloud2 ,'Extrapolate',true);
tform = pcregisterndt(ptCloud2,ptCloud1,0.5 );
tform = invert(tform);
RT=tform.T;
save( ['regis_',int2str(i),'.mat'],'RT');
end









ptCloud = pcread('teapot.ply');
gridStep = 0.1;
ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);

figure;
pcshow(ptCloudA);
A = [cos(pi/6) sin(pi/6) 0 0; ...
    -sin(pi/6) cos(pi/6) 0 0; ...
            0         0  1 0; ...
            5         5 10 1];
tform1 = affine3d(A);
ptCloudTformed = pctransform(ptCloud,tform1);
R=A(1:3,1:3)*A(1:3,1:3);
T=A(4,1:3).*2;
AA=zeros(4,4);
AA(4,4)=1;
AA(1:3,1:3)=R;
AA(4,1:3)=T;
tform2 = affine3d(AA);
ptCloudTformed2 = pctransform(ptCloud,tform2);
% pcshow(ptCloudTformed2);
% title('Transformed Teapot');
tform = pcregistericp(ptCloudTformed2,ptCloudTformed,'Extrapolate',true);
tform2 = invert(tform);

filename=['Hokuyo_0.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
ptCloud1 = pointCloud(data);
filename=['Hokuyo_2.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
ptCloud2 = pointCloud(data);
tform = pcregisterndt(ptCloud2,ptCloud1,0.5 );
k=1;
filename='icpList.csv';
pose = importdata(filename);
R=[pose.data(k+1,1:3);pose.data(k+1,5:7);pose.data(k+1,9:11)];
T=[pose.data(k+1,4);pose.data(k+1,8);pose.data(k+1,12)];
RR=R*tform.T(1:3,1:3)';
TT=T+tform.T(4,1:3)';
k=3;
filename='icpList.csv';
pose = importdata(filename);
R_true=[pose.data(k+1,1:3);pose.data(k+1,5:7);pose.data(k+1,9:11)];
T_true=[pose.data(k+1,4);pose.data(k+1,8);pose.data(k+1,12)];
errorR=RR-R;
errorT=TT-T;


filename=['PointCloud1.csv'];
cloud = importdata(filename);
data=cloud.data(:,2:4);
choose=(data(:,1)>2).*(data(:,1)<8).*(data(:,2)>10).*(data(:,2)<16).*(data(:,3)>9).*(data(:,3)<15);
data=data(choose==1,:);
x=data(:,1);
y=data(:,2);
z=data(:,3);
save target2.mat x y z





for k=1:1:19
filename='icpList.csv';
pose = importdata(filename);
R=[pose.data(k+1,1:3);pose.data(k+1,5:7);pose.data(k+1,9:11)];
T=[pose.data(k+1,4);pose.data(k+1,8);pose.data(k+1,12)];
noise=0.1;
R1_error=unifrnd(-1,1)*noise;
R2_error=unifrnd(-1,1)*noise;
R3_error=unifrnd(-1,1)*noise;
T1_error=unifrnd(-1,1)*noise;
T2_error=unifrnd(-1,1)*noise;
T3_error=unifrnd(-1,1)*noise;
thetax=atan(R(3,2)/R(3,3))*(1+R1_error);
thetay=atan(-1*R(3,1)/sqrt(R(3,2).^2+R(3,3).^2))*(1+R2_error);
thetaz=atan(R(2,1)/R(1,1))*(1+R2_error);
RR=zeros(3,3);
RR(1,1)=cos(thetay)*cos(thetaz);
RR(1,2)=sin(thetax)*sin(thetay)*cos(thetaz)-cos(thetax)*sin(thetaz);
RR(1,3)=cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz);
RR(2,1)=cos(thetay)*sin(thetaz);
RR(2,2)=sin(thetax)*sin(thetay)*sin(thetaz)+cos(thetax)*cos(thetaz);
RR(2,3)=cos(thetax)*sin(thetay)*sin(thetaz)-sin(thetax)*cos(thetaz);
RR(3,1)=-sin(thetay);
RR(3,2)=sin(thetax)*cos(thetay);
RR(3,3)=cos(thetax)*cos(thetay);
R=RR;
T(1,1)=T(1,1)*(1+T1_error);
T(2,1)=T(2,1)*(1+T1_error);
T(3,1)=T(3,1)*(1+T1_error);
%filename=['RT_',int2str(k),'.mat'];
save (['RT_',int2str(k),'.mat'], 'R', 'T');
end

cp=CP(:,target_find);
mu=mu1(target_find,:);
sigma=Sigma1(:,:,target_find);
save target_frame2.mat cp mu sigma;
aaa=1;