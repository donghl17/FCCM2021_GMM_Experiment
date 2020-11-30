%%In this roc_calculation method, all the sampled points are regarded as
%%positive in prediction (P). We compare these points to ground_truth input
%%dataset. If sampled points are close to ground_point, they will be
%%observed positive points (TP+FN). If the sampled points are too far
%%away from the ground_truth points, they will be the observed negative
%%points (FP+TN). Then, we iterate the possibility threshold to get all the result on roc curve.
%%In this way, we can draw the roc curve.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data input
SIGMA=[]; %covariance
MU=[]; %mean
CP=[]; %weights
num=10000; %number of the sampled points 
load('data_total.mat'); % we can get ground_truth points in "data_total"
%     MU=[MU;mu1];
%     SIGMA=cat(3,SIGMA,Sigma1);
%     CP=[CP,GMModel1.ComponentProportion];

%mengtecalo sampling
pointnum=round(CP.*num);
result=zeros(num,3);
occupancy=zeros(num,1);%ground truth
pro=zeros(num,1); %occupancy rpobability for all the sampled points
ptr=1;
for p=1:1:length(CP)%sampling from every component
    mu = MU(p,:);
    sigma = SIGMA(:,:,p);
    R = chol(sigma);
    temp = repmat(mu,pointnum(1,p),1) + randn(pointnum(1,p),3)*R;
    [L,useless]=size(temp);%Length is the number of the array, useless=3
    if L==0 %no point is sampled in current Gaussian component
        continue;
    else
        totalpdf=zeros(L,1); %occupancy probability for each sampled point in current component
    for pp=1:1:length(CP)%calculate occupancy probability
        Npdf=mvnpdf(temp,MU(pp,:),SIGMA(:,:,pp));
        totalpdf=totalpdf+Npdf*CP(1,pp);     
    end
    pro(ptr:ptr+L-1,1)=totalpdf;
    result(ptr:ptr+L-1,:)=temp;
    
    for pp=1:1:L
        tempdata=data_total-temp(pp,:);
        distance=sqrt(min(tempdata(:,1).^2+tempdata(:,2).^2+tempdata(:,3).^2)); %calculate the distance from sampled point to the nearest ground_truth point
        if distance<0.1
            occupancy(ptr,1)=1;
            ptr=ptr+1;
%             occupancy=[occupancy;1];
        else
            occupancy(ptr,1)=0;
            ptr=ptr+1;
%             occupancy=[occupancy;0];
        end   
    end
    end  
end
final=[result,pro,occupancy];%sampled_points+occupancy_probability+occupancy_status
final=sortrows(final,4); %probability
pcount=numel(find(occupancy==1)); % all observed positive points (TP+FN)
fcount=length(occupancy)-pcount; %all observed positive points (FP+TN)
roc=zeros(length(final),2);
px=0;%p_predict/f_real
py=0;%p_predict/p_real
for p=length(final):-1:1 %iterate possibility threshold, record roc results. https://www.jianshu.com/p/2ca96fce7e81
    if final(p,5)==1 
        py=py+1; %TP++
    else
        px=px+1; %FP++
    end
    roc(p,1)=px/fcount; 
    roc(p,2)=py/pcount; 
end
plot(roc(:,1),roc(:,2));






