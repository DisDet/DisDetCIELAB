% main function
%Copyright
%{Bei Cao;Yunfeng Ge; Kaili Chen; Geng Liu; Qian Chen; Weixiang Chen;
% geyunfeng@cug.edu.cn;
% China University of Geosciences (Wuhan);
% Cite as:
% Rock Joints Detection from 3D point clouds Based on CIELAB Color Space
% Keywords: Rock joints, semi-automatic identification, color space, point clouds, laser scanning
%}
clear 
close all
clc

disp('Welcome to Use ¡¶Color Space-based Intelligent Measurement Software for Geometrical Parameters of discontinuities¡·!');
disp('Please select point cloud data:');
%% Import data

[FileName,PathName,~] = uigetfile({'*.xls';'*.xlsx';'*.txt';'*.dat';'*.pcd'},'Select point cloud data');
if(isempty(FileName) || length(FileName)<=1)
    error('No data file is selected');
end
switch(FileName(end-2:end))  
    case 'dat'
        pcData=load([PathName,FileName]);
    case 'xls'
        pcData=xlsread([PathName,FileName]);
    case 'lsx'
        pcData=xlsread([PathName,FileName]);
    case 'txt'
        FileID=fopen([PathName,FileName]);
        cell_data=textscan(FileID,'%f %f %f');
        pcData=[cell_data{1},cell_data{2},cell_data{3}];
        fclose(FileID);
    case 'pcd'
        data_temp=pcread([PathName,FileName]);
        pcData=data_temp.Location;
    otherwise
        error('Only support .dat/.xls/.xlsx/.pcd£¡');
end

if(~ismatrix(pcData)) 
    error('Data dimensions are inconsistent!');
elseif(size(pcData,2)<3) 
    error('Data must be stored separately(x,y,z)£¡');
else
    pcData=pcData(:,1:3);
    disp('Data import succeeded£¡');
    disp(['Your data dimension is:',num2str(size(pcData,1)),'x',num2str(size(pcData,2))]);
end
%Data import complete
%% Point cloud feature calculation
%Point normal 
k=input('Please enter the number of neighbor points k(recommended value :40-60):'); 
if isempty(k)  %Setting defaults
    k=60;
end
ptCloud=pointCloud(pcData(:,1:3));
pitNormal=pcnormals(ptCloud,k);
pcData(:,4:6)=pitNormal;

%lab values
pcData(:,7:9)=rgb2lab(pcData(:,4:6));

%Point curvature
[s,v]=CovarianceMatrix(ptCloud,k);
pcData(:,10)=s;

%% To remove edges according to point curvature.
ii=1;
while (ii==1)
    prompt = {'Please enter the point curvature threshold r1 (points with point curvature higher than r1 will be excluded):'};
    dlgtitle = 'Input';
    dims = [1 50];
    answer = inputdlg(prompt,dlgtitle,dims);
    r1 = str2num(answer{1});  %#ok<*ST2NM>
    if isempty(r1)  
    r1=6;
    end
    ql_pcData=pcData(find(pcData(:,10)<=r1),:); %#ok<FNDSB> 
    figure;
    pcshow(ql_pcData(:,1:3),ql_pcData(:,10)) 
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',14);
    xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
    ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
    zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
    axis equal;
    view(-115,15);
    pause;
    answer = questdlg('Whether the result of Edge removing is satisfactory?');
    switch answer
    case 'Yes'   
        ii = 0;
    case 'No'
        ii = 1;  
    end
end
disp('Edge removing is complete  £¡');
%% Grouping

rest_pcData=ql_pcData;
n=input('Please enter the number of discontinuity set:'); 
J_pcData = cell(n,1);
%%
%Repeat this step
i=input('Please enter the discontinuity set that you are currently extracting:'); 
a=1;
while a==1
    m=input('Please enter the current number of times to extract this discontinuity set:');
    N=10*i+m;
    prompt = {'Please enter the color thresholdr2:'};
    dlgtitle = 'Input';
    dims = [1 35];
    answer = inputdlg(prompt,dlgtitle,dims);
    r2 = str2num(answer{1});  %#ok<*ST2NM>
    if isempty(r2)  
        r2=9;
    end
    [color,J]=find_similar_color(rest_pcData,r2); %Jis a set of data selected this time
    pause;
    answer = questdlg('Whether the results of the selected discontinuity set are satisfactory?', ...
        'Menu', ...
        'Yes','No','Need to eliminate errors','Need to eliminate errors');
    switch answer
        case 'Yes'  
            J(:,11)=ones(size(J,1),1)*i;
            [a1,b1]=ismember(J(:,1:3),rest_pcData(:,1:3),'rows');
            rest_pcData(b1,:)=[];
            Joint{m}=J;
            eval(['J',num2str(N),'=J',';']);
            eval(['rest_pcData',num2str(N),'=rest_pcData',';']);
            figure;
            pcshow(rest_pcData(:,1:3),rest_pcData(:,4:6))
            view(-115,15)
            grid on;
            set(gca,'fontname','Times New Roman','fontsize',14);
            xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis equal;
            pause
            answer = questdlg('Whether to supplement the discontinuity data of this set?');
            switch answer
                case 'Yes'
                    a=1;
                case 'No'
                    a=0;
                    for ii=1:m
                        J_pcData{i}=[J_pcData{i};Joint{ii}]  ;
                    end
            end
        case 'No' %If you are not satisfied with the result, you re-assign r2 and re-select the points to get a discontinuity set data
            a=1;
        case 'Need to eliminate errors' 
            b=1;
            c=1;
            while b==1
                prompt = {'Please enter color thresholds r3 between error data:'};
                dlgtitle = 'Input';
                dims = [1 35];
                answer = inputdlg(prompt,dlgtitle,dims);
                r3 = str2num(answer{1});  %#ok<*ST2NM>
                [color,q_J]=find_similar_color(J,r3);
                [a2,b2]=ismember(q_J,J,'rows');
                eval(['q_J',num2str(c),'=q_J',';']);
                J(b2,:)=[]; 
                figure;
                pcshow(J(:,1:3),J(:,4:6))
                view(-115,15)
                grid on;
                set(gca,'fontname','Times New Roman','fontsize',14);
                xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
                ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
                zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
                axis equal;
                pause
                answer = questdlg('Whether the error data has been eliminated?');
                switch answer
                    case 'Yes'  
                        b=0;
                        d=0;
                    case 'No'  
                        b=1;
                        c=c+1;
                        d=0;
                    case 'Cancel'
                        b=0;
                        d=1;
                end
            end
            if  d==1&&b==0
                a=1;
            else
                J(:,11)=ones(size(J,1),1)*i;
                [a1,b1]=ismember(J(:,1:3),rest_pcData(:,1:3),'rows');
                rest_pcData(b1,:)=[];
                Joint{m}=J;
                eval(['J',num2str(N),'=J',';']);
                eval(['rest_pcData',num2str(N),'=rest_pcData',';']);
                figure;
                pcshow(rest_pcData(:,1:3),rest_pcData(:,4:6))
                view(-115,15)
                grid on;
                set(gca,'fontname','Times New Roman','fontsize',14);
                xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
                ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
                zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
                axis equal;
                pause
                answer = questdlg('Whether to supplement the discontinuity data of this set?');
                switch answer
                    case 'Yes'
                        a=1;
                    case 'No'
                        a=0;
                        for ii=1:m
                            J_pcData{i}=[J_pcData{i};Joint{ii}]  ;
                        end
                end
            end
    end
end


%% Clustering

for ii=1:n
nvecinos=k+1;
[m,~]=size(J_pcData{ii});
density=[];
if nvecinos > m
    nvecinos=m;
    [~,dist]=knnsearch(J_pcData{ii}(:,1:3),J_pcData{ii}(:,1:3),'NSMethod','kdtree','distance','euclidean','k',nvecinos);
    data=dist(:,nvecinos); %tomamos todos las distancias del 4? vecino
else
    [~,dist]=knnsearch(J_pcData{ii}(:,1:3),J_pcData{ii}(:,1:3),'NSMethod','kdtree','distance','euclidean','k',nvecinos);
    if m<5
        data=dist(:,m); 
    else                                                               
        data=dist(:,5);
    end
end
data=unique(data,'sorted'); 
eps=mean(data)+2*std(data);

d=1;
while d==1
    prompt = {'Please input ppcluster:'};
    dlgtitle = 'Input';
    dims = [1 50];
    answer = inputdlg(prompt,dlgtitle,dims);
    ppcluster = str2num(answer{1});  %#ok<*ST2NM>
    if isempty(ppcluster) 
        ppcluster=20;
    end
    J_pcData{ii}(:,12)=f_dbscan( J_pcData{ii}(:,1:3) , eps, ppcluster);
    J=J_pcData{ii}(find(J_pcData{ii}(:,12)~=0),:);%#ok<FNDSB>
    figure;
    pcshow(J(:,1:3),J(:,4:6))
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',14);
    xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
    ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
    zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
    axis equal;
    pause
    answer = questdlg('Whether the clustering results are satisfactory?');
    switch answer
        case 'Yes'  
            d=0;
            eval(['J',num2str(ii),'=J',';']);
        case 'No'  
            d=1;
    end
end
end
disp('Clustering compelted£¡');
%% Orientation

for ii=1:n
    m=max(J_pcData{ii}(:,12));
    jointData=cell(m,1);
    orientation_J=zeros(m,2);
    for jj=1:m
        jointData{jj}=J_pcData{ii}(find(J_pcData{ii}(:,12)==jj),:); %#ok<FNDSB>
        pc=pca(jointData{jj}(:,1:3));
        normal=pc(:,3)';
        [dd,dip] = OrientationM( normal(:,1),normal(:,2),normal(:,3));
        orientation_J(jj,1)=dd;
        orientation_J(jj,2)=dip;
    end 
    eval(['jointData',num2str(ii),'=jointData',';']);
    eval(['orientation_J',num2str(ii),'=orientation_J',';']);
end
%% Plot
%Draw each discontinuity set
for ii=1:n
    J=J_pcData{ii}(find(J_pcData{ii}(:,12)~=0),:); %#ok<FNDSB>
    figure;
    pcshow(J(:,1:3),J(:,4:6))
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',14);
    xlabel(gca,'X (mm)','fontname','Times New Roman','fontsize',16 );
    ylabel(gca,'Y (mm)','fontname','Times New Roman','fontsize',16 );
    zlabel(gca,'Z (mm)','fontname','Times New Roman','fontsize',16 );
    view(-115,15);
    axis equal;
end
% One color per discontinuity set
figure;
for ii=1:n
    pcshow(J_pcData{ii}(:,1:3),J_pcData{ii}(:,11))
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',14);
    xlabel(gca,'X (mm)','fontname','Times New Roman','fontsize',16 );
    ylabel(gca,'Y (mm)','fontname','Times New Roman','fontsize',16 );
    zlabel(gca,'Z (mm)','fontname','Times New Roman','fontsize',16 );
    hold on;
    axis equal;  
end
view(-115,15);

% One color per discontinuity
numJoint=0;   
for ii=1:n
    m=max(J_pcData{ii}(:,12));
    numJoint=numJoint+m;
end
cx=rand(numJoint,1);
cy=rand(numJoint,1);
cz=rand(numJoint,1);
a=1;
figure;
for ii=1:1:n
    m=max(J_pcData{ii}(:,12));
    for jj=1:m
    j=J_pcData{ii}(find(J_pcData{ii}(:,12)==jj),:); %#ok<FNDSB>    
    pcshow(j(:,1:3),[cx(a,:),cy(a,:),cz(a,:)])
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',14);
    xlabel(gca,'X (mm)','fontname','Times New Roman','fontsize',16 );
    ylabel(gca,'Y (mm)','fontname','Times New Roman','fontsize',16 );
    zlabel(gca,'Z (mm)','fontname','Times New Roman','fontsize',16 );
    view(-115,15);
    axis equal;
    hold on;
    a=a+1;
    end
end

   
   

   
   
   
   
   
   
   
   
   
   
   
   
