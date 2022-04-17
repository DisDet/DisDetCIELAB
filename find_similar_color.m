
function [color,p]=find_similar_color(data,m)
%data 1-3:XYZ; 4-6:point normal 7-9:lab values; 10:discontinuity set ID
%m threshold
ii=1;
[n1,~]=size(data);
color=data;

point = getpointsXYZ(data,1) ;
[~,b1]=ismember(point,data(:,1:3),'rows');
C1=color(b1,7:9);
    for kk=1:n1
        C2=color(kk,7:9);
        de=deltaE2000(C1,C2);
        R=real(de);
        X=imag(de);
        if X~=0
            color(kk,11)=ii;
        elseif R<=m
           color(kk,11)=ii;
            else
             color(kk,11)=ii+1;
        end 
    end
n=size(find(color(:,11)==1),1);
p=zeros(n+1,11);
p(1,:)=color(b1,:);
p(2:end,:)=color(find(color(:,11)==1),:); %#ok<*FNDSB>
figure;
pcshow(p(:,1:3),p(:,4:6));
view(60,30)
grid on;
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
axis equal;
  

