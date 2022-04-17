function pos = getpointsXYZ(data,n) 
%Copyright {阿昆的科研日常 https://blog.csdn.net/qq_26447137/article/details/95502478}  
hFigure= figure;
pcshow(data(:,1:3),data(:,4:6));
view(0,90);
grid on;
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
axis equal;
datacursormode on 
dcm_obj = datacursormode(hFigure);
pos = zeros(n,3);
for i = 1:n    
    disp('Click line to display a data tip, then press Return.')  
   % Wait while the user does this.   
   pause        
   c_info = getCursorInfo(dcm_obj);  
   pos(i,:) = c_info.Position;
end
end
