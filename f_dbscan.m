function [ T ] = f_dbscan( A , eps, ppcluster)
% [ T, eps ] = f_dbscan( A , npb, ppcluster)
% Input:
% - A: matriz con las coordenadas de los puntos 
% - eps: radio para b¨²squeda de vecinos   
% - ppcluster: n m¨ªnimo de puntos por cl¨²ster 
% Output:
% - T: cl¨²sters asignados a cada vecino T=zeros(n,1); [n,d]=size(A); 
%   Copyright (C) {2015}  {Adri¨¢n Riquelme Guill, adririquelme@gmail.com}



[n,d]=size(A);
h=waitbar(0,['Cluster analysis in process. ',num2str(n),' points. Please wait']);

minpts=d+1; 
T=zeros(n,1);   
maxcluster=1; %
[idx, ~] = rangesearch(A,A,eps);
for i=1:n
    NeighborPts=idx{i};
    if length(NeighborPts)>=minpts 
      
        cv=T(NeighborPts); 
        mincv=min(cv); 
        mincv2=min(cv((cv>0))); 
        maxcv=max(cv);
        if maxcv==0
            caso=0; 
        else
            if maxcv==mincv2
                caso=1; 
            else
                caso=2; 
            end
        end
        switch caso
            case 0
                T(NeighborPts)=maxcluster; 
                % T(i)=maxcluster;
                maxcluster=maxcluster+1; 
            case 1
                if mincv==0
                    T(NeighborPts(cv==0))=mincv2;
                end
            case 2
                T(NeighborPts(cv==0))=mincv2;
                b=cv(cv>mincv2);
                [~,n1]=size(b);
                aux=0;
                for j=1:n1
                    if b(j)~=aux
                        T(T==b(j))=mincv2;
                        aux=b(j);
                    end
                end
        end
    else
    end
    waitbar(i/n,h);
end
%% 

if sum(T)==0 
%If the output is empty, that is, the cluster is not found, no action is performed
else
    T2=T;
    cluster=unique(T2,'sorted');
    cluster=cluster(cluster>0); % 
    [ nclusters,~]=size(cluster);
    % Count the number of points that belong to each cluster
    A=zeros(2,nclusters);
    numeroclusters=zeros(1, nclusters);
    for ii=1:nclusters
        numeroclusters(ii)=length(find(T2(:,1)==cluster(ii,1)));
    end
    A(2,:)=cluster; A(1,:)=numeroclusters;   
    [~,IX]=sort(A(1,:),'descend'); A=A(:,IX);
    % findclusters with more than n points
    n=ppcluster;
    I=find(A(1,:)>n);
    J=find(A(1,:)<=n);
    for ii=1:length(J)
        T(T2==A(2,J(ii)))=0;
    end
    for ii=1:length(I)
        T(T2==A(2,I(ii)))=ii;
    end
end
close(h);