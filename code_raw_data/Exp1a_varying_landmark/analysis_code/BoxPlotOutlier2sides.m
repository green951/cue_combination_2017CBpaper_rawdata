function [output_left,output_right,outliers_left,outliers_right]=BoxPlotOutlier2sides(responses_left, responses_right,quar)

[r1 c]=size(responses_left);
centers1=mean(responses_left);
if c==2
    for i=1:r1
        dist_left(i)=sqrt((responses_left(i,1)-centers1(1))^2+(responses_left(i,2)-centers1(2))^2);
    end
else
    for i=1:r1
        dist_left(i)=abs(    responses_left(i,1)-centers1(1)  );
    end   
end

[r2 c]=size(responses_right);
centers2=mean(responses_right);
if c==2
    for i=1:r2
        dist_right(i)=sqrt((responses_right(i,1)-centers2(1))^2+(responses_right(i,2)-centers2(2))^2);
    end
else
    for i=1:r2
        dist_right(i)=abs(    responses_right(i,1)-centers2(1)  );
    end
end

%side_index=[ones(1,r1),ones(1,r2)*2];
dist=[dist_left,dist_right];
[dist_sort,index]=sort(dist);
%side_index_sort=side_index(index);

% %%%%% old method of finding 1st and 3rd quartile
% quartile3=dist_sort(round(length(dist)*.75));
% quartile1=dist_sort(round(length(dist)*.25));

%%%%% probably more appropriate method of finding quartiles %%%%%% 
if mod(length(dist_sort),2)==0  % if the list length is even
     a=length(dist_sort)/2;
    if mod(a,2)==1
        quartile1=dist_sort((a+1)/2);
        quartile3=dist_sort((a+1)/2+a);
    elseif mod(a,2)==0
        quartile1=(dist_sort(a/2)+dist_sort((a+2)/2))/2;
        quartile3=(dist_sort(a/2+a)+dist_sort((a+2)/2)+a)/2;
    end
elseif mod(length(dist_sort),2)==1   % if the list length is odd
    n=floor(length(dist_sort)/4);  % len=4*n+remain
    remain=mod(length(dist_sort),4);
    if remain==1
        quartile1=0.25*dist_sort(n)+0.75*dist_sort(n+1);
        quartile3=0.75*dist_sort(3*n+1)+0.25*dist_sort(3*n+2);
    elseif remain==3
        quartile1=0.75*dist_sort(n+1)+0.25*dist_sort(n+2);
        quartile3=0.25*dist_sort(3*n+2)+0.75*dist_sort(3*n+3);
    end
end
        

band=(quartile3-quartile1)*quar;
upper=quartile3+band;
lower=quartile1-band;

index=1;
for i=1:r1
    if dist_left(i)<upper
        %if (dist(i)>lower || dist(i)==lower)  && (dist(i)<upper || dist(i)==upper)
        output_left(index,:)=responses_left(i,:);
        index=index+1;
    end
end

index=1;
for i=1:r2
    if  dist_right(i)<upper
        %if (dist(i)>lower || dist(i)==lower)  && (dist(i)<upper || dist(i)==upper)
        output_right(index,:)=responses_right(i,:);
        index=index+1;
    end
end

[r,c1]=size(output_left);
outliers_left=r1-r;
[r,c1]=size(output_right);
outliers_right=r2-r;
end
