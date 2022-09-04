function [output,deletenum,index_keep]=BoxPlotOutlierHe(he,quar)

% output  -- cleaned data
% deletenum -- how many outliers are deleted
% index_keep  -- ids of retained datapoints 



 [r c]=size(he);
% for i=1:r
%     dist(i)=sqrt(responses(i,1)^2+responses(i,2)^2);
% end

[he_sort,index]=sort(he);
quartile3=he_sort(round(length(he)*.75));
quartile1=he_sort(round(length(he)*.25));
band=(quartile3-quartile1)*quar;  % interquartile range
upper=quartile3+band;
lower=quartile1-band;

index=1;
for i=1:r
    if he(i)>lower && he(i)<upper
        output(index)=he(i);
        index_keep(index,1)=i;
        index=index+1;  
    end
end
output=output';
[r1,c1]=size(output);
deletenum=r-r1;  
end
