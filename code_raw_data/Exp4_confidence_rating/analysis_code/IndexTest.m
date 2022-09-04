function [mean_ind,Nsub,p,t]=IndexTest(data,num)

%identify outliers

quar=3;
data_sort=sort(data);
quartile3=data_sort(round(length(data)*.75));
quartile1=data_sort(round(length(data)*.25));
band=(quartile3-quartile1)*quar;
lower=quartile1-band;
upper=quartile3+band;

index=1;
for i=1:length(data)
    if data(i)>lower  && data(i)<upper 
    
        data_cleaned(index,:)=data(i);
        index=index+1;
    end
end


[h,p,ci,stats] = ttest(data_cleaned,num);
t=stats.tstat;
mean_ind=mean(data_cleaned);
Nsub=length(data_cleaned);

end