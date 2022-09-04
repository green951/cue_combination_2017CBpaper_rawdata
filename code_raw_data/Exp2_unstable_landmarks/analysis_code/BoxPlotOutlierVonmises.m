function [output, deletenum]=BoxPlotOutlierVonmises(head_dirs)

% refer to the book 'statistical analysis of circular data' by NI Fisher
sample_size=[3 4 5 6 7 8 9 10 12 15 20];
critical_values=[.661 .73 .731 .714 .680 .648 .611 .583 .528 .464 .387];

R= sqrt( sum(cos(head_dirs))^2+ sum(sin(head_dirs))^2);
n=length(head_dirs);
index=0;
for i=1:length(head_dirs)
    data=[];
    data=head_dirs;
    data(i)=[];
    R(i)= sqrt( sum(cos(data))^2+ sum(sin(data))^2);
    M(i)=(R(i)-R+1)/(n-R);

    % find the closes sample sizes in the table
    if length(find(sample_size==n))==1
        z=find(sample_size==n);
        cutpoint(i)=critical_values(z);
    elseif length(find(sample_size==n))==0
        x=find(sample_size>n);
        y=find(sample_size<n);
        lower=y(length(y));
        upper=x(1);
        
        cutpoint(i)=critical_values(lower)+ (critical_values(upper)-critical_values(lower)) *(n-lower)/(upper-lower);
    end  
    
    if M(i)<cutpoint(i)
        index=index+1;
        output(index)=head_dirs(i);
    end
    
end
deletenum=length(head_dirs)-length(output);