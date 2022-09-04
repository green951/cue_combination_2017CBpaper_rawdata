clear all;
close all;
mu1s=0:.02:2;
%mu2s=mu1s(    Shuffle(   1:1:length(mu1s)   )       );


std1=.8;
std2=.5;
%std_opt=sqrt(std1^2*std2^2/(std1^2+std2^2));
w1=std2^2/(std1^2+std2^2);
conflict=2;
Ntrial=length(mu1s);

iterations=1000;

for n=1:iterations,
    
    for i=1:Ntrial     
        cue1(i)=randn()*std1+mu1s(i);
        cue2(i)=randn()*std2+0;
        cue3(i)=cue1(i)*w1+cue2(i)*(1-w1);  % combination condition, optimal integration       
        cue4(i)=(cue1(i)+conflict)*w1+cue2(i)*(1-w1); % conflict condition, optimal integration     
    end
   
    std1_act(n)=std(cue1);
    std2_act(n)=std(cue2);
    std3_act(n)=std(cue3);
    std4_act(n)=std(cue4);
    rr_act(n)=std2_act(n)^2/(std2_act(n)^2+std1_act(n)^2);
    d1_comb(n)=abs(mean(cue1)-mean(cue3));
    d2_comb(n)=abs(mean(cue2)-mean(cue3));
    d1_conf(n)=abs(mean(cue1)-mean(cue4));
    d2_conf(n)=abs(mean(cue2)-mean(cue4));
    rp_comb(n)=d2_comb/(d1_comb+d2_comb);
    rp_conf(n)=d2_conf/(d1_conf+d2_conf);
    std_opt(n)=sqrt(std1_act(n)^2*std2_act(n)^2/(std1_act(n)^2+std2_act(n)^2));
    
end

std_means=[mean(std1_act); mean(std2_act); mean(std3_act); mean(std4_act); mean(std_opt)];
weight_means=[mean(rp_comb);mean(rp_conf);(mean(rp_comb)+mean(rp_conf))/2; mean(rr_act)];


