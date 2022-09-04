
 
cd C:\Users\xiaoli\Documents\experiments\CueCombination\pilot_rawdata

%% proximity-based predctions

clear all;
dimension=0;  % 0-2d space; 1- x axis; 2- z axis; 99- 24 subjects

if dimension==0
% 2d space
SdVision=[];
SdMotion=[];
SdComb=[];
SdConf=[];
visionW_pr=[];


elseif dimension==1
% x axis
SdVision=[];
SdMotion=[];
SdComb=[];
SdConf=[];
visionW_pr=[];
motionW_pr=[];

elseif dimension==2
% z axis
SdVision=[];
SdMotion=[];
SdComb=[];
SdConf=[];
visionW_pr=[];
motionW_pr=[];

    end
%% group subjects based on relative reliabilities
RR=SdMotion.^2 ./ (SdVision .^2+SdMotion.^2); % landmark relative reliability
IDs=1:length(SdMotion);
[RR_sort,idx]=sort(RR);   % sort by landmark relative reliability
SdVision_sort=SdVision(idx);
SdMotion_sort=SdMotion(idx);
SdComb_sort=SdComb(idx);
%SdConf_sort=SdConf(idx);
%visionW_sort=visionW_pr(idx);
%motionW_sort=motionW_pr(idx);

SdVision_g1=SdVision_sort(1:length(SdVision)/2);
SdVision_g2=SdVision_sort(length(SdVision)/2+1:length(SdVision));

SdMotion_g1=SdMotion_sort(1:length(SdVision)/2);
SdMotion_g2=SdMotion_sort(length(SdVision)/2+1:length(SdVision));

SdComb_g1=SdComb_sort(1:length(SdVision)/2);
SdComb_g2=SdComb_sort(length(SdVision)/2+1:length(SdVision));

%SdConf_g1=SdConf_sort(1:8);
%SdConf_g2=SdConf_sort(9:16);

SdOpt_g1=sqrt(SdVision_g1.^2 .* SdMotion_g1.^2./(SdVision_g1.^2+SdMotion_g1.^2));
SdOpt_g2=sqrt(SdVision_g2.^2 .* SdMotion_g2.^2./(SdVision_g2.^2+SdMotion_g2.^2));
RR_g1=RR_sort(1:length(SdVision)/2);
RR_g2=RR_sort(length(SdVision)/2+1:length(SdVision));

% visionW_g1=visionW_sort(1:8);
% visionW_g2=visionW_sort(9:16);
% motionW_g1=motionW_sort(1:8);
% motionW_g2=motionW_sort(9:16);
%%
% weight is calculated based on relative distance proximity
w_vision=visionW_pr;
w_motion=motionW_pr;

SdMotion=SdMotion;
SdVision=SdVision;
w_vision_var=SdMotion.^2 ./ (SdMotion.^2 + SdVision.^2);
w_motion_var=1-w_vision_var;
% experimental data
for i=1:length(SdMotion)
sd_comb_pred(i)=sqrt(w_vision(i)^2*SdVision(i)^2+w_motion(i)^2*SdMotion(i)^2);
sd_comb_opt(i)=sqrt(SdVision(i)^2*SdMotion(i)^2/(SdVision(i)^2+SdMotion(i)^2));
end

% simulated data
w_vision_sim=0:0.02:1;
w_motion_sim=1-w_vision_sim;

for j=1:length(w_vision_sim)
    w_vision_adj=w_vision_var-(mean(w_vision_var)-w_vision_sim(j));
    mean(w_vision_adj)
    w_motion_adj=1-w_vision_adj;
    for i=1:length(SdMotion)
        sd_comb_sim(i,j)=sqrt(w_vision_adj(i)^2*SdVision(i)^2+w_motion_adj(i)^2*SdMotion(i)^2);
    end
    eb(j)=std(sd_comb_sim(:,j))/sqrt(length(sd_comb_sim(:,j))); % standard error of the mean
end

figure(2)
title('higher vision variance')
%errorbar(w_vision_sim,mean(sd_comb_sim,1),eb,'b-','LineWidth',1);

x=w_vision_sim;
y=mean(sd_comb_sim,1);
err=eb;
patch([x fliplr(x)],[y+err fliplr(y-err)],[0.8 0.8 0.8]);
hold on
plot(x,y,'k-');hold on

errorbar(mean(w_vision),mean(SdComb),std(SdComb)/sqrt(length(SdComb)),'r*','LineWidth',2);
hold on
herrorbar(mean(w_vision),mean(SdComb),std(w_vision)/sqrt(length(w_vision)));
hold on

errorbar(mean(w_vision),mean(SdConf),std(SdConf)/sqrt(length(SdConf)),'b*','LineWidth',2);
hold on
herrorbar(mean(w_vision),mean(SdConf),std(w_vision)/sqrt(length(w_vision)));

hold on
y=0:0.01:1;
plot(mean(w_vision_var)*ones(length(y)),y,'b-');

% plot(mean(w_vision),mean(SdComb),'r*','LineWidth',2);
% hold on
% plot(mean(w_vision),mean(SdConf),'b*','LineWidth',2);

legend('simulated', 'combination condition','conflict condition')
xlabel('landmark weight')
ylabel('SD')
axis([0 1 0 1.1]);

%%
% weight is calculated based on relative distance proximity
w_vision=[visionW_g1;visionW_g2];
w_motion=[motionW_g1;motionW_g2];

SdMotion=[SdMotion_g1;SdMotion_g2];
SdVision=[SdVision_g1;SdVision_g2];
w_vision_var=SdMotion.^2 ./ (SdMotion.^2 + SdVision.^2);
w_motion_var=1-w_vision_var;
% experimental data
[r c]=size(SdMotion);
for j=1:r
    for i=1:c
        sd_comb_pred(j,i)=sqrt(w_vision(j,i)^2*SdVision(j,i)^2+w_motion(j,i)^2*SdMotion(j,i)^2);
        sd_comb_opt(j,i)=sqrt(SdVision(j,i)^2*SdMotion(j,i)^2/(SdVision(j,i)^2+SdMotion(j,i)^2));
    end
end

figure(3)
title('variability (predicted and actual)');
x=[1 ,2];
ymean=[mean(sd_comb_opt(1,:)),mean(sd_comb_opt(2,:))];
plot(x,ymean,'k-','LineWidth',2);
hold on;
SEs=[std(sd_comb_opt(1,:))/sqrt(length(sd_comb_opt(1,:))),std(sd_comb_opt(2,:))/sqrt(length(sd_comb_opt(2,:)))];
ylower=ymean-SEs;
yupper=ymean+SEs;
plot(x,ylower,'k--','LineWidth',1.5);
hold on;
plot(x,yupper,'k--','LineWidth',1.5);
hold on;
errorbar(x,[mean(SdComb_g1),mean(SdComb_g2)],[std(SdComb_g1)/sqrt(length(SdComb_g1)),std(SdComb_g2)/sqrt(length(SdComb_g2))],'r*','LineWidth',2);
hold on
errorbar(x,[mean(SdConf_g1),mean(SdConf_g2)],[std(SdConf_g1)/sqrt(length(SdConf_g1)),std(SdConf_g2)/sqrt(length(SdConf_g2))],'b*','LineWidth',2);

figure(4)
title('vision weight (predicted and actual)');
x=[1,2];
ymean=[mean(w_vision_var(1,:)),mean(w_vision_var(2,:))];
plot(x,ymean,'k-','LineWidth',2);
hold on;
SEs=[std(w_vision_var(1,:))/sqrt(length(w_vision_var(1,:))),std(w_vision_var(2,:))/sqrt(length(w_vision_var(2,:)))];
ylower=ymean-SEs;
yupper=ymean+SEs;
plot(x,ylower,'k--','LineWidth',1.5);
hold on;
plot(x,yupper,'k--','LineWidth',1.5);
hold on;
errorbar(x,[mean(visionW_g1),mean(visionW_g2)],[std(visionW_g1)/sqrt(length(visionW_g1)),std(visionW_g2)/sqrt(length(visionW_g2))],'r*','LineWidth',2);




%% variance based predictions
for i=1:length(SdMotion_z)
% weight is calculated based on relative distance proximity
w_vision(i)=visionW_var_z(i);
w_motion(i)=motionW_var_z(i);

sd_comb_pred(i)=sqrt(w_vision(i)^2*SdVision_z(i)^2+w_motion(i)^2*SdMotion_z(i)^2);
end

% simulated data
w_vision_sim=0:0.02:1;
w_motion_sim=1-w_vision_sim;

for j=1:length(w_vision_sim)
    for i=1:length(SdMotion_z)
    sd_comb_sim(i,j)=sqrt(w_vision_sim(j)^2*SdVision_z(i)^2+w_motion_sim(j)^2*SdMotion_z(i)^2);
    end
    eb(j)=std(sd_comb_sim(:,j))/sqrt(length(sd_comb_sim(:,j)));
end

figure(2)
x=w_vision_sim;
y=mean(sd_comb_sim,1);
err=eb;
patch([x fliplr(x)],[y+err fliplr(y-err)],[0.8 0.8 0.8]);
hold on
plot(x,y,'k-');hold on

errorbar(mean(w_vision),mean(SdComb_z),std(SdComb_z)/sqrt(length(SdComb_z)),'r*','LineWidth',2);
hold on
herrorbar(mean(w_vision),mean(SdComb_z),std(w_vision)/sqrt(length(w_vision)));
hold on
errorbar(mean(w_vision),mean(SdConf_z),std(SdConf_z)/sqrt(length(SdConf_z)),'b*','LineWidth',2);
hold on
herrorbar(mean(w_vision),mean(SdConf_z),std(w_vision)/sqrt(length(w_vision)));
legend('simulated', 'combination condition','','conflict condition','')
xlabel('landmark relative weight (based on variances)')
ylabel('SD of double-cues condition')

