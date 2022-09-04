 clear all;
 close all;
%  analyzing disorientation
cd C:\Users\xiaoli\Documents\experiments\CueCombination\rawdata
  IDs=[9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
 Names={'Kathleen' 'David' 'Alyson' 'Chien' 'Jessica' 'Michelle' 'Chen22' 'Emma' 'Morgan' 'Nicole' 'Zach' 'Bailey' 'Kyle' 'Brooks' 'Harding' 'Matthew' 'Torren'};

% IDs=[ 21 22 23 24 25];
% Names={'Kyle' 'Brooks' 'Harding' 'Matthew' 'Torren'};


index1=0;  % vision
index2=0;  % motion
index3=0;  % combination
index4=0;  % conflict

index11=0; %vision, first pointing
index12=0; %vision, last pointing

for i=1:length(IDs),
    id=IDs(i);
    name=Names{i};
    filename=['VR_CueCombination_Facing' num2str(id) '_' name '.dat'];
    data=textread(filename);
   
    trialNum=data(:,2);
    condition=data(:,3);
    CorAngle=data(:,4);
    RespAngle=data(:,5);
    
    DispAngle=RespAngle-CorAngle;
    for n=1:length(DispAngle),
        if DispAngle(n)>180
            DispAngle(n)=DispAngle(n)-360;
        elseif DispAngle(n)<-180
            DispAngle(n)=DispAngle(n)+360;
        end
    end
    
    
    vision=find(condition==1);
    motion=find(condition==2);
    comb=find(condition==3);
    conf=find(condition==4);
    
    % calculate sd for each subject
    AngleVision=DispAngle(vision);
    AngleMotion=DispAngle(motion);
    AngleComb=DispAngle(comb);
    AngleConf=DispAngle(conf);
    
    
    
    % concatenating trials across subjects
    % pick up the first landmark-pointng response
    for j=1:length(vision),
        % all pointing responses
        index1=index1+1
        AngleVisionAll(index1)=DispAngle(vision(j));
        
        % only the first response in a trial
        if vision(j)==1 || (vision(j)>1 && trialNum(vision(j))~=trialNum(vision(j)-1))
            index11=index11+1;
            AngleVisionAllFirst(index11)=DispAngle(vision(j));
            vision(j)
         end
        
        %only the last response in a trial
        if vision(j)==length(trialNum) || (vision(j)<length(trialNum) && trialNum(vision(j))~=trialNum(vision(j)+1))
            index12=index12+1;
            AngleVisionAllLast(index12)=DispAngle(vision(j));
            vision(j)
        end
       
       
       
    end
    
    for j=1:length(motion),
        index2=index2+1;
        AngleMotionAll(index2)=DispAngle(motion(j));
    end
    
    for j=1:length(comb),
        index3=index3+1;
        AngleCombAll(index3)=DispAngle(comb(j));
    end
    
    for j=1:length(conf),
        index4=index4+1;
        AngleConfAll(index4)=DispAngle(conf(j));
    end
      
end
    
%%
figure(1)
subplot(2,2,1)
 [f,xi] = ksdensity(AngleVisionAll);
 plot(xi,f,'LineWidth',2);
 title('vision')
 axis([-200,200,0,3/100]);
 
 subplot(2,2,2)
 [f,xi] = ksdensity(AngleMotionAll);
 plot(xi,f,'LineWidth',2);
  title('motion')
  axis([-200,200,0,3/100]);
 
 subplot(2,2,3)
 [f,xi] = ksdensity(AngleCombAll);
 plot(xi,f,'LineWidth',2);
  title('combination')
  axis([-200,200,0,3/100]);
 
 subplot(2,2,4)
 [f,xi] = ksdensity(AngleConfAll);
 plot(xi,f,'LineWidth',2);
  title('conflict')
  axis([-200,200,0,3/100]);
  
  %%  
  figure(2)
%   subplot(1,3,1)
%   [f,xi]=ksdensity(AngleVisionAll);
%  plot(xi,f,'LineWidth',2);
%  title('vision (all)')
%  axis([-200,200,0,2/100]);
%  
%  
%  subplot(1,3,2)
%   [f,xi]=ksdensity(AngleVisionAllFirst);
%  plot(xi,f,'LineWidth',2);
%  title('vision (first pointing)')
%  axis([-200,200,0,2/100]);
 
%  subplot(1,3,3)
  [f,xi]=ksdensity(AngleVisionAllLast);
 plot(xi,f,'LineWidth',2);
 title('vision (last pointing)')
 axis([-200,200,0,2/100]);