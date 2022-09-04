function output=LandmarkRelativeProximity(ResponseCenter,TargetLandmark)

output=zeros(1,3);
%RichConfDist1(i)=sqrt((RichConfCenter(i,1)-0)^2+(RichConfCenter(i,2)-0)^2); % response center to self-motion defined target location
%RichConfDist2(i)=sqrt((RichConfCenter(i,1)-targx_conf(i))^2+(RichConfCenter(i,2)-targz_conf(i))^2);% response center to landmark defined target location
%prox_lm(i)=RichConfDist1(i)/(RichConfDist1(i)+RichConfDist2(i)); %landmark relative proximity

dist2d1=sqrt( (ResponseCenter(1)-0)^2 +(ResponseCenter(2)-0)^2);
dist2d2=sqrt( (ResponseCenter(1)-TargetLandmark(1))^2 +(ResponseCenter(2)-TargetLandmark(1))^2);
output(1)=dist2d1/(dist2d1+dist2d2);

distx1=abs(ResponseCenter(1)-0);
distx2=abs(ResponseCenter(1)-TargetLandmark(1));
output(2)=distx1/(distx1+distx2);

distz1=abs(ResponseCenter(2)-0);
distz2=abs(ResponseCenter(2)-TargetLandmark(2));
output(3)=distz1/(distz1+distz2);
end