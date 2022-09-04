function [outputx outputz]=CoordinateTransform(respx,respz,targx,targz,originx,originz)


% [respx,respz]=CoordinateTransform_new(respx_org,respz_org,originx,originz,targx+jitterx,targz+jitterz);
    
 
% respx and respz are in a spatial coordinate with (0,0) as the origin
% this function transforms respx and respz into a spatial coordinate with
% (targx, targz) as the origin

[r c]=size(respx);
outputx=zeros(r,c);
outputz=zeros(r,c);

% for i=1:r  % cannot figure out how this works, but it works
%     angle_targ=radial2degree(atan2(targz(i),targx(i)));
%     angle_resp=atan2(respz(i),respx(i))*180/pi;
%     dist=sqrt(respx(i)^2+respz(i)^2);
%     outputx(i)=dist*cos(degree2radial((90+angle_resp-angle_targ)));
%     outputz(i)=dist*sin(degree2radial((90+angle_resp-angle_targ)))-sqrt(targx(i)^2+targz(i)^2);
% end

for i=1:r  % this works and I can tell how it works, first do the transaltion and then do the rotation
    respx_prime=respx(i)-targx(i);
    respz_prime=respz(i)-targz(i);
    
    angle1=radial2degree(atan2(respz_prime-originz,respx_prime-originx));
    
    ang_targ=radial2degree(atan2(targz(i)-originz,targx(i)-originx));
    
    angle=degree2radial(angle1+90-ang_targ);
    dist=sqrt(respx_prime^2+respz_prime^2);
    outputx(i)=dist*cos(angle);
    outputz(i)=dist*sin(angle);
end

end