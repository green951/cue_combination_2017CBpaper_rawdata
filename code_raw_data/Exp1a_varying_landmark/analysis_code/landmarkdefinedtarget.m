function output=landmarkdefinedtarget(xz,rotation,dimen)

[r c]=size(xz);
for i=1:r
    
    dist1(i)=sqrt(xz(i,1)^2+xz(i,2)^2);
    dist2(i)=2*dist1(i)*sin(degree2radial(rotation/2));
    if dimen==1
        output(i,1)=dist2(i)*cos(degree2radial(-rotation/2));
        output(i,2)=dist2(i)*sin(degree2radial(-rotation/2));
    elseif dimen==2
        output(i,1)=dist2(i)*cos(degree2radial(-rotation/2));
    elseif dimen==3
        output(i,1)=dist2(i)*sin(degree2radial(-rotation/2));
    end
    
    
end