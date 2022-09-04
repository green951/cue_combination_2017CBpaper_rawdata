
dist1=3.4;
dist2=4.5;

% pos1=[[130,dist1],[70,dist2],[40,dist1]]
% pos2=[[150,dist1],[100,dist2+.5],[30,dist1]]
% pos3=[[110,dist1],[60,dist2],[20,dist2]]
% pos4=[[165,dist2],[105,dist2+.5],[50,dist1]]
% pos5=[[135,4],[90,4],[45,4]]
% pos10=[[145,4.2],[90,2.5],[35,4.2]]

pos(1,:,1)=[130,dist1];
pos(2,:,1)=[70,dist2];
pos(3,:,1)=[40,dist1];

pos(1,:,2)=[150,dist1];
pos(2,:,2)=[100,dist2+.5];
pos(3,:,2)=[30,dist1];

pos(1,:,3)=[110,dist1];
pos(2,:,3)=[60,dist2];
pos(3,:,3)=[20,dist2];

pos(1,:,4)=[165,dist2];
pos(2,:,4)=[105,dist2+.5];
pos(3,:,4)=[50,dist1];

pos(1,:,5)=[135,4];
pos(2,:,5)=[90,4];
pos(3,:,5)=[45,4];

pos(1,:,10)=[145,4.2];
pos(2,:,10)=[90,2.5];
pos(3,:,10)=[35,4.2];

num1=[1 2 3 4];
num2=[6 7 8 9];
for i=1:3
    for j=1:4
    pos(i,1,num2(j))=180-pos(i,1,num1(j));
    pos(i,2,num2(j))=pos(i,2,num1(j));
    end
end

for n=1:10
    for i=1:3
        x=pos(i,2,n)*cos(pos(i,1,n)*pi/180);
        z=pos(i,2,n)*sin(pos(i,1,n)*pi/180);
        pos2(i,:,n)=[x z];
    end
end

figure(1)        
for n=1:10
    subplot(2,5,n)
    plot(pos2(:,1,n),pos2(:,2,n),'s-')
    hold on
    axis([-5,5,1,5])
    title(['landmarks positions:' num2str(n-1)])
end
    
    
    
    

