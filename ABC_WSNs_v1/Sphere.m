function coveragarea = Sphere(x,Trange,area)
Nodes_pos = reshape(x,[numel(x)/2,2]);
%distribute pts points
pts=10000;
for i=0 : 99
    for j= 1 : 100
        xa(i*100+j) = i+1;
        ya(i*100+j) = j;
    end
end
pointspos=(area(1)/100).*[xa(:) ya(:)];
coveredpt=zeros(pts,1);
for ii=1:pts
    for jj=1:size(Nodes_pos,1)
        dist = sqrt((pointspos(ii,1)-Nodes_pos(jj,1))^2+(pointspos(ii,2)-Nodes_pos(jj,2))^2);
        if dist<Trange || dist==Trange
            coveredpt(ii)=1;
            break
        end
    end
end
%%
coveragarea=1/(numel(find(coveredpt==1))/pts);


