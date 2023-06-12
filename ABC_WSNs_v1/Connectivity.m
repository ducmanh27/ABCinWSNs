function check = Connectivity(pop,numberNodes,Rcom)
pop1 = pop(1:numberNodes);
pop2 = pop((numberNodes +1):(2*numberNodes));
pop = [pop1(:),pop2(:)];
adj_matrix = zeros(numberNodes,numberNodes);
for i=1:numberNodes
    for j=1:numberNodes
        if (((pop(i,1)-pop(j,1))^2+(pop(i,2)-pop(j,2))^2)<=(Rcom)^2)
            adj_matrix(i,j) = 1;
        end
    end
end
for i=1:numberNodes
    adj_matrix(i,i) = 0;
end
G= graph(adj_matrix);
v = dfsearch(G,1);
if (size(v,1)==numberNodes)
    check = 1;
else
    check = 0;
end