Plist={};

for i = 1:size(initialPWLModel.Triangles,1)
    for j=1:size(initialPWLModel.Triangles,2)
    Plist(initialPWLModel.Triangles(i,j),end+1)={initialPWLModel.Triangles(i,:)};
    end
end
Neighbour_points={};
for i=1:size(Plist,1)
    P_i=Plist(i,~cellfun('isempty',Plist(i,:)));
    P_i=unique(horzcat(P_i{:}));
    Neighbour_points(i)={P_i(find(P_i~=i))};
end
Pdist=[];
for i=1:size(Neighbour_points,2)
    Pdist(i)=mean(pdist2(initialPWLModel.Points(i,1:2),initialPWLModel.Points(Neighbour_points{i}(:),1:2)));
end

