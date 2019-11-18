function [] = check_model(debug_mode, delta_Vds, delta_Vgs, F, P)


if (debug_mode == 1)
    mos_check = fopen('MOS_check_model.log', 'w');
end;
counter=0;
n = 0;
m = 0;
loop_counter = 0;
triangle_of_nr=0;

%%% parsing model.dat %% can't find this file
fid=fopen('model.dat');
while(1)
    newline=fgetl(fid); %% 读取文件中的行，并删除文件中的换行符
    if (newline == -1)
        break;
    end
    %get 'INPUT1 from model.dat (2D)'
    if regexp(newline, 'voltage d s ')  %%匹配正则表达式，所有匹配项的开始和结束索引.在newline里搜索voltage d s
        [start_idx, end_idx, extents, matches, tokens, names, splits]= regexp(newline, ' ');
        %nr_of_faces = splits{4};
        for i=1:length(splits)
            if regexp(splits{i}, '(')
                n = n+1;
                %2D
                mydata(n,1) = str2num(regexprep(splits{i+0}, '.*\((.*)', '$1'));
                mydata(n,2) = str2num(splits{i+1});
                mydata(n,3) = str2num(splits{i+2});
                mydata(n,4) = str2num(regexprep(splits{i+3}, '(.*)\).*', '$1'));
            end;
        end;
    end;
    %get 'INPUT2 from model.dat (2D)'
    if regexp(newline, '^[(')
        [start_idx, end_idx, extents, matches, tokens, names, splits]= regexp(newline, ' ');
        %nr_of_faces = splits{4};
        for i=1:length(splits)
            if regexp(splits{i}, '(')
                m = m+1;
                %3D
                threedim_data(m,1) = str2num(regexprep(splits{i+0}, '.*\((.*)', '$1'));
                threedim_data(m,2) = str2num(splits{i+1});
                threedim_data(m,3) = str2num(regexprep(splits{i+2}, '(.*)\).*', '$1'));
            end;
        end;
    end;
end;
%clear end_idx extents fid i matches n names newline splits start_idx tokens ;

%%%DEBUGOUT
if (debug_mode == 1)
    fprintf(mos_check, 'DEBUGOUT: coordinates of vertices\n');
    fprintf(mos_check, '\tX\tY\t\n');
    fprintf(mos_check, '===============\n');
end;
%%% computating intersections of triangle
triangle_counter=0;
for data=1:3:length(mydata)
triangle_counter=triangle_counter+1;
    %%% free x,y
    x=[];
    y=[];
    
	%%% get gradient and y-intercept
    for k=1:3
        m(k) = mydata(data+(k-1),2);
        b(k) = mydata(data+(k-1),3);
    end;
    
    %%% computating 1. intersection of triangle
    % between line m1 & m2
    if (m(1) == Inf) && (m(2) ~= Inf)
        x(1) = mydata(data,3);
        y(1) = m(2)*x(1)+b(2);
    elseif (m(2) == Inf) && (m(1) ~= Inf)
        x(1) = mydata(data+1,3);
        y(1) = m(1)*x(1)+b(1);
    else
        x(1) = (b(2)-b(1)) / (m(1)-m(2));
        y(1) = m(1)*x(1)+b(1);
    end;
    
    %%% computating 2. intersection of triangle
	% between line m2 & m3
    if (m(2) == Inf) && (m(3) ~= Inf)
        x(2) = mydata(data+1,3);
        y(2) = m(3)*x(2)+b(3);
    elseif (m(3) == Inf) && (m(2) ~= Inf)
        x(2) = mydata(data+2,3);
        y(2) = m(2)*x(2)+b(2);
    else
        x(2) = (b(3)-b(2)) / (m(2)-m(3));
        y(2) = m(2)*x(2)+b(2);
    end;
    
    %%% computating 3. intersection of triangle
    % between line m3 & m1
    if (m(3) == Inf) && (m(1) ~= Inf)
        x(3) = mydata(data+2,3);
        y(3) = m(1)*x(3)+b(1);
    elseif (m(1) == Inf) && (m(3) ~= Inf)
        x(3) = mydata(data,3);
        y(3) = m(3)*x(3)+b(3);
    else
        x(3) = (b(1)-b(3)) / (m(3)-m(1));
        y(3) = m(3)*x(3)+b(3);
    end;

    %%% computation of centroid of the current triangle for numbering
    % for numbering of triangle
    sx=sum(x)/3;
    sy=sum(y)/3;

    %%% computation of centroid of the current triangle-edge for numbering
    sx_edge(1) = (x(1)+x(2))/2;
    sy_edge(1) = (y(1)+y(2))/2;
    sx_edge(2) = (x(2)+x(3))/2;
    sy_edge(2) = (y(2)+y(3))/2;
    sx_edge(3) = (x(3)+x(1))/2;
    sy_edge(3) = (y(3)+y(1))/2;
    
    %%% plot triangle
    % copy first point to plot triangle
    x(4) = x(1);
    y(4) = y(1);
    % plot
    plot(x,y,'Color','blue','LineWidth',0.5);hold on;

    %get coordiantes (x,y,z) of vertices of triangles
    counter=counter+1;
    x_check(counter,:) = x;
    y_check(counter,:) = y;
    z_check(counter,:) = (1/(threedim_data(counter,1)))*x + (threedim_data(counter,2)) + (threedim_data(counter,3))*y;

    
    if (debug_mode == 1)
        fprintf(mos_check, '\n%g',triangle_of_nr);
        fprintf(mos_check, '\t%g\t%g\t\n',x(1),y(1));
        fprintf(mos_check, '\t%g\t%g\t\n',x(2),y(2));
        fprintf(mos_check, '\t%g\t%g\t\n',x(3),y(3));
    end;

    %%% numbering of triangle (beginning at zero)
    text(sx,sy,sprintf('%d',triangle_of_nr),'FontSize',14,'FontWeight','bold');
    triangle_of_nr = triangle_of_nr+1;
    
    %%% numbering of edge (upper or lower limit)
    space_size_x = delta_Vds*.005;
    space_size_y = delta_Vgs*.02;
    if (mydata(data,2) == Inf)
        if (mydata(data+1,1) == 1)
            text(sx_edge(3)-space_size_x,sy_edge(3),sprintf('%d',1));
        else
            text(sx_edge(3)+space_size_x,sy_edge(3),sprintf('%d',0));
        end;
    else
        if (mydata(data,1) == 1)
            text(sx_edge(3),sy_edge(3)-space_size_y,sprintf('%d',1));
        else
            text(sx_edge(3),sy_edge(3)+space_size_y,sprintf('%d',0));
        end;
    end;
    
    if (mydata(data+1,2) == Inf)
        if (mydata(data+1,1) == 1)
            text(sx_edge(1)-space_size_x,sy_edge(1),sprintf('%d',1));
        else
            text(sx_edge(1)+space_size_x,sy_edge(1),sprintf('%d',0));
        end;
    else
        if (mydata(data+1,1) == 1)
            text(sx_edge(1),sy_edge(1)-space_size_y,sprintf('%d',1));
        else
            text(sx_edge(1),sy_edge(1)+space_size_y,sprintf('%d',0));
        end;
    end;
    
    if (mydata(data+2,2) == Inf)
        if (mydata(data+2,1) == 1)
            text(sx_edge(2)-space_size_x,sy_edge(2),sprintf('%d',1));
        else
            text(sx_edge(2)+space_size_x,sy_edge(2),sprintf('%d',0));
        end;
    else
        if (mydata(data+2,1) == 1)
            text(sx_edge(2),sy_edge(2)-space_size_y,sprintf('%d',1));
        else
            text(sx_edge(2),sy_edge(2)+space_size_y,sprintf('%d',0));
        end;
    end;
    
    %%% TODO %%%
    %set array of cartesian coordiantes (x,y) of vertices of triangles
    %NOTE: needed by function polybool to catch overlapping triangles
    %check_overlap_x(triangle_counter,:)=x;
    %check_overlap_y(triangle_counter,:)=y;
    
end;

    %%% check overlapping triangles (::TODO::)
    %check_overlap_x,check_overlap_y; % (DEBUGOUT)
    % function 'polybool'

if (debug_mode == 1)
    fclose(mos_check);
end;

%save result as .png image 
xlabel('Uds')
ylabel('Ugs')
title('Check-Model: Uds-Ugs triangles in 2dim-layer')
print -dpng check_model.png;

%delete redundant column
x_check(:,4)=[];
y_check(:,4)=[];
z_check(:,4)=[];

%plot recomputated
figure();
stem3 (x_check, y_check, z_check, 'DisplayName', 'x_check, y_check, z_check'); hold on;

tri_check = delaunay(x_check,y_check);
p_check = trisurf(tri_check,x_check,y_check,z_check);
set(p_check,'Visible','off');

%write csv file
[F_check,P_check] = reducepatch(p_check,1);
csvwrite('MOS_model_check.csv',P);

%figure;
p=patch('Faces',F,'Vertices',P);
set(p,'FaceColor','w','EdgeColor','b','LineWidth',0.5)

xlabel('Uds')
ylabel('Ugs')
zlabel('I (Uds,Ugs)')
title('MOS characteristic (reading by model.dat)')
print -dpng MOS_characteristic_check.png;

end
