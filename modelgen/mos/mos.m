%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOSFET modeling for PRAISE-pwl-models %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate a reduced PWL model starting with a sampled DC-sweep output (CSV file)
% using the surface simplification method 'reducepatch()'.
% 
% Inputs:
% 	settings.m:		Some modeling information are needed.
%	saber_output.csv:	Sampled nonlinear device
%
% Outputs:
%	data/model.xml:		reduced PWL model
%	data/model.dat:		reduced PWL model (obsolete)
%
% Example:
%	-
%
% Needs:
%	programed and tested with Matlab R2012a
%
% Author: Lukas Lee <lukas.lee@ims.uni-hannover.de> 
%         Don't hesitate to contact me in case of problems
% Last update: 15-Jan-2016  
% Last revision: 1.1

if (exist('settings.m', 'file'))
	settings;
end

if ~(exist('element_name'))
	error('Variable "element_name" is not defined!')
end
if ~(ischar(element_name))
    error('Variable "element_name" is not a string!')
end

if ~(exist('reduce_factor'))
	error('Variable "reduce_factor" is not defined!')
end
if ~(isfloat(reduce_factor))
    error('Variable "reduce_factor" is not a number!')
end

if ~(exist('Vgs_min'))
	error('Variable "Vgs_min" is not defined!')
end
if ~(isfloat(Vgs_min))
    error('Variable "Vgs_min" is not a number!')
end

if ~(exist('Vgs_max'))
	error('Variable "Vgs_max" is not defined!')
end
if ~(isfloat(Vgs_max))
    error('Variable "Vgs_max" is not a number!')
end

mos_warning = fopen('MOS_warnings.log', 'w');
fprintf(mos_warning, '\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  set global variables       %
debug_mode                 = 1; %   0=off, 1=on
delete_points              = 0; %   0=off, 1=on
re_triangulation           = 0; %   0=off, 1=on
praise_initial             = 1; %   0=off, 1=on
del_perpendicular_triangl  = 1; %   0=off, 1=on
round_values               = 0; %   0=off, 1=on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTES %%%
% debug_mode:
% - call function check_model.m and write some debug information
% delete_points:
% - delete vertices of triangles which causes neg. Rds (neg. gradient in Uds)
% re_triangulation:
% - if some vertices are deleted or user has modifed the coordinates of
% vertices, the surface has to be remodeled
% praise_initial:
% - search for initial triangle, which includes the coordiante (0,0)

if (debug_mode == 1)
    mos_log = fopen('MOS_debugout.log', 'w');
    fprintf(mos_log, 'MOSFET modeling settings:\n',debug_mode);
    fprintf(mos_log, '=========================\n',debug_mode);
    fprintf(mos_log, '%g : debug_mode\n',debug_mode);
    fprintf(mos_log, '%g : delete_points\n',delete_points);
    fprintf(mos_log, '%g : re_triangulation\n',re_triangulation);
    fprintf(mos_log, '%g : praise_initial\n',praise_initial);
    fprintf(mos_log, '%g : del_perpendicular_triangl\n\n\n',del_perpendicular_triangl);
end;

%%% get characteristic data from Saber output %%%
% read Saber ouput
input=csvread('data/saber_output.csv',1,0);
% parsing range of V_GS
%Vgs_min = ...;
%Vgs_max = ...;
%Vgs_min=0;     %(for user DEBUG)
%Vgs_max=20;    %(for user DEBUG)
delta_Vgs = (Vgs_max-Vgs_min);
nr_of_steps = (length(input(1,:))-2);
Vgs_step = delta_Vgs / nr_of_steps;

for nr_points=1:length(input(1,:))-1
   x(:,nr_points) = input(:,1);
end;
% parsing range of V_DS
Vds_min = min(x(:,1));
Vds_max = max(x(:,1));
delta_Vds = Vds_max - Vds_min; % necessary for check_model.m
for nr_points=1:length(input(:,1))
    for step=0:nr_of_steps
        y(nr_points,step+1) = Vgs_min + (step*Vgs_step);
    end;
end;

for nr_points=1:length(input(:,1))
    for step=1:nr_of_steps+1
        z(nr_points,step) = input(nr_points,step+1);
    end;
end;

%%% use debug_dummy_mosfet_characteristic %%%
%x=[]; y=[]; z=[]; 
%x=[0,0,0;1,1,1;2,2,2;];
%y=[0,1,2;0,1,2;0,1,2;];
%z=[0,0,0;0,2,8;0,2,9;];

%%% triangulation (delaunay) %%%
tri = delaunay(x,y);
%%% get patch of characteristic %%%
p = trisurf(tri,x,y,z);
%%% output format for elaboration %%%
%p = surf(x,y,z);
set(p,'Visible','off');

if (debug_mode == 1)
    % hold on triangulated surface (org. saber input) for overlapping with reduced model
    hold on;
    set(p,'Visible','off');
end;

%%% reduce patch %%%
%%%reduce_factor = ...;
%reduce_factor = 0.50;     %(for user DEBUG)
if (debug_mode == 1)
    [F_org,P] = reducepatch(p,reduce_factor,'verbose');
    %[F_org,P] = reducepatch(p,reduce_factor,'fast');
    % NOTE: verbose = prints progress messages to command window as the computation progresses
    % NOTE: fast    = assumes the vertices are unique and does not compute shared vertices.
else
    [F_org,P] = reducepatch(p,reduce_factor);
end;

% %%% create 'MOS_model.csv' if not already exist %%%
% if (exist('MOS_model.csv') == 0)
%     %%% delete points on the right border (Vds_max), which cause neg. Rds %%%
%     if (delete_points == 1)
%         k=1;
%         del_points=0;
%         for i=1:length(P(:,1))
%             if ((P(i,1) == Vds_max) && (P(i,2) < Vgs_max) && (P(i,2) > Vgs_min))
%                 del_points(k)=i;
%                 k=k+1;
%             end;
%         end;
%         %delete column (point)
%         if (del_points ~= 0)
%             for n=length(del_points):-1:1
%                 P(del_points(n),:) = [];
%             end;
%         end;
%     end;
%     csvwrite('MOS_model.csv',P);
% % otherwise read existing 'MOS_model.csv'
% else
% 	fprintf(mos_warning, 'WARNING: generated model read from csv-file: "MOS_model.csv" \n');
%     P = csvread('MOS_model.csv');
% end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re_triangulation 
if (re_triangulation == 1)

    % set new coordinates(x,y,z)
    x_new = P(:,1)';
    y_new = P(:,2)';
    z_new = P(:,3)';

    % new triangulation
    tri_new = delaunay(x_new,y_new);
    p_new = trisurf(tri_new,x_new,y_new,z_new);
    set(p_new,'Visible','on');
    [F_new,P_new] = reducepatch(p_new,1);

    %set new faces(F)&vertices(P)
    F_org=F_new;
    P=P_new;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unusedVertexCounter=0;
%%% delete perpendicular triangles
if (del_perpendicular_triangl == 1)
    j=1;
    minimum_bbox_size = 0;  %criterion for minimal bounding box size of each triangle (set by user)
    facecounter = 1;
    %%% check geometry of each triangle; delete perpendicular triangles
    for (face=1:length(F_org(:,1)))
        triple = F_org(face,:);
        %get three vertices of current triangle
        for i=1:3
            X(i) = P(triple(i),1);
            Y(i) = P(triple(i),2);
            Z(i) = P(triple(i),3);
        end;
        %create bounding box of current triangle
        max_x = max(X)-min(X);
        max_y = max(Y)-min(Y);
        if ((max_x > minimum_bbox_size) && (max_y > minimum_bbox_size))
            %delete current face (row of F)
            F(facecounter,:) = F_org(face,:);
            facecounter = facecounter+1;
        else
            %deleted_faces(j)=face; %(Debugout)
            
            %identify unused vertex
            for iter1=1:3
                unusedVertex=0;
                currentVertex= F_org(face,iter1);
                for (iter2=1:length(F_org(:,1)))
                   % searching for using vertex
                   if ( ~isempty(find(F_org(iter2,:)==currentVertex)) )
                       unusedVertex=unusedVertex+1;
                   else
                   end
                end
                if (unusedVertex==1)
                    unusedVertex=0;
                    unusedVertexCounter=unusedVertexCounter+1;
                    VertexToDelete(unusedVertexCounter)=currentVertex;
                end
            end
            clear iter1 iter2
           
            fprintf(mos_warning, 'WARNING: delete perpendicular triangle with nr. %g\n',face-1);
            j=j+1;
        end;
    end;
else
    F=F_org;
end;

%delete vertex and correct F_org
for iterVertex=1:length(VertexToDelete)
    %delete vertex
    P(VertexToDelete(iterVertex),:)=[];
    fprintf(mos_warning, 'WARNING: delete vertex with nr. %g\n',VertexToDelete(iterVertex));

    %correction of F_org
    for (iter2=1:length(F(:,1)))
        for iter1=1:3
           if ( F(iter2,iter1) > VertexToDelete(iterVertex) )
               F(iter2,iter1)=F(iter2,iter1)-1;
           end
        end
    end
end





%%% set number of faces/planes %%%
nr_of_faces = length(F(:,1));

%%% check if 'user input' == nr_of_faces %%%
if (reduce_factor ~= nr_of_faces && reduce_factor ~= 1)
    fprintf(mos_warning, 'WARNING: nr. of triangles has been changed to %g \n',nr_of_faces);
end;

%%% show reduced graphic %%%
p=patch('Faces',F,'Vertices',P);
%set(p,'FaceColor','w','EdgeColor','b','LineWidth',0.5)

grid on
set(gca,'LineWidth',1)
xlim([min(P(:,1)) max(P(:,1))])
xlabel('U_{DS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
ylim([min(P(:,2)) max(P(:,2))])
ylabel('U_{GS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
zlim([min(P(:,3)) max(P(:,3))])
zlabel('I_{DS}(A)', 'FontSize', 18, 'FontName', 'Helvetica');
p=patch('Faces',F,'Vertices',P);
view(-35, 20);
set(p,'FaceColor',[0.8,0.8,0.8],'EdgeColor','black','LineWidth',1)

%%% fill neighbour matrix with initial value "-1" %%%
neighbour=zeros(nr_of_faces,3)-1;

%%% get all coordinates for current triangle %%%
for (face=1:nr_of_faces)
    %triple = sort(F(face,:));
    triple = (F(face,:));
    for i=1:3
        X(i) = P(triple(i),1);
        Y(i) = P(triple(i),2);
        Z(i) = P(triple(i),3);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% get neighbour triangles

    % append first vertex for pairwise operation
    triple(1,4)=triple(1,1);
    counter = 0;
    
    % for each vertices of current triangle
    for pair=1:3
        % get the first pair of vertices of current triangle
        current_nodes(1) = triple(1,pair);
        current_nodes(2) = triple(1,pair+1);
        % search current pair in other triangles
        for i=1:length(F(:,1))
            % skip current triangle (own triangle)
            if (i~=face)
                % check if current pair found in another triangle
                for j=1:3
                    if (F(i,j)==current_nodes(1)) || (F(i,j)==current_nodes(2))
                        counter=counter+1;
                    end;
                end;
                % if current pair found in another triangle => neighbour found
                if (counter == 2) && (neighbour(face,pair)==-1)
                    neighbour(face,pair) = i-1;
                end;
                counter=0;
            end;
        end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% - computation of gradient (m) in 2 dim. (Uds,Ugs)
    %%% - define upper/lower limit for triangle egde
    
        % append a vertex for pair of vertices (P1&P2, P2&P3, P3&P1)
        X(4)=X(1);
        Y(4)=Y(1);

        % computation of centroid of the current triangle
        sx=(X(1)+X(2)+X(3))/3;
        sy=(Y(1)+Y(2)+Y(3))/3;

        % define 2d rotation matrix
        rot_matrix = [0,-1;
                      1,0;];

        % sort X coordinates of current vertex (increasing from left->right)
        if (X(pair)>X(pair+1))
            P_left  = [X(pair+1);Y(pair+1);];
            P_right = [X(pair);Y(pair);];
        elseif (X(pair)<X(pair+1))
            P_left  = [X(pair);Y(pair);];
            P_right = [X(pair+1);Y(pair+1);];
        else
            % if X coordinates are equal, sort by Y coordinates (increasing from bottom up)
            if (Y(pair) < Y(pair+1))
                P_left  = [X(pair);Y(pair);];
                P_right = [X(pair+1);Y(pair+1);];
            else
                P_left  = [X(pair+1);Y(pair+1);];
                P_right = [X(pair);Y(pair);];
            end;
        end;
        current_vektor = P_right-P_left;
        %NOTE: x_coord=current_vektor(1); y_coord=current_vektor(2)

        % calculating centroid of each triangle egde
        sx_edge = (P_left(1)+P_right(1))/2;
        sy_edge = (P_left(2)+P_right(2))/2;

        %special case: if m=Inf
        if (current_vektor(1) == 0)
            %set gradient 
            m(face,pair) = Inf;
            %set border
            if (sx_edge < sx)
                limit(face,pair) = 0;
            else
                limit(face,pair) = 1;
            end;
        %special case: if m=0
        elseif (current_vektor(2) == 0)
            %set gradient
            m(face,pair) = 0;
            %set border
            if (sy_edge < sy)
                limit(face,pair) = 0;
            else
                limit(face,pair) = 1;
            end;
        else
            %set gradient
            m(face,pair) = (Y(pair+1)-Y(pair)) / (X(pair+1)-X(pair));
            
            %set border (rotate current gradient, 90 degree anti-clockwise)
            current_vektor_ortho = rot_matrix*current_vektor;
            %define coordiante (x,y) for checkpoint for "inpolygon"
            drift = 0.00001; % 0.001%
            in_check_point(1) = sx_edge+drift*current_vektor_ortho(1);
            in_check_point(2) = sy_edge+drift*current_vektor_ortho(2);

            in = inpolygon(in_check_point(1),in_check_point(2),X,Y);

            % coordiante is outside of triangle
            if (in == 0)
                limit(face,pair) = 1;
            % coordiante is inside of triangle
            else
                limit(face,pair) = 0;
            end;
        end;
    end;
            %check if the coordiante (0,0) within triangle for "initial-PRAISE-triangle"
            initial_praise_triangle(face) = inpolygon(0,0,X,Y);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% computation of y-intercept in 2 dim. (Uds,Ugs)                 %%%

    % y-intercept
    if (m(face,1) == Inf)
        b(face,1) = X(1);
    else
        b(face,1) = Y(1)-(m(face,1)*X(1));
    end;
    if (m(face,2) == Inf)
        b(face,2) = X(2);
    else
        b(face,2) = Y(2)-(m(face,2)*X(2));
    end;
    if (m(face,3) == Inf)
        b(face,3) = X(3);
    else
        b(face,3) = Y(3)-(m(face,3)*X(3));
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% computation of gradient and z-intercept in 3 dim. (Uds,Ugs,I)  %%%
    % paper: 2dim pwl (Jimenez-Fernandez) % 

    %page(322) / (7)
    K=[X(1),Z(1),Y(1);
       X(2),Z(2),Y(2);
       X(3),Z(3),Y(3);];
    delta_ny = det(K);
    
    %page(322) / (8)
    L=[Y(1),X(1),1;
       Y(2),X(2),1;
       Y(3),X(3),1;];
    delta = det(L);

    %page(322) / (6) (z-intercept)
    G(face) = delta_ny/delta;

    %page(322) / (4) (gradient of Ugs)
    A(face) = ( Z(1)*X(2) - Z(1)*X(3) - X(1)*Z(2) + X(1)*Z(3) + Z(2)*X(3) - X(2)*Z(3) ) / delta;

    %page(322) / (5) (gradient of Uds)
    B(face) = ( Y(1)*Z(2) - Y(1)*Z(3) - Z(1)*Y(2) + Z(1)*Y(3) + Y(2)*Z(3) - Z(2)*Y(3) ) / delta;
    
%     %check for neg. Rds
%     for i=length(B(1,:))
%         if (B(i) < 0)
%             negative_rds(i) = 1;
%         else
%             negative_rds(i) = 0;
%         end;
%     end;
%     flag_neg_rds = sum(negative_rds);
end;

% %end MATLAB script if negative Rds is detected
% if (flag_neg_rds > 0)
%     fprintf(mos_warning, 'ERROR: neg. Rds\n');
%     fprintf(mos_warning, 'ERROR: file ''model.dat'' and ''model.xml'' NOT created!\n');
%     fclose(mos_log);
%     fclose(mos_warning);
%     return;
% end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% search for PRAISE initial triangle "0" %%%
%fill default output_sequence (beginning with 1)
for i=1:nr_of_faces
    output_sequence(i)= i;
end;
initial_triangle = 0;
if (praise_initial == 1)
    maxA=1e12;  %set high value for default
	for i=1:length(F(:,1))
        if (initial_praise_triangle(i) == 1)
            %search for smallest vccs (A)
            if (abs(A(i)) < maxA)
                maxA = abs(A(i));
                initial_triangle = i;
            end;
        end;
	end;

    %swapping initial triangle with first triangle
    if (initial_triangle > 1)
        output_sequence(1,initial_triangle)=1;
        output_sequence(1,1)=initial_triangle;

	Fbackup1 = F(1,1);
	Fbackup2 = F(1,2);
	Fbackup3 = F(1,3);
	F(1,1) = F(initial_triangle,1);
	F(1,2) = F(initial_triangle,2);
	F(1,3) = F(initial_triangle,3);
	F(initial_triangle,1) = Fbackup1;
	F(initial_triangle,2) = Fbackup2;
	F(initial_triangle,3) = Fbackup3;

        %adapating neighbours because of swapped triangles
        for i=1:length(neighbour(:,1))
            for j=1:3
                if (neighbour(i,j) == initial_triangle-1)
                    neighbour(i,j) = 0;
                elseif (neighbour(i,j) == 0)
                    neighbour(i,j) = initial_triangle-1;
                end;
            end;
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Debug_out:
%F,P,A,B,G,m,b,limit,neighbour,initial_triangle,output_sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% round values
if (round_values == 1)
    for i=1:length(A)
        A(i)=eval(sprintf('%.12f',A(i)));
        B(i)=eval(sprintf('%.12f',B(i)));
        G(i)=eval(sprintf('%.12f',G(i)));
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create outputfile
fid = fopen('model.dat', 'w');
    % HEADER
    fprintf(fid, 'MODULENAME:\n');
    fprintf(fid, '%s\n', element_name);
    fprintf(fid, 'NODES:\n');
    fprintf(fid, 's d g\n');
    fprintf(fid, 'ELEMENTS:3\n');
    fprintf(fid, 'r.ID p:d m:s = rnom=ID\n');
    fprintf(fid, 'i.ID p:d m:s = dc=ID, ac=ID\n');
    fprintf(fid, 'vccs.ID vp:g vm:s p:d m:s = k=ID\n');

    % INPUT1
    fprintf(fid, 'INPUT1:\n');
    fprintf(fid, 'voltage d s %d [',nr_of_faces);
    for i=1:nr_of_faces
            for j=1:3
                fprintf(fid, '(%g ', limit(output_sequence(i),j));
                fprintf(fid, '%g ', m(output_sequence(i),j));
                fprintf(fid, '%g ', b(output_sequence(i),j));
                fprintf(fid, '%g)', neighbour(output_sequence(i),j));
                if (j==3)
                  if (i~=nr_of_faces)
                    fprintf(fid, ';');
                  end;
                else
                    fprintf(fid, ',');
                end;
            end;
            if (i==nr_of_faces)
                fprintf(fid, ']\n');
            end;
    end;

    % INPUT2
    fprintf(fid, 'INPUT2:\n');
    fprintf(fid, 'voltage g s\n');
    
    % PWLMODEL
    fprintf(fid, 'PWLMODEL:%d\n', nr_of_faces);
    fprintf(fid, '[');
    
    for i=1:nr_of_faces
%         if ((1/B(output_sequence(i))) < 0)
%             fprintf(mos_warning, 'ERROR: neg. Rds in triangle %g\n',i-1);
%             %NOTE: "-1" because beginning numbering of faces in PRAISE at "0"
%             %NOTE: in MATLAB the nr of faces/triangles beginning with "1" 
%         end;
        fprintf(fid, '(%g ', 1/B(output_sequence(i)));
        fprintf(fid, '%g ', G(output_sequence(i)));
        fprintf(fid, '%g)', A(output_sequence(i)));
        if (i==nr_of_faces)
            fprintf(fid, ']\n');
        else
            fprintf(fid, ',');
        end;
    end;
fclose(fid);

fxml = fopen('model.xml', 'w');
fsaber = fopen('model.ai_dat', 'w');
	fprintf(fxml, '<praise>\n');
	fprintf(fxml, '<pwlmodels>\n');
	fprintf(fxml, '<element name="%s" variant="VARIANT_NAME">\n', element_name);
	fprintf(fxml, '	<model2dim name="ids">\n');

% 	fprintf(fxml, '	<nodes>\n');
% 	fprintf(fxml, '		<node name="s"/>\n');
% 	fprintf(fxml, '		<node name="d"/>\n');
% 	fprintf(fxml, '		<node name="g"/>\n');
% 	fprintf(fxml, '	</nodes>\n');

% 	fprintf(fxml, '	<inputs>\n');
% 	fprintf(fxml, '		<input1 type="voltage" node1="d" node2="s"/>\n');
% 	fprintf(fxml, '		<input2 type="voltage" node1="g" node2="s"/>\n');
% 	fprintf(fxml, '	</inputs>\n');

	% new XML format
 	fprintf(fxml, '		<input1 type="voltage" node1="d" node2="s"/>\n');
 	fprintf(fxml, '		<input2 type="voltage" node1="g" node2="s"/>\n');
 	fprintf(fxml, '		<output type="current" node1="d" node2="s"/>\n');

% 	fprintf(fxml, '	<elements>\n');
% 	fprintf(fxml, '		<element id="0" name="r.ID p:d m:s = rnom=ID"/>\n');
% 	fprintf(fxml, '		<element id="1" name="i.ID p:d m:s = dc=ID, ac=ID"/>\n');
% 	fprintf(fxml, '		<element id="2" name="vccs.ID vp:g vm:s p:d m:s = k=ID"/>\n');
% 	fprintf(fxml, '	</elements>\n');

	fprintf(fxml, '		<vertices>\n');
	for i=1:size(P,1)
		x1 = P(i,1);
		x2 = P(i,2);
		y  = P(i,3);
		fprintf(fxml, '			<vertex id="%d" x1="%g" x2="%g" y="%g"/>\n', (i - 1), x1, x2, y);
        fprintf(fsaber, '%g %g %g\n', x1, x2, y);
	end;
	fprintf(fxml, '		</vertices>\n');

	fprintf(fxml, '		<triangles>\n');
	for i=1:size(F,1)
		vertex0 = F(i,1) - 1;
		vertex1 = F(i,2) - 1;
		vertex2 = F(i,3) - 1;
		fprintf(fxml, '			<triangle id="%d" vertex0="%d" vertex1="%d" vertex2="%d"/>\n', (i - 1), vertex0, vertex1, vertex2);
	end;
	fprintf(fxml, '		</triangles>\n');

    fprintf(fxml, '	</model2dim>\n');
	fprintf(fxml, '</element>\n');
	fprintf(fxml, '</pwlmodels>\n');
	fprintf(fxml, '</praise>\n');
fclose(fxml);
fclose(fsaber);

fprintf(mos_warning, 'Please set VARIANT_NAME manually in outfile "data/model.xml"!\n');
fprintf(mos_warning, '\n');

%save as .png image 
xlabel('Uds')
ylabel('Ugs')
zlabel('I (Uds,Ugs)')
title('MOS characteristic')
print -dpng MOS_characteristic.png;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 2dim-Layer (Uds, Ugs)
% note: plot triangles by coordinates

if (debug_mode == 1)
    figure;
    fprintf(mos_log, 'DEBUGOUT: 2dim-coordinates of vertices\n');
    fprintf(mos_log, '\tX\tY\t\n');
    fprintf(mos_log, '===============\n');

    for (face=1:nr_of_faces)
        %%% free X,Y
        X=[];
        Y=[];

        triple = F(face,:);
        for i=1:3
            X(i) = P(triple(i),1);
            Y(i) = P(triple(i),2);
        end;

        X(4)=X(1);
        Y(4)=Y(1);
        plot(X,Y,'Color','red','LineWidth',0.5);hold on;

        %%% DEBUGOUT
        if (debug_mode == 1)
            fprintf(mos_log, '\n%g',face-1);
            fprintf(mos_log, '\t%g\t%g\t\n',X(1),Y(1));
            fprintf(mos_log, '\t%g\t%g\t\n',X(2),Y(2));
            fprintf(mos_log, '\t%g\t%g\t\n',X(3),Y(3));
        end;
    end;

    % check 2dim-Layer (Uds, Ugs)
%    check_model(debug_mode, delta_Vds, delta_Vgs, F, P);
    % NOTE: after "check_model.m" all lines have to be blue
    % if red lines of triangles are still visible, the written model.dat is wrong
    fclose(mos_log);
end;
fclose(mos_warning);

