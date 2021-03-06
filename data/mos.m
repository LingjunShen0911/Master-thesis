%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOSFET modeling for PRAISE PWL models %%%
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
%
% Example:
%	-
%
% Needs:
%	programed and tested with Matlab R2014b
%
% Author: Lukas Lee <lukas.lee@ims.uni-hannover.de> 
%         Don't hesitate to contact me in case of problems
% Last update: 29-Jan-2016  
% Last revision: 1.2

clear;

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
praise_initial             = 1; %   0=off, 1=on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTES %%%
% debug_mode:
% - call function check_model.m and write some debug information
% vertices, the surface has to be remodeled
% praise_initial:
% - search for initial triangle, which includes the coordiante (0,0)

if (debug_mode == 1)
    mos_log = fopen('MOS_debugout.log', 'w');
    fprintf(mos_log, 'MOSFET modeling settings:\n',debug_mode);
    fprintf(mos_log, '=========================\n',debug_mode);
    fprintf(mos_log, '%g : debug_mode\n',debug_mode);
    fprintf(mos_log, '%g : praise_initial\n',praise_initial);
end;

%%% get characteristic data from Saber output %%%
% read Saber ouput
input=csvread('saber_output.csv',1,0); %%csvread:从saber_output文件中读取第1行往下第0列往后的数据集。just read the number not the words

%%input是以步长0.68 equidistance sampling

% parsing range of V_GS
%Vgs_min = ...;
%Vgs_max = ...;
%Vgs_min=0;     %(for user DEBUG)
%Vgs_max=20;    %(for user DEBUG)
delta_Vgs = (Vgs_max-Vgs_min);
nr_of_steps = (length(input(1,:))-2);   %%列的长度 27-2=25
Vgs_step = delta_Vgs / nr_of_steps;     %%将Vgs均分 12/25

for nr_points=1:length(input(1,:))-1    %%1:27=27 27-1
   x(:,nr_points) = input(:,1);         %% input(:,1):第一列的所有数值，x的值就是就是将前述的第一列重复26遍.x行一致

end;

% parsing range of V_DS
Vds_min = min(x(:,1));  %% -2
Vds_max = max(x(:,1));     %% 15
delta_Vds = Vds_max - Vds_min; %% necessary for check_model.m =17

for nr_points=1:length(input(:,1))  %% 1:27
    for step=0:nr_of_steps      %% 0:25
        y(nr_points,step+1) = Vgs_min + (step*Vgs_step); %% at first y(1,1)= -2+(1*0.4)以步长0.4前进，将y轴也均等分为26份。y列一致
    end;
end;

%%parsing range of I_ds

for nr_points=1:length(input(:,1)) 
    for step=1:nr_of_steps+1        %% 1:26
        z(nr_points,step) = input(nr_points,step+1);   %% z的26*26= input(1:26,2:27)因为input(1:26,1)是Vds的值
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

if (debug_mode == 1)
    figure;
    p = trisurf(tri,x,y,z,'Visible','on');
    grid on
    xlabel('Uds')
    ylabel('Ugs')
    zlabel('I (Uds,Ugs)')
    title('Sampled model')
    az = -35; %% 负y轴绕z轴-35度
    el = 35;  
    view(az, el); %%view为视点指定
    %save as .png image 
    print -dpng MOS_sampled.png;
else
    p = trisurf(tri,x,y,z,'Visible','off');  %%将三角图画完后返回一个填充句柄
end


%%% reduce patch %%%
%reduce_factor = 18;     %(for user DEBUG)

%% if reducefactor大于1，表示三角形的目标数，如果小于1，为减为原来的百分之r数。

if (debug_mode == 1)
    [F,P] = reducepatch(p,reduce_factor,'verbose');  %% 函数表义为缩减句柄。p为原三角形数,r为减少到某个数
    %[F_org,P] = reducepatch(p,reduce_factor,'fast');
    % NOTE: verbose = prints progress messages to command window as the computation progresses
    %%在计算过程中将进度消息输出至命令行窗
    % NOTE: fast    = assumes the vertices are unique and does not compute shared vertices.
else
    [F,P] = reducepatch(p,reduce_factor); %%F为返回面，P为返回顶点
end;

%%% set number of faces/triangles %%%
nr_of_faces = length(F(:,1));  %%面的数量
 
%%% check if 'user input' == nr_of_faces %%%
if (reduce_factor ~= nr_of_faces)
    fprintf(mos_warning, 'WARNING: nr. of triangles has been changed to %g \n',nr_of_faces);
end;

if (debug_mode == 1)
    %%% show reduced graphic %%%    %%画出r=输入值的简化模型
    figure;
    p=patch('Faces',F,'Vertices',P); %%创建多边形
    set(p,'FaceColor','w','EdgeColor','b','LineWidth',0.5)
    grid on
    xlabel('Uds')
    ylabel('Ugs')
    zlabel('I (Uds,Ugs)')
    title('Reduced model')
    az = -35;
    el = 35;
    view(az, el);
    %save as .png image 
    print -dpng MOS_reduced.png;
end

%%% write XML ouput file %%%
fxml = fopen('model.xml', 'w');    %%打开文件
fsaber = fopen('model.ai_dat', 'w');
	fprintf(fxml, '<praise>\n');    %%将数据写入文本文件
	fprintf(fxml, '<pwlmodels>\n');
	fprintf(fxml, '<element name="%s" variant="VARIANT_NAME">\n', element_name);
	fprintf(fxml, '	<model2dim name="ids">\n');

 	fprintf(fxml, '		<input1 type="voltage" node1="d" node2="s"/>\n');
 	fprintf(fxml, '		<input2 type="voltage" node1="g" node2="s"/>\n');
 	fprintf(fxml, '		<output type="current" node1="d" node2="s"/>\n');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 2dim-Layer (Uds, Ugs)
% note: plot triangles by coordinates

if (debug_mode == 1)
    figure;
    fprintf(mos_log, 'DEBUGOUT: 2dim-coordinates of vertices\n');
    fprintf(mos_log, '\tX\tY\t\n');
    fprintf(mos_log, '===============\n');

    for (face=1:length(F(:,1)))
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
    grid off
    xlabel('Uds')
    ylabel('Ugs')
    zlabel('I (Uds,Ugs)')
    title('Check: 2D projection of reduced model')
    az = 0;
    el = 90;
    view(az, el);
    %save as .png image 
    print -dpng MOS_2D_Projection.png;
    
    % check 2dim-Layer (Uds, Ugs)
%(obsolete) %        check_model(debug_mode, delta_Vds, delta_Vgs, F, P); 
%     % NOTE: after "check_model.m" all lines have to be blue
%     % if red lines of triangles are still visible, the written model.dat is wrong
    fclose(mos_log);
end;
fclose(mos_warning);

