function [vertex, newRevertCount] = move_Point(PWLModel, relativeSchrittweite, weight, vertexIndex, Interpolated_Ref, Typ, Modellgrenzen)
%MOVE_POINT Random motion of a selected vertex
% move
% Inputs:
% p     -   Reduzierungsfaktor der Range
% temp  -   Kopie des aktuellen Punkts muss von außen kommen, damit sich
%           der Vektor nicht selber überschreibtmove
% 
% Author: Lukas Lee, Gerald Koroa
% Last Edit: 28.10.16


    Xmax=Modellgrenzen(1);
    Xmin=Modellgrenzen(2);
    Xstep=Modellgrenzen(3);
    Ymax=Modellgrenzen(4);
    Ymin=Modellgrenzen(5);
    Ystep=Modellgrenzen(6);
    Zmax=Modellgrenzen(7);
    Zmin=Modellgrenzen(8);

    % keep old position of selected vertex
    oldVertex=PWLModel.Points(vertexIndex,:);
if (strcmp(Typ,'MOD'))
    %------------------------------------------------X
    if (or(PWLModel.Index(vertexIndex)==0,PWLModel.Index(vertexIndex)==2))
        while 1
            % Modified with weight version
            Richtung = 1;
            Zufall = rand(1);
            if Zufall < weight(1)
                Richtung = -1;    
            end
            Zufall = rand(1);
            add=relativeSchrittweite.X*Zufall*Richtung;

            % move selected vertex
            vertex(1)=oldVertex(1)+add;
            % update PWL model
            PWLModel.Points(vertexIndex,1)=vertex(1);
 
        %     %FOR DEBUG
        %     surf(x,y,z)

            % if vertex movement exceed limits
            if (or(vertex(1)<Xmin,vertex(1)>Xmax));
                % undo moving vertex
                PWLModel.Points(vertexIndex,1)=oldVertex(1);
                vertex(1)=oldVertex(1);
                PWLModel.RevertCount=PWLModel.RevertCount+1;
                disp('WARNING: revert vertex move on x axis!')
            else
                break
            end
        end
    else
        vertex(1)=oldVertex(1);
    end
    
    %-----------------------------------------------Y
    if (or(PWLModel.Index(vertexIndex)==0,PWLModel.Index(vertexIndex)==1))
        while 1
            % Modified with weight version
            Richtung = 1;
            Zufall = rand(1);
            if Zufall < weight(4)
                Richtung = -1;    
            end
            Zufall = rand(1);
            add=relativeSchrittweite.Y*Zufall*Richtung;

            % move selected vertex
            vertex(2)=oldVertex(2)+add;
            % update PWL model
            PWLModel.Points(vertexIndex,2)=vertex(2);
   
        %     %FOR DEBUG
        %     surf(x,y,z)
        
            % if vertex movement exceed limits
            if (or(vertex(2)<Ymin,vertex(2)>Ymax));
                % undo moving vertex
                PWLModel.Points(vertexIndex,2)=oldVertex(2);
                vertex(2)=oldVertex(2);
                PWLModel.RevertCount=PWLModel.RevertCount+1;
                disp('WARNING: revert vertex move on y axis!')
            else
                break
            end
        end
    else
        vertex(2)=oldVertex(2);
    end
    
    %-----------------------------------------------Z
    % create bounding box
    Xvalues=[vertex(1)-0.5, vertex(1), vertex(1)+0.5];
	Yvalues=[vertex(2)-0.5, vertex(2), vertex(2)+0.5];
    [MeshX,MeshY] = meshgrid(Xvalues,Yvalues);
    % find Z-axis value with interpolation
    MeshZ=Interpolated_Ref(MeshX,MeshY);
    % find out the minimum and maximum value
    Zmin=min(min(MeshZ));
    Zmax=max(max(MeshZ));
    Ztoleranz=Zmax-Zmin;

    Zufall=(2*rand(1)-1);
    add=Ztoleranz*Zufall;
    vertex(3)=MeshZ(2,2)+add;
    

else 
    % full random version
    %-----------------------------------------------X
    if (or(PWLModel.Index(vertexIndex)==0,PWLModel.Index(vertexIndex)==2))
        while 1

            Zufall=(2*rand(1)-1);
            add=relativeSchrittweite.X*Zufall;

            % move selected vertex
            vertex(1)=oldVertex(1)+add;
            % update PWL model
            PWLModel.Points(vertexIndex,1)=vertex(1);

            % if vertex movement exceed limits
            if (or(vertex(1)<Xmin,vertex(1)>Xmax));
                % undo moving vertex
                PWLModel.Points(vertexIndex,1)=oldVertex(1);
                vertex(1)=oldVertex(1);
                PWLModel.RevertCount=PWLModel.RevertCount+1;
                disp('WARNING: revert vertex move on x axis!')
            else
                break
            end
        end
    else
        vertex(1)=oldVertex(1);
    end
    
    %-----------------------------------------------Y
    if (or(PWLModel.Index(vertexIndex)==0,PWLModel.Index(vertexIndex)==1))
        while 1
            
            Zufall=(2*rand(1)-1);
            add=relativeSchrittweite.Y*Zufall;

            % move selected vertex
            vertex(2)=oldVertex(2)+add;
            % update PWL model
            PWLModel.Points(vertexIndex,2)=vertex(2);

            % if vertex movement exceed limits
            if (or(vertex(2)<Ymin,vertex(2)>Ymax));
                % undo moving vertex
                PWLModel.Points(vertexIndex,2)=oldVertex(2);
                vertex(2)=oldVertex(2);
                PWLModel.RevertCount=PWLModel.RevertCount+1;
                disp('WARNING: revert vertex move on y axis!')
            else
                break
            end
        end
    else
        vertex(2)=oldVertex(2);
    end
    %-----------------------------------------------Z
    Zufall=(2*rand(1)-1);
    add=relativeSchrittweite.Z*Zufall;
    vertex(3)=oldVertex(3)+add;
end

    newRevertCount=PWLModel.RevertCount;

end
    

