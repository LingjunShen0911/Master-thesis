function [Points, Lines, Triangles,variantName] = read_xml(inputFilename)
% reads a xml-model
%
% Syntax:  
%    [Points, Lines, Triangles] = read_xml(inputFilename)
%
% Inputs:
%   inputFilename: filename of a xml-file
%
% Outputs:
%   Points:     Array of all Points, Point(n,x) = xth (1-3) coordinate of point n 
%   Lines:      Array of all Lines,  Lines(n,x) = Index of xth (1-2) point of line n 
%   Triangles:  Array of all Triangles, Triangles(n,x) = Index of xth (1-3)
%               Line of Triangle n
%   
% Example
%   read_xml('../examples/model.xml')
%
% Author: Enno Röhrig
% Written: November-2014
% Last update: 6-January-2015
% Last revision: 1.1

try
   xml = xmlread(inputFilename);
catch
   error('Failed to read XML file %s.',inputFilename);
end

%try
   
   
   Vertices = xml.getElementsByTagName('vertex');
   Points = zeros(Vertices.getLength,3);
   for vertex_nr = 1:Vertices.getLength
       
        vertex=Vertices.item(vertex_nr-1);
        Attributes = vertex.getAttributes;
        Points(vertex_nr,1) = str2double(Attributes.item(1).getValue);
        Points(vertex_nr,2) = str2double(Attributes.item(2).getValue);
        Points(vertex_nr,3) = str2double(Attributes.item(3).getValue);

   end
   
   
   TrianglesNodes = xml.getElementsByTagName('triangle');
   
   Triangles = zeros(TrianglesNodes.getLength,3);
   Lines = zeros(length(Points(:,1))*5,2);
   LinesFrom  = zeros(length(Points(:,3)),5);          %2d structure: Index=startpoint, values=endpoints from all lines
   LinesIndex  = zeros(length(Points(:,3)),5);     %2d structure: LinesHelperIndex(p1,i) = Index from the Linie (p1 to  LinesHelper(p1,i)) in the Matrix "Lines"
   maxIndex = 0;
   P = zeros(3);
   
   for triangle_nr = 1:TrianglesNodes.getLength
        triangle=TrianglesNodes.item(triangle_nr-1);
        Attributes = triangle.getAttributes;
        P(1) = str2num(Attributes.item(1).getValue)+1;
        P(2) = str2num(Attributes.item(2).getValue)+1;
        P(3) = str2num(Attributes.item(3).getValue)+1;
%         P=sort(P);
        
        start_nr = [1,1,2];                         %all kombination of points, which represent a line
        end_nr   = [2,3,3];                          %1-2, 1-3, 2-3
        for l=1:3   
            P1=P(start_nr(l));
            P2=P(end_nr(l));
            LineExist=0;
            
            To=1;
            while To <= length(LinesFrom(P1,:)) && LinesFrom(P1,To) ~= 0
                
                if LinesFrom(P1,To) == P2
                    LineExist = 1;
                    break;
                end
                To=To+1;
            end
            
            if LineExist == 1
%                 Triangles(triangle_nr, l) = LinesIndex(P(start_nr(l)),To);
            else
                maxIndex=maxIndex+1;
%                 Triangles(triangle_nr, l) = maxIndex;
                
                LinesFrom(P1,To)  = P2;
                LinesIndex(P1,To) = maxIndex;
               
                Lines(maxIndex,1) = P1;
                Lines(maxIndex,2) = P2;
            end
        Triangles(triangle_nr, l)=P(l,1);
       end
   end
   Lines = Lines(1:maxIndex,:);
   
%catch
 %  error('Unable to parse XML file %s.',inputFilename);
%end

variantName=xml.getElementsByTagName('element').item(0).getAttributes.item(1).getValue;
return;            
end





