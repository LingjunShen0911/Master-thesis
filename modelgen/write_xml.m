function [] = write_xml(Points, Triangles, outputFilename, elementName, variantName) 
% writes a xml-model
%
% Syntax:  
%   write_xml(Points, Lines, Triangles, outputFilename)
%
% Inputs:
%   Points:     Array of all Points, Point(n,x) = xth (1-3) coordinate of point n 
%   Lines:      Array of all Lines,  Lines(n,x) = Index of xth (1-2) point of line n 
%   Triangles:  Array of all Triangles, Triangles(n,x) = Index of xth (1-3)
%               Line of Triangle n
%   outputFilename: filename of a gts-file
%
% Outputs:
%   xml-file
%   
% Example
%   write_xml([1,1,1;2,2,2;3,3,3],[1,2;2,3;1,3],[1,2,3],'../examples/model.xml')
%
% Author: Enno Röhrig
% Written: November-2014
% Last update: 6-January-2015
% Last revision: 1.1

fxml = fopen(outputFilename, 'w');
	fprintf(fxml, '<praise>\n');
	fprintf(fxml, '<pwlmodels>\n');
	fprintf(fxml, ['<element name="%s" variant="' variantName '">\n'], elementName);
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

    for (i = 1:size(Points(:,1)))
		x1 = Points(i,1);
		x2 = Points(i,2);
		y  = Points(i,3);
        fprintf(fxml, '			<vertex id="%d" x1="%.16g" x2="%.16g" y="%.16g"/>\n', (i - 1), x1, x2, y);
	end;
	fprintf(fxml, '		</vertices>\n');

	fprintf(fxml, '		<triangles>\n');
	for (i = 1:1:size(Triangles(:,1)))
		vertex0 = Triangles(i,1) -1;
		vertex1 = Triangles(i,2) -1;
		vertex2 = Triangles(i,3) -1;
%         if(vertex2 == vertex1 || vertex2==vertex0) 
%             vertex2 = Lines(Triangles(i,2),2)-1;
%         end
		fprintf(fxml, '			<triangle id="%d" vertex0="%d" vertex1="%d" vertex2="%d"/>\n', (i - 1), vertex0, vertex1, vertex2);
	end;
	fprintf(fxml, '		</triangles>\n');

    fprintf(fxml, '	</model2dim>\n');
	fprintf(fxml, '</element>\n');
	fprintf(fxml, '</pwlmodels>\n');
	fprintf(fxml, '</praise>\n');
fclose(fxml);

end

