function creat_vtk_file_quad(coord,elem,P,S,step)


ext = '.vtk';
fname = 'out';
step = num2str(step);

nnode = size(coord,1);
nelem = size(elem,1);
sumnodebyelem = sum(sum(elem(:,1:4) ~= 0));
countriang = sum(elem(:,4) == 0);
countquad = sum(elem(:,4) ~= 0);
type_elem = zeros(size(elem,1),1);
fname_vtk = [fname '00' step ext];
fid = fopen(fname_vtk,'w'); % Output 'w'riting file
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n',nnode);
data1 = [coord(:,1:3)'];
fprintf(fid,'%26.16E %26.16E %26.16E \n',data1);
fprintf(fid,'\r\n');
fprintf(fid,'CELLS %i %i \r\n\r\n',nelem,sumnodebyelem + nelem);
numdata = 4*ones(countquad,1);
data2 = ...
        [numdata';(elem(countriang + 1:nelem,1) - 1)';...
        (elem(countriang + 1:nelem,2) - 1)';...
        (elem(countriang + 1:nelem,3) - 1)';...
        (elem(countriang + 1:nelem,4) - 1)'];
    fprintf(fid,'%i %i %i %i %i \r\n',data2);
    %Definition of element type (7 is a code to quadrangle)
    type_elem(countriang + 1:nelem) = 7;  
fprintf(fid,'\r\n');

fprintf(fid,'CELL_TYPES %i\n',nelem);
fprintf(fid,'%i\n',type_elem);
fprintf(fid,'CELL_DATA %i\n',nelem);
%Pressure
fprintf(fid,'SCALARS Pressure float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'\t %26.16E\n',P);
%Saturation
fprintf(fid,'SCALARS Saturation float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'\t %26.16E\n',S);
fclose(fid);
end