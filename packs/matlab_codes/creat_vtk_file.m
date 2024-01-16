function creat_vtk_file(coord,elem,P,S,step)

ext = '.vtk';
fname = 'out';
step = num2str(step);

nnode = size(coord,1);
nelem = size(elem,1);

fname_vtk = [fname '00' step ext];
fid = fopen(fname_vtk,'w'); % Output 'w'riting file
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n',nnode);
data1 = [coord(:,1:3)'];
fprintf(fid,'%26.16E %26.16E %26.16E \n',data1);
fprintf(fid,'CELLS %i %i\n',nelem,4*nelem);
numdata = 3*ones(nelem,1);
data2 = [numdata';(elem(:,1)-1)';(elem(:,2)-1)';(elem(:,3)-1)'];
fprintf(fid,'%i %i %i %i\n',data2);
fprintf(fid,'CELL_TYPES %i\n',nelem);
type_elem = 5*ones(nelem,1);            % 5 = triangles
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