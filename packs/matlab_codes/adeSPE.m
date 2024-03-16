%adequacao para rodar SPE
% wells = [        6490           1         301           1           0           1;           
% 	       1           2           0           0         501           1;
%        12981           3           0           0         502           1;
%        13200           4           0           0         503           1;
%          220           5           0           0         504           1];
%      

layer = 60;
 sizeX = 1;
 sizeY = 1;
 %sizeX = 2200;
 %sizeY = 1200;
meshX = 220;
meshY = 60;
Rotacionado = 0;  
triangular = 0;
%%stretching
tflag = max(coord);
if (tflag(1) < sizeX) & (tflag(2) < sizeY)
    coord(:,1) = coord(:,1) * sizeX ;
    coord(:,2) = coord(:,2) * sizeY ;
    centelem(:,1) = centelem(:,1) * sizeX ;
    centelem(:,2) = centelem(:,2) * sizeY ;
end
% end
%bedge(:,end-1:end)= 201;
%superFolderMKSPE
SPEField;


  