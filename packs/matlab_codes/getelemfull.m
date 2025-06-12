%--------------------------------------------------------------------------
%                           FUNCTION "getelemfull"
%--------------------------------------------------------------------------
%This is a modified version of the  original getelem function. 
%It reads the elements in a MSH file.
%This modification allows reading partioned meshes.
%--------------------------------------------------------------------------
%                                   INPUT
%--------------------------------------------------------------------------
% veryfymshfile - pointer to MSH file
% nnode - number of nodes
% numwell - number of wells
% well - information on well - matrix with the flags as stated in start.dat
%--------------------------------------------------------------------------
%                                   OUPUT
%--------------------------------------------------------------------------
%
% elem - matrix containing all elements (triangles and quadrangles only)
%        NODE1  NODE2  NODE 3 NODE 4  FLAG
% Node 1 to 4 takes the node that describe the material. if a triangle
% is described, NODE4 takes 0. FLAG is the material flag
% nbe - number of edge-type entities located in the boundary
% nelem - number of elements
% nodelim - node-type entity on the boundary
% intnode - node-type entities located inside the mesh
% flaglim - vector with internal node flags
% inboundedge - edges-type entity inside the domain
% meshtype- true: multiscale mesh false: fem mesh
% bcflag - matrix with boundary conditions flag from start.dat
% elemloc - vector with the partion of each element
% npar - number of partions
% pboundel coarseelem - list all elements inside a coarse cell
% ex: coarseelem (3) = list of elements on the boundary of partion 3
% ghostelem - list elements neighbor n neighbor to m
% ex: list: ghostelem(3,4) = list of elements of 4 that are neighbors of 3

function [elem,nbe,nelem,nodelim,intnode,flaglim,inboundedge, meshtype, elemloc,npar,coarseelem ,ghostelem ] = ...
    getelemfull(verifymshfile,nnode,numwell,well)
%Identify on the *.msh file a like-element type. Actualy, the *.msh
%gives entities which hold points, edges and elements (like-elements).
%This entities will all of them be reads.
%Open the *.msh file.
readmsh = fopen(verifymshfile);
        
%"nent" is the number of all entities understud as element in *.msh file
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
%Attribute the data to "nent" (number of entity, that is: points, lines 
%and elements)
nent = getmshdata{1};
%Once we have "nent", both like-element and external edges entities may be 
%located. "meshtype" is the type of element used in domain discretization.
%2 ==> triangles; 3 ==> quadrangles. But that will be used also to define
%how many points entities there on the *.msh file

%captures 13 integers  (max number of flags in a row)
getmshdata = textscan(readmsh,'%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%*[^\n]',nent);
fclose(readmsh);
%Attribute the data to "auxmat"
auxmat = cell2mat(getmshdata);
clear getmshdata 
%deletes all columns filled with NaN only 
auxmat = auxmat(:,not(all(isnan(auxmat),1)));

%% gets the number of each type entities
% ntri = number of triangles
% nquad = number of quadrangles
% nline = number of lines
% npoint = number of points

ntri = sum(auxmat(:,2) == 2);
nquad = sum(auxmat(:,2) == 3);
nline = sum(auxmat(:,2) == 1);
npoint = sum(auxmat(:,2) == 15);
%finding the number of 2d elements
nelem = ntri + nquad;

%finding the max number of flags among all elements
nflagmax = max(auxmat(:,3));
%checking the how many partions there are

if nflagmax > 2
   
    meshtype = true;
    npar = max(auxmat(:,7));
    
    %allocate space for ghostelem
    %pboudnel - list all elements in a coarse cell
    % ex: coarseelem (3) = list of elements on the boundary of partion 3
    coarseelem  = cell(npar,1);

    %ghostelem - list elements neighbor n neighbor to m
    %ex: list: ghostelem(3,4) = �list of elements of 4 that are neighbors of 3
    ghostelem  = cell(npar);

    %allocating space for elemloc
    elemloc = zeros(nelem,1);
    
else
    meshtype = false;
    npar = 0;
    coarseelem  = cell(0);
    ghostelem  = cell(0);
    elemloc = zeros(0);

    
    
end


if nquad >0
    interval = [-3,-2];
else
    interval = [-2 ,-1];
end


%capture the edges inside the domain  
inboundedge = auxmat(logical(auxmat(:,2) ~= 15 & auxmat(:,4) > 1000),end + interval(1):end + interval(2));
%sorts inboundedge
inboundedge = sort(inboundedge,2);
%define the number of edge-type entities on the bondary.
nbe = nline - size(inboundedge,1);
%Get the amount of internal node (scalar)
intnode = sum(logical(auxmat(:,2) == 15 & auxmat(:,4) > 1000));
%We exclude any node inside the domain.    
nodelim = npoint - intnode;
%vector with all flags associated with external nodes
flaglim = auxmat(logical(auxmat(:,2) == 15 & auxmat(:,4) < 1000),4);

%Allocate space for elem
elem = zeros(nelem,5);

celem = 1;


%to be rewritten to improve efficency
for index = 1: size(auxmat(:,1))
    auxvec = auxmat(index,:);
    pe = 0; el1 = 0; el2  = 0; el3 = 0; el4 = 0;
    %check each line line for the element-type described
    % change here to pick information on internal boundaray
    %debug
    
    if auxvec(2) == 2 | auxvec(2) == 3
            pe = auxvec(4);
            %num of flags 
            fcount = auxvec(3);
            %setting min interval
            intmin = 3 + fcount +1 ;
            %malha puramente triangular
            if ntri > 0 & nquad == 0
                intmax = intmin + 2;
                fcell = num2cell([auxvec(intmin:intmax)]);
                [el1,el2,el3] = fcell{:};
                
                %if mesh is a multiscale mesh
                if meshtype == true
                    elemloc(celem) = auxvec(intmin-auxvec(6));
                    for i = 7: (7 +auxvec(6)-1)
                        if i == 7
                            coarseelem {auxvec(7)}(end+1) = celem;
                        else
                            ghostelem{auxvec(7),abs(auxvec(i))}(end+1) = celem;
                        end
                    end

                end
   
            %malha triangular + quad 
            else
                intmax = intmin + 3;
                %[intmin intmax]
              
                %auxvec
                t=auxvec(intmin:intmax);
                el1 = t(1);
                el2 = t(2);
                el3 = t(3);
                el4 = t(4);
                %fcell = num2cell([auxvec(intmin:intmax)]);
                %[el1,el2,el3,el4] =  fcell{:};
               
                %if mesh is a multiscale mesh
                if meshtype == true
                    elemloc(celem) = auxvec(7);
                    for i = 7: (7 +auxvec(6)-1)
                        if i == 7
                            coarseelem {auxvec(7)}(end+1) = celem;
                        else
                            ghostelem{auxvec(7),abs(auxvec(i))}(end+1) = celem;
                        end
                    end
                end
                

            end
        
           
            elem(celem,:) = [el1,el2,el3,el4,pe];
            celem = celem +1;
    end
    
end


end
