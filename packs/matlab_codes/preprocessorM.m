%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 30/11/2011
%Modify data: 11/01/2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza; Luiz Eduardo; Tulio Cavalcante; Artur Castiel
%--------------------------------------------------------------------------
%Goals:
%The main goal of this numerical routine is the creation of a mesh in a...
%determinate domain by gmsh and after that, load the matrix of
%coordinate, elements and both iner and boundary edge. All of them will be
%used to produce the numerical statment of each method studed by time

%--------------------------------------------------------------------------
%This FUNCTION execute the gmsh with input data about domain's geometry and 
%options parameter. After that, generate an output file with informations 
%about mesh, which will be read by this numerical routine 
%--------------------------------------------------------------------------

function [coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
    normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,kmap,wells,elemloc,npar,coarseelem ,ghostelem,...
    multiCC, pmethod,mshfile, keymsfv,interptype,keypath1,smethod,xyrz,r0,symaxe,coarseratio,auxcvfactor,...
    nonlinparam,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,mainpath,resfolder,benchkey,...
    klb,limiterflag,rowposit] = preprocessorM 
global osMode
%--------------------------------------------------------------------------
%It reads the data file (from "Start.dat")

[head,benchkey,pathstore,numcase,complemfile,phasekey,xyrz,r0,symaxe,...
    keymsfv,coarseratio,pmethod,auxcvfactor,interptype,nonlinparam,smethod,...
    multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,lsexp,...
    limiterflag,keyfile,keypath1,mshfile,keypath2,basename,keymesh,numparx,...
    transflabelx,numpary,transflabely,stdclflag,maxelemsize,minelemsize,...
    keyelem,keyalgor,keyrefine,nref,scalemesh,keyopen,keydens,densname,...
    dens,keyvisc,viscname,visc,satlimit,keyperm,klimits,hipermname,...
    numlayer,kmap,keypor,pormap,porname,keygravity,g,keycapil,ncaplcorey,...
    bcflag,numbcnode,bcnodepointer,courant,totaltime,numwell,well,multiCC] = ...
    getdatafile;

%--------------------------------------------------------------------------
%Generating or acessing the folder  

%Verify if the path (with diretory) exists.
%It does not exist. Must be created.
if exist(char(pathstore),'dir') == 0
    %Create a folder to put Results
    if strcmp(osMode, 'linux')
        command = ['mkdir ' char(pathstore)];
        %It calls system
    elseif strcmp(osMode, 'windows')
        command = ['md ' char(pathstore)];        
    end    
    system(command);
end  %End of IF

%Verify if ther is a complement file. If it there exists, a folder will be
%created inside the "storepath".
%There is NOT ("0" was written)
if strcmp(complemfile,'0')
    %Attribute to "innerfolderpath" the "pathstore". There is no subfolder   
    mainpath = char(pathstore);
    %Define "resfolder"
    resfolder = 'Results';
 
    if strcmp(osMode, 'linux')
        %Define an auxiliary command (for RESULTS):
        commandres = [mainpath '/' 'Results'];
        %Define an auxiliary command (for IMAGES):
        commandimg = [mainpath '/' 'Images'];
        %Define an auxiliary command (for TABLES):
        commandtbl = [mainpath '/' 'Tables'];
    elseif    strcmp(osMode, 'windows')
        commandres = [mainpath '\' 'Results'];
        %Define an auxiliary command (for IMAGES):
        commandimg = [mainpath '\' 'Images'];
        %Define an auxiliary command (for TABLES):
        commandtbl = [mainpath '\' 'Tables'];
    end
%There is subfolder. 
else
    %Create a subfolder to put Results:
    %Define an auxiliary command
    auxcommand = [char(pathstore) '/' char(complemfile)];
    %Verify if it already exists
    if exist(auxcommand,'dir') == 0
        %Define command
        if strcmp(osMode, 'linux')
        command = ['mkdir ' auxcommand];
        elseif strcmp(osMode, 'windows')
        %It calls system
        command = ['md ' auxcommand];
        end
        system(command);
    end  %End of IF

     
    if strcmp(osMode, 'linux')
        %Attribute to "innerfolderpath" the "pathstore" + a subfolder
        mainpath = [char(pathstore) '/' char(complemfile)];
        %Define "resfolder"
        resfolder = ['Results_' char(complemfile)];
        %Define an auxiliary command (for RESULTS):
        commandres = [mainpath '/' resfolder];
        %Define an auxiliary command (for IMAGES):
        commandimg = [mainpath '/' 'Images_' char(complemfile)];
        %Define an auxiliary command (for TABLES):
        commandtbl = [mainpath '/' 'Tables_' char(complemfile)];
        %It calls system
    elseif strcmp(osMode, 'windows')
        %Attribute to "innerfolderpath" the "pathstore" + a subfolder
        mainpath = [char(pathstore) '\' char(complemfile)];
        %Define "resfolder"
        resfolder = ['Results_' char(complemfile)];
        %Define an auxiliary command (for RESULTS):
        commandres = [mainpath '\' resfolder];
        %Define an auxiliary command (for IMAGES):
        commandimg = [mainpath '\' 'Images_' char(complemfile)];
        %Define an auxiliary command (for TABLES):
        commandtbl = [mainpath '\' 'Tables_' char(complemfile)];
    end    
    
end  %End of IF

%Create folders to store results, images and tables.
%Once it does not exist, we will create it.
if exist(commandres,'dir') == 0
    %Define command
    if strcmp(osMode, 'linux')
        command = ['mkdir ' commandres];
        %It calls system
    elseif strcmp(osMode, 'windows')
        command = ['md ' commandres];        
    end    
    %It calls system
    system(command);
end  %End of IF

%Once it does not exist, we will create it.
if exist(commandimg,'dir') == 0
    %Define command
    if strcmp(osMode, 'linux')
        command = ['mkdir ' commandimg];
    elseif strcmp(osMode, 'windows')
        command = ['md ' commandimg];
    end
        %It calls system
    system(command);
end  %End of IF

%Once it does not exist, we will create it.
if exist(commandtbl,'dir') == 0
    %Define command
    if strcmp(osMode, 'linux')
        command = ['mkdir ' commandtbl];
    elseif strcmp(osMode, 'windows')   
        command = ['md ' commandtbl];
    end
        %It calls system
    system(command);
end  %End of IF

%--------------------------------------------------------------------------
%When we have a *msh file: 
    
%You have the mesh file. Once the mesh already exists, gmsh doesn't 
%need to be used.
if strcmp(keyfile,'y')  
    %In this option is necessary to write a path.
    %Attribute to "verifymshfile" the path and the name of file.
    if strcmp(osMode, 'linux')
       verifymshfile = [char(keypath1) '\' char(mshfile)];
    elseif strcmp(osMode, 'windows')
       verifymshfile = [char(keypath1) '\' char(mshfile)];
    end
    %Warning to user!
    if exist(verifymshfile,'file') == 0
        disp(' ');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('>>> WARNING!');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('The *.msh file was not found!!!');
        disp('Verify either the path or the name attributed to file.');
        disp('---------------------------------------------------------------');
    end  %End of internal IF
        
%We do not have a *msh file:
%In this case, gmsh is used to produce the mesh. However, it is not
%executed (User Interface) for such finality.
elseif strcmp(keyfile,'n')  
    %Attribute to "verifygeofile" the path AND the file's name (with *.geo)
    if strcmp(osMode, 'linux')
       verifygeofile = [char(keypath2) '\' char(basename)];
    elseif strcmp(osMode, 'windows')
       verifygeofile = [char(keypath2) '/' char(basename)];
    end
    %Warning to user!
    if exist(verifygeofile,'file') == 0
        disp(' ');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('>>> WARNING!');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('The *.geo file not found!!!');
        disp('Check the path or the name attributed to file.');
        disp('---------------------------------------------------------------');
    end  %End of internal IF
    
    %Obtain the name of GEOMETRY file which will be stored in the new 
    %folder created.
    %Get only the base name (without the extension *.geo)
    getbasename = char(basename);
    getbasename = getbasename(1:length(getbasename) - 4);
    
    %Obtain the name of GEO file to be stored in the new folder.
    geofile = [getbasename '.geo']';
    %Obtain the name of MESH file also stored in the new folder. 
    mshfile = [getbasename '.msh']';
    %Obtain the name of the OPTION file.
    optfile = [getbasename '.opt']';
        
    %----------------------------------------------------------------------
    %Create and guard the new GEOMETRY and OPTIONS file.
    if strcmp(osMode, 'linux')
        fS = '/';
    elseif strcmp(osMode, 'linux')
        fS = '\';
    end

    %Open the file "*.geo" with its respective path
    geometry = fopen([mainpath fS geofile'],'w');

    %Open the file "*.opt" with its respective path
    options = fopen([mainpath fS optfile'],'w');
    
    %Get the mesh generated. It sends the mesh to data structure function.
    %See "getcoord", "getelem" and others.
    verifymshfile = [mainpath fS mshfile'];
    
    %Attribute this informations to "*.geo" created
    for i = 1:length(head);
        %Write in the geometry file
        fprintf(geometry,'%s\r\n',char(head(i))');
        %Write in the options file
        fprintf(options,'%s\r\n',char(head(i))');
    end  %End of FOR
       
    %Read the BASE geometry file and copy the date to new geometry file 
    %inside project folder created.

    %Get the date and hour to knoledge of user
    [date] = clock;
       
    %Write the date and hour in the geometry file
    fprintf(geometry,'//Create date: %u/%u/%u;\thour: %d:%dh',...
        date(1),date(2),date(3),date(4),date(5));
    %Feature the geometry file
    fprintf(geometry,'\r\n\r\n%s\r\n',char(head(1))');
    fprintf(geometry,'//This file has CAD parameters. It is related to building of domain\r\n\r\n');
       
    %Write the date and hour in the options file
    fprintf(options,'//Create date: %u/%u/%u;\thour: %d:%dh',...
        date(1),date(2),date(3),date(4),date(5));
    %Feature the options file 
    fprintf(options,'\r\n\r\n%s\r\n',char(head(1))');
    fprintf(options,'//This file has OPTIONS parameters. Configure GMSH environment\r\n\r\n');

    %OBS.: These files will be closed when each one is complete. 
    %The command "fclose" will be used.        
        
    %----------------------------------------------------------------------
    %Parameters of MESH:

    %If "keymesh" receives "n", configure unstructured mesh
    if strcmp(keymesh,'n')
        %Build the "*.opt" file with unstructured data:
            
        %MAX ELEMENT SIZE:
        %Write a head of function. Using "//" the gmsh understand as a
        %comment.
        fprintf(options,'//Maximum mesh element size:\r\n//1e+22 ==> [Default]\r\n\r\n');
        %Write the function necessary to gmsh configuration. 
        fprintf(options,'Mesh.CharacteristicLengthMax = %f;\r\n\r\n',...
            maxelemsize);
            
        %MIN ELEMENT SIZE:
        %This parameter is write again in the new "*.opt"
        %Write head of parameter
        fprintf(options,'//Minimum mesh element size:\r\n//0 ==> [Default]\r\n\r\n');
        %Write the function
        fprintf(options,'Mesh.CharacteristicLengthMin = %f;\r\n\r\n',...
            minelemsize);
            
        %OBS.: Other parameters are out of this loop because they are 
        %of general interest (used in both structured and unstructured 
        %mesh). 
    end  %End of IF

    %Type of mesh:
    
    %The element is either a triangle or a quadrangle
    if keyelem == 0
        %Write the head of the function
        fprintf(options,'//Type of mesh used:\r\n//0 ==> [AllTriangles]; 1 ==> [AllQuadrangles]\r\n\r\n');
        %Write the function
        fprintf(options,'Mesh.SubdivisionAlgorithm = %u;\r\n\r\n',...
            keyelem);
    
    %The element is quadrangular (unstructured mesh)
    elseif keyelem == 1 && strcmp(keymesh,'n')
        %Write the head of the function related to triangle
        %building
        fprintf(options,'//Type of mesh used:\r\n//0 ==> [AllTriangles]; 1 ==> [AllQuadrangles]\r\n\r\n');
        %Write the function related to triangle building
        fprintf(options,'Mesh.SubdivisionAlgorithm = 1;\r\n\r\n');
        
    %The option choosen was "BOTH". Will be created both mesh (triang. and 
    %quadrangular). The algorithm will write the function as trinangle and 
    %in sequence to write a new routine (Mesh.RecombineAll)
    elseif keyelem == 2
        %Write the head of the function related to triangle
        %building
        fprintf(options,'//Type of mesh used:\r\n//0 ==> [AllTriangles]; 1 ==> [AllQuadrangles]\r\n\r\n');
        %Write the function related to triangle building
        fprintf(options,'Mesh.SubdivisionAlgorithm = 0;\r\n\r\n');
        %Write the new function "Mesh.RecombineAll"
        fprintf(options,'//Apply recombination algorithm to all surfaces:\r\n//0 ==> [Default]\r\n\r\n');
        %Write the new function
        fprintf(options,'Mesh.RecombineAll = 1;\r\n\r\n');
    %Disp a mesage to user.
    elseif keyelem > 2
        disp('Incorrect number related to mesh type! Try again');
    end  %End of IF (type of element)

    %OBS.: The variables above will be used to rewrite the new "*.geo" in 
    %the end of section "rename physical properties". Case the option 
    %structured mesh be choose, some command line will be written in the 
    %new "*.geo" (see line 460)
        
    %Type of algorithm used:
    %[1] ==> "MeshAdapt"; [5] ==> "Delaunay"; [6] ==> "Frontal"
    
    %Write the head of this function
    fprintf(options,'//Type of algorithm used:\r\n//1 ==> [MeshAdapt]; 5 ==> [Delaunay]; 6 ==> [Frontal]\r\n\r\n');
    %Function statment:
    fprintf(options,'Mesh.Algorithm = %u;\r\n\r\n',keyalgor);
        
    %Type of mesh refinement
    
    %You can choose [0] ==> "No split"; [1] ==> "Automatic
    %head of function
    fprintf(options,'//Type of remeshing algorithm:\r\n//0 ==> [No split]; 1 ==> [Automatic]\r\n\r\n');
    %Statment of function
    fprintf(options,'Mesh.RemeshAlgorithm = %u;\r\n\r\n',keyrefine);
        
    %Apply a SCALE FACTOR in MESH SIZE
    
    %This last option apply a scale in mesh size (const. in all mesh).
    %head of function
    fprintf(options,'//Factor applied to all element sizes:\r\n//1 ==> [Default]\r\n\r\n');
    %Write function
    fprintf(options,'Mesh.CharacteristicLengthFactor = %f;\r\n\r\n',...
        scalemesh);

    %----------------------------------------------------------------------
    %Rename physical propertie (both line and surface), in the new "*.geo".

    %Open the ORIGINAL *.geo file:
    readgeo = fopen(verifygeofile);

    %"getgeodata" reads each row of *.geo and store it as a string array. 
    getdata = textscan(readgeo,'%s','delimiter',';');
    %Attribute the data to "getgeodata"
    getgeodata = getdata{1};
    
    %Close the *.geo file
    fclose(readgeo);
    
    %"lengthgeodata" get the length of "getgeodata" vector. This 
    %information is used in the "for" below. 
    lengthgeodata = length(getgeodata);

    %Initialize the parameter "iplane", "iedge" and "ipoint". This is 
    %the order which the permeability and boundary conditions are 
    %distributed respectivaly.
    iplane = 1;
    iedge = 1;
    ipoint = 1;
    iline = 1;
    igeopoints = 1;  %Is a counter of coordinates which are accumulated
        
    %"igetloop" is a counter of howmany loops there in the *.geo file.
    %This is important to associate each surface to nodes which
    %constitute them. Below, the parameter is initialized.
    igetloop = 1;

    %Compare the nodes whose each boundary condition in attributed with 
    %node which constitute each line in *.geo. Depending of comparation, 
    %it changes the label of edge (or line) in *.geo
    %Swept the *.geo file
        
    %----------------------------------------------------------------------
    %Catch the numbers in each ROW evaluated
            
    for i = 1:lengthgeodata
        %Initialize "linelabel" to each line of *.geo evaluated
        linelabel = 0;
        %"rowgeo" take the whole command line evaluated in the *.geo
        rowgeo = char(getgeodata(i));
        %"takeletter" takes all characteres of the current row (*.geo) 
        takeletter = rowgeo;
        %"isanumber" recognaze only the characters which are numbers.
        %That has the same size of "takeletter" (FANTASTIC COMMAND)
        isanumber = isstrprop(takeletter,'digit');
        %"isnumberpos" takes the position of vector "isanumber" whose
        %number identify a numeric character.
        isnumberpos = find(isanumber(:) == 1);
        %"isdotpos" finds the position of a point in each line
        %of *.geo evaluated
        isdotpos = find(rowgeo == '.');
        %"ilabel" is the amount of numbers in each command line
        ilabel = 1;
        %"iletter" is a counter of letters in each row of *.geo 
        %evaluated 
        iletter = 1;
        %Initialize the sign of number.
        sign = 1;

        %Swept all lines of *.geo
        while iletter <= length(rowgeo)
            %"getletter" catch each character of line evaluated.
            %The next step is verify if this letter is or not a
            %number. For this, the functions "str2num" and "isempty"
            %will be used.
            getletter = str2num(rowgeo(iletter));
            %The letter evaluated does not is a number or is "i",
            %which would mean a complex number. In these cases 
            %"iletter" is increased
            if isempty(getletter) || ...
                    strcmp(rowgeo(iletter),'i')
                %Increment "iletter"                       
                iletter = iletter + 1;
            %Change the value of the sign when exists a
            %negative number.
            if strcmp(rowgeo(iletter),'-')
                sign = -1;
            end  %End of IF
            %If the function "ismepty" does not equal to one, "getletter"
            %received a not null value. The letter evaluated is a number 
            elseif isempty(getletter) == 0 
                %"ahead" receives the same value that "iletter"
                ahead = iletter;
                %Swept the other letters to evaluate if the next
                %letters is also a number
                while ahead <= length(rowgeo) 
                    %The field of eval. is between "iletter" and "ahead"
                    getletter = str2num(rowgeo(iletter:ahead));
                    %In this verification, maybe a ahead letter may
                    %be not a number. In this case, the procedure
                    %is aborted with "break"
                    if isempty(getletter) || ...
                            strcmp(rowgeo(ahead),',')
                        break;
                    end  %End of IF
                    %"getnumber" receives the value catched by 
                    %"getletter"
                    getnumb = getletter;
                    %"ahead" is encreased to evaluate other letter
                    ahead = ahead + 1;
                end  %End of WHILE
                %Fill "linelabel". This parameter store all numbers
                %inside each line evaluated.
                linelabel(ilabel) = getnumb*sign;
                %Jump to next number
                ilabel = ilabel + 1;
                iletter = ahead + 1;
                %"sign" is recovered again
                sign = 1;
            end  %end of IF
        end  %End of external WHILE
        
        %------------------------------------------------------------------
        %Begins to evaluate each command line of original *.geo
        
        %Create more parameter "cl" as a function of wells seted in
        %"Start.dat". Just have sense modify this parameter if the
        %mesh is unstructured, for this the second condition in the
        %decision command.
        if strcmp(rowgeo(1:2),'cl') && strcmp(keymesh,'n') 
            %Print the information about "cl1"
            fprintf(geometry,'//"cl1" corresponds to element size attributed in "Start.dat";\r\n');
            %Print the original element size statment ("cl1"). "stdflag"
            %is a parameter which store the value attributed to
            %standard "cl" in the "Start.dat"
            fprintf(geometry,'cl1 = %f;\r\n',stdclflag);
            %If there is no well (numwell = 0), other variable "cl"
            %does not will set.
            if numwell > 0
                for iwell = 1:numwell
                    fprintf(geometry,'//"cl%u" corresponds to element size surrounding the well center %u;\r\n',...
                        iwell + 1,iwell);
                    fprintf(geometry,'cl%u = %f;\r\n',...
                        iwell + 1,well(iwell,4));
                end  %End of FOR
            end  %End of IF
            %Jump a row in the new *.geo
            fprintf(geometry,'\r\n');

        %Modify the command lines associated with each POINT. Change 
        %the parameter "cl1" which means the element size associate 
        %to point evaluated.
        %OBS.: The modification will be associated to set up 
        %reported to "well" above. The parameter "well" define the 
        %"cl" in each point which have well.
        %IMPORTANT: if the mesh is structured ("keymesh" == 'n') 
        %or there is no wells
        %(numwell = 0), this section does not work.
        elseif strcmp(rowgeo(1:5),'Point') && strcmp(keymesh,'n') %??????????????????????????????
            %Fill the matrix "geopoints". This will store the
            %coordinates of the points in order to use below to do
            %compare when a well position is evaluated.
            geopoints(igeopoints,:) = linelabel(2:4);
            %If there is well
            if numwell > 0
                %Rewrite the command line associated to nodes
                %Swept the matrix "well"
                for iwell = 1:numwell    
                    if sum(well(iwell,2) == linelabel(2)) && ...
                            sum(well(iwell,3) == linelabel(3)) 
                        %"clflag" receives either the original value to
                        %"cl1" or the corrected value
                        clflag = well(iwell,1) + 1;
                        %"break" is used to getout of the loop is the 
                        %condition is attended. 
                        break;
                    %If does not exist change in the line command 
                    %evaluated, "clflag" receive the original value 
                    %which is the fiveth value of "linelabel" vector.
                    else 
                        clflag = linelabel(5);
                    end  %End of internal IF
                end  %End of FOR
            %Does not exist well, so in not necessary create a
            %different mesh size associated to node. In this case
            %"clflag" receives its own value.
            else
                clflag = linelabel(5);                        
            end  %End of IF
            %The line is printed with or without the last numeric 
            %value corrected
            fprintf(geometry,'Point(%u) = {%f, %f, %f, cl%u};\r\n',...
                linelabel(1),linelabel(2),linelabel(3),...
                linelabel(4),clflag);
            %Increment "igeopoints"
            igeopoints = igeopoints + 1;
            
        %Verify the command LINE
        elseif strcmp(rowgeo(1:5),'Line(')
            %"storepoint" stores the points which constitute each
            %line evaluated.
            storepoint(iline,:) = linelabel(2:3);
            iline = iline + 1;
            %Print the command line reported to "Lines"
            fprintf(geometry,'%s;\r\n',char(rowgeo));
                
        %Print the PHYSICAL POINT:
        %Print part of command line. The rest will be printed
        %inside of internal for.
        %"numbcnode" is the amount of boundary condition applyed 
        %only to nodes. Its value is obtained in the line 338 
        elseif strcmp(rowgeo(1:2),'cl') ~= 1 && ...
                strcmp(rowgeo(1:10),'Physical P') && ...
                ipoint <= numbcnode
            %Print part of command line. The rest will be printed
            %according with the data present in the row.
            fprintf(geometry,'Physical Point(%u) = {',...
                bcnodepointer(ipoint));
            %Verify in the row evaluated the peace after the "key"
            %So, this peace is added to first part of row
            letterlimit = find(rowgeo == '{');
            fprintf(geometry,sprintf('%s;\r\n',...
                rowgeo(letterlimit + 1:length(rowgeo))));
            ipoint = ipoint + 1;

        %Modify the PHYSICAL LINES:
        %catch the command lines which report a physical line and 
        %the lines which is associated to it.
        %Identify the command line where to be written "Physical L"
        elseif strcmp(rowgeo(1:2),'cl') ~= 1 && ...
                strcmp(rowgeo(1:10),'Physical L') && ...
                iedge <= size(bcflag,1)
            %Print part of command line. The rest will be printed
            %according with the data present in the row.
            fprintf(geometry,'Physical Line(%u) = {',...
                bcflag(iedge,1));
            %Verify in the row evaluated the peace after the "key"
            %So, this peace is added to first part of row
            letterlimit = find(rowgeo == '{');
            fprintf(geometry,sprintf('%s;\r\n',...
                rowgeo(letterlimit+1:length(rowgeo))));
            iedge = iedge + 1;
                    
        %Modify the physical SURFACES:
        %Catch the line reported to "plane surface" and put its
        %label in increasing order to associated to label of each
        %permeability seted in "Start.dat"
        elseif strcmp(rowgeo(1:2),'cl') ~= 1 && ...
                strcmp(rowgeo(1:10),'Physical S') && ...
                iplane <= numlayer
            %Print part of command line. The rest will be printed
            %according with the data present in the row.
            fprintf(geometry,'Physical Surface(%u) = {',iplane);
            %Verify in the row evaluated the peace after the "key"
            %So, this peace is added to first part of row
            letterlimit = find(rowgeo == '{');
            fprintf(geometry,sprintf('%s;\r\n',...
                rowgeo(letterlimit+1:length(rowgeo))));
            %Increment "iplane"
            iplane = iplane + 1;
            
        %Verify LINE LOOP:
        %Verify which "Line Loops" constitute each surface. This
        %option is just relevant when the mesh structured option is set.
        elseif strcmp(rowgeo(1:7),'Line Lo')
            %Fill "getloop". For more then one materials "getloop"
            %will be a matrix
            getloop(igetloop,1:length(linelabel)) = linelabel;
            %Add 1 to "igetloop"
            igetloop = igetloop + 1;
                    
            %After catch the informations write the command line 
            %referred to "Line Loop"
            fprintf(geometry,'%s;\r\n',char(rowgeo));

        %Write the lines which have no change:
        %There is no any one change. The command line are simply 
        %copyed. this occur mainly in both "Point" and "Planar 
        %Surphace" command.
        %Avoid rewrite the lines which were changed.
        else                     
            %Write again the command lines in the new *.geo file 
            %without any intervenction in the command line
            fprintf(geometry,'%s;\r\n',char(rowgeo));
        end  %End of external IF
    end  %End of FOR which swepts each row of *.geo

    %----------------------------------------------------------------------??????????????????????????????????????????????????????????????????????????????
    %Set WELLS out of corners
            
    %Initialize "endline". Stores how many points define the well line
    endline = 0;
    %jump one line in *.geo file
    fprintf(geometry,'\r\n');
    %"contwellin" is a counter which quantify how many wells are inside
    %of damain. The wells markeds in the corners does not came to this
    %parameter.
    if numwell > 0 && strcmp(keymesh,'n')
        %When the well is a POINT:
        contwellin = 1;
        for iwell = 1:numwell
            %Verify if well coordinate (x and y) is among the point
            %coordinate in geo file (points which constitute the 
            %domain). Verify also if the well is by a point or by a
            %line (seventh column of "well" define this)
            if all(ismember(well(iwell,2:3),geopoints(:,1:2))) == 0 
                %Write the command line in the new *.geo to create a
                %new point. This will be used to refine the mesh 
                %surrounding it. 
                fprintf(geometry,'Point(%u) = {%f, %f, 0, cl%u};\r\n',...
                    size(geopoints,1) + contwellin,well(iwell,2),...
                    well(iwell,3),iwell + 1);
                %Write the Physical group reported to point above
                %created.
                fprintf(geometry,'Physical Point(%u) = {%u};\r\n',...
                    1000 + contwellin,size(geopoints,1) + contwellin);
                %Write the last command. This reports the point with
                %its surface. The surface's number is coutch from
                %matrix "well" (first column)
                if well(iwell,9) == 0
                    fprintf(geometry,'Point {%u} In Surface {%u};\r\n',...
                        size(geopoints,1) + contwellin,...
                        getloop(well(iwell,1),1));
                end  %End of IF
                %Increments "endline" if "well(iwell,9)" == 3 
                %(point which define the source line)
                endline(iwell) = ...
                    (size(geopoints,1) + contwellin)*(well(iwell,9) == 3);
                %Add 1 to "contwellin"
                contwellin = contwellin + 1; 
            end  %End of internal IF
        end  %End of FOR

        %When the well is a LINE:      %???????????????????????????   
        %Redifine "endline"
        endline = endline(logical(endline ~= 0));
        %Initialize "contline"
        contline = contwellin;
        %More than one point has flag to build the line source
        if sum(logical(endline ~= 0)) >= 2
            %Write the command line in the new *.geo to create a
            %new LINE. This will be used to refine the mesh 
            %surrounding it. 
            fprintf(geometry,'Line(%u) = {%f, %f};\r\n',...
                size(storepoint,1) + 1,endline(1),endline(2));
            %Write the Physical group reported to point above
            %created.
            fprintf(geometry,'Physical Line(%u) = {%u};\r\n',...
                1000 + contline,size(storepoint,1) + 1);
            %jump one line again
            fprintf(geometry,'\r\n');
            %Write the last command. This reports the LINE with
            %its surface. The surface's number is coutch from
            %matrix "well" (first column)
            fprintf(geometry,'Line {%u} In Surface {%u};\r\n',...
                size(storepoint,1) + 1,...
                getloop(well(iwell,1),1));
        end  %End of IF (well as a line)                    
    end  %End of external IF ("numwell" > 0)
            
    %----------------------------------------------------------------------
    %Set parameters about STRUCTURED mesh building
        
    %Before the *.geo to be closed, some lines are added. These are 
    %related to insertion of command "transfine" in the end of new *.geo
        
    if strcmp(keymesh,'y')
        %jump one line
        fprintf(geometry,'\r\n');
        %Insert the command line referred to structured mesh em 
        %direction "x"
        %Write the command "Transfine" to each mate of line in the
        %direction "x"
        for ix = 1:numparx
            fprintf(geometry,'Transfinite Line {%u,%u} = %u Using Progression %f;\r\n',...
                transflabelx(ix,3),transflabelx(ix,4),transflabelx(ix,1),...
                transflabelx(ix,2));
        end  %End of FOR
        %Write the command "Transfine" to each mate of line in the
        %direction "y"
        for iy = 1:numpary
            %Insert the command line referred to structured mesh em 
            %direction y
            fprintf(geometry,'Transfinite Line {%u,%u} = %u Using Progression %f;\r\n',...
                transflabely(iy,3),transflabely(iy,4),transflabely(iy,1),...
                transflabely(iy,2));
        end  %End of FOR
        %Swept the amount of layers catched of "Start.dat". "numlayer"
        %is also the amount of rows in the "getloop" matrix
        for iloop = 1:numlayer
            %Same initializators:
            ilist = 1;  %Is a counter to be used below
            orderpoints = 0;  %Storages all points in each row evaluated
            looplinelist = 0;  %Stores the definitive sequence of points
            %Put in increasing order the lines which constitute
            %each loop evaluated.
            for jloop = 2:size(getloop,2)
                %Accumulate the points verified in each line 
                orderpoints = ...
                    [orderpoints storepoint(abs(getloop(iloop,jloop)),:)];
                %The third FOR will organize the points in
                %increasing order. 
                for itransfer = 2:length(orderpoints)
                    %This condition avoids repetead number and null
                    %numbers.
                    if orderpoints(itransfer) ~= 0 && ...
                            sum(orderpoints(itransfer) == ...
                            looplinelist) == 0
                        looplinelist(ilist) = ...
                            orderpoints(itransfer);
                        ilist = ilist + 1;
                    end  %End of IF
                end  %End of third FOR
            end  %End of second FOR
 
            looplinelist = sort(looplinelist);
            %Do "transfine in the surface"
            fprintf(geometry,'Transfinite Surface {%u} = {%u,%u,%u,%u};\r\n',...
                getloop(iloop,1),looplinelist(1),looplinelist(2),...
                looplinelist(3),looplinelist(4));
            %Verify if the element is a quadrangle. This information is
            %caught from "Start.dat". Case the element is a quadrangle
            %a new command line is written ("Recombine"). This is
            %associated to surface, whose number is caught from "getloop" 
            %with row of number "getlooprow"
            if keyelem == 1
                %jump one line
                fprintf(geometry,'\r\n');
                %Write the command
                fprintf(geometry,'Recombine Surface {%u};\r\n',...
                    getloop(iloop,1));
            end  %End of IF
        end  %End of FOR
    end  %End of IF

    %Close all text files opened 
    fclose(geometry);
    fclose(options);

    %----------------------------------------------------------------------
    %Execulte gmsh by SYSTEM, creating a mesh without refinement
        
    %Attribute to variable "command" the statment to inicialize gmsh.
    %"sprintf" changes variables in text. This is necessary to use the
    %"system" command.
    %"char()" converts a string in a chain of char
    command = ['gmsh ' mainpath '\' optfile' ' ' mainpath '\' geofile' ...
        ' -2 -o ' mainpath '\' char(mshfile')];
    %Begin gmsh throughout system prompt
    system(command);
    
    %----------------------------------------------------------------------
    %Mesh refinement:
    
    %Building a loop of refinement. The process repeat until "i" reach
    %the number "nref"
        
    %Just came to loop if will there the purpose of refine the mesh ...
    %(where nref > 0)
    if nref > 0
        for i = 1:nref  %Loop to refine
            %Attribute to variable "command" the inicialization with
            %refinement
            command = ['gmsh ' mainpath '\' geofile' ' -2 ' mainpath ...
                '\' char(mshfile)' ' -refine']; 
            %Begin gmsh throughout DOS ptompt
            system(command);
        end  %End of FOR
    end  %End of internal IF
                    
else  %In the case where there is an error in starting parameter
    disp('Incorrect letter about *msh file! Try again');
end  %End of BIG EXTERNAL IF

%--------------------------------------------------------------------------
%Opening of gmsh's user interface
    
%This option open the gmsh file after building mesh with the data above
if strcmp(keyopen,'y') && strcmp(keyfile,'n')
    %Build the command
    command = ['gmsh ' mainpath '\' char(mshfile')];
    %Open gmsh throughout system prompt
    system(command);
%This option open the gmsh file with the mesh that you have
elseif strcmp(keyopen,'y') && strcmp(keyfile,'y')
    %Build the command
    command = ['gmsh ' verifymshfile];
    %Open gmsh throughout system prompt
    system(command);
%Gives a message for the user.    
elseif strcmp(keyopen,'y') ~= 1 && strcmp(keyopen,'n') ~= 1
    disp('Incorrect letter about gmsh opening! Try again');
end  %End of IF

%--------------------------------------------------------------------------
%Guard "Start.dat" in the new folder created (this is done when a folder is
%generated

%Building the "command" to be used in the dos's prompt
command = ['copy start.dat ' char(pathstore)];
%sprintf('copy Start.dat %s\\%s',char(keypath2),foldername);
%Transfer "Start.dat" using the system call command "copy".
system(command);

%--------------------------------------------------------------------------
%IT BUILDS THE DATA STRUCTURE 
%--------------------------------------------------------------------------

%Jump a line in the matlab's prompt
disp(' ');
%Call atention to data to be generated
disp('---------------------------------------------------');
disp('>> It is generating the Data Structure...');
disp(' ');

%--------------------------------------------------------------------------
%"coord" matrix - The coordinates (x,y,z) of each point in a carts. 

[coord,nnode] = getcoord(verifymshfile);
%Verify if we have a "rz" coordinate system
if strcmp(xyrz,'rz')
    coord(:,1) = r0 + coord(:,1);
end  %End of IF

%Gives the information of "coord" generated
disp('"coord" was generated!');

%--------------------------------------------------------------------------
%"elem" matrix - The points that constitute each element
% meshtype- true: multiscale mesh false: fem mesh
[elem,nbe,nelem,nodelim,intnode,flaglim,inboundedge,meshtype,elemloc,npar,coarseelem ,ghostelem] = ...
    getelemfull(verifymshfile,nnode,numwell,well);
    
%Gives the information of "elem" generated
disp('"elem" was generated!');

%--------------------------------------------------------------------------
%"centelem" matrix - It get the centroid coordinate for each control volume 

%Fill the matrix "centelemcoord"
centelem = getcentelem(coord,elem);

%Gives the information of "elem" generated
disp('"centelem" was generated!');

[ elem ] = reordnodelem( elem, centelem, coord );
%
%--------------------------------------------------------------------------
%"elemarea" vector - It get the area for each control volume 

%Fill the vector "elemarea"
%When we have a 'xy' coordinate system
if strcmp(xyrz,'xy')
    elemarea = calcelemarea(coord,elem);
%When we have a 'rz' coordinate system
else
    elemarea = calcelemarea_polarsys(coord,elem,centelem(:,1));
end  %End of IF

%Gives the information of "elem" generated
disp('"elemarea" was generated!');

%--------------------------------------------------------------------------
%"bedge" and "inedge" - It get inform. about boundary and internal edges

[bedge,inedge,klb] = getinfoedgefull(coord,elem,nnode,nbe,nelem,nodelim,...
    intnode,flaglim,verifymshfile,bcflag,meshtype);


% ref1 = inedge(:,3) == 0;
% ref2 = inedge(:,4) == 0;
% ref = ref1 | ref2;
% inedge(ref,:) = [];
% 
% aa = [inedge(ref1,1:3) 201*ones(sum(ref1),1),201*ones(sum(ref1),1) ];
% bb = [inedge(ref2,1:2) inedge(ref2,4) 201*ones(sum(ref2),1),201*ones(sum(ref2),1) ];
% bedge = [bedge; aa ; bb];
% inedge(ref,:) = [];
%Gives the information of "bedge" and "inedge" generated
disp('"bedge" was generated!');
disp('"inedge" was generated!');

%--------------------------------------------------------------------------
%"normals" is a matrix with the normal vectors (normal to left element)

%Get the normals
%When we have a 'xy' coordinate system
if strcmp(xyrz,'xy')
    normals = calcnormals(coord,centelem,bedge,inedge);
%When we have a 'rz' coordinate system
else
    normals = calcnormals_polarsys(coord,centelem,bedge,inedge);
end  %End of IF

%Gives the information of "normals" generated
disp('"normals" was generated!');





%--------------------------------------------------------------------------
%"esurn" vectors - It reports the amount of ELEMENTS surronding each POINT.

%This is constituted by two vectors as a linked list. The first one ...
%("esurn1") contain the number of each element associated with the node ...
%listeds. The second one ("esurn2") contain numbers related to position...
%(on "esurn1" vector) of start of counter.

%Swept the "Elem" matriz from 1 until Nelem (rows) and inside it swept from
%1 until the number of nodes that constitute each element (columns)

%Create and initialize "esurn2"
%nnode + 1 rows and 1 column. Begins from 0 until "nnode"
esurn2 = zeros(nnode + 1,1);   

%"ielem" is a counter of elements (swept the rows of "elem")
for ielem = 1:nelem
    %Get "elemcontents"
    elemcontents = elem(ielem,1:4);
    %Get the amount of edges by element.
    amountedges = sum(elemcontents ~= 0);
    %Get a vector with the vertices for the element evaluated. 
    vertices = elemcontents(1:amountedges);
    %Get the number of vertices for this element.
    numvert = length(vertices);
    %This loop swept all nodes who constitute each element (columns of 
    %"elem" matrix). 
    for inode = 1:numvert    
        %"nodeval" is the receives the number of each node. 
        %The elem(1,1) has nodeval = 1, elem(1,2) has nodeval = 5 
        %(see basic exemple in Marcio's notebook) 
        nodeval = vertices(inode);
        %Add one to position which corresponds to number of node. We
        %have to consider that "esurn2" start with value 0 for a better
        %use as will be saw in the future (see Lohner, 2001 - Chap 2)
        esurn2(nodeval + 1) = esurn2(nodeval + 1) + 1;
    end  %End of FOR (in each element)
end  %End of FOR (in all "elem" matrix)

%"esurnqnt" recives "esurn2" before this be recalculated in order to 
%represent an acumulate value (function "cumsum") 
esurnqnt = esurn2;    
%Create a vector with the amount accumulated from "esurn2" already 
%obtained (esurn2cum(i) = esurn2(i) + esurn2(i-1))
esurn2 = cumsum(esurn2);
esurn2aux = esurn2;

%Create and initialize "esurn1"
%The number of rows of "esurn1" is equal to amount accumulated of elements
%in each node
esurn1 = zeros(esurn2(nnode + 1),1);
%Initialize "nodestore" which control if the node evaluated 
%already was verifyed. If the node in the "elem" was already evaluated 
%inode is added in one.
nodestore = 0;
%"inodestore" is a counter of "nodestore"
inodestore = 1;
%Initialize "pointvertdone"
pointvertdone = zeros(size(coord,1),1);

%Fill "esurn1" (list of elements surrounding each node)
%Swept from 1 until the last element
for ielem = 1:nelem  
    %Get "elemcontents"
    elemcontents = elem(ielem,1:4);
    %Get the amount of edges by element.
    amountedges = sum(elemcontents ~= 0);
    %Get a vector with the vertices for the element evaluated. 
    vertices = elemcontents(1:amountedges);
    %Get the number of vertices for this element.
    numvert = length(vertices);
    %Swept each node of row's "elem" evaluated
    for inode = 1:numvert
        %"nodeval" receives the index of node on "elem" matrix
        nodeval = vertices(inode);

        %Evaluate "pointvertdone"
        if pointvertdone(nodeval) == 0
            %Just one element concors to node evaluated.
            %In general this occur in the domain corner.
            if esurnqnt(nodeval + 1) == 1
                %"store" obtain the value (of "esurn2aux" vector) which 
                %corresponds to position whose number is "nodeval" + 1.
                %Remenbering "esurn2" begins from 0 (for that +1)
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "Elem" row in the "esurn1"
                %whose the position is termined by "strore" index
                esurn1(store) = ielem;
                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - 1;
            %Two or more elements concor to node evaluated and the node 
            %evaluated to be inside domain. In general it is what happen in 
            %the domain.
            elseif esurnqnt(nodeval + 1) > 1 && ...
                    sum(nodeval == nodestore) ~= 1
                %Initialize "elemorder". This parameter receives the
                %elements which try to node evaluated in increasing order.
                elemorder = zeros(1,esurnqnt(nodeval + 1));
                
                %----------------------------------------------------------
                %"nodeval" over boundary
                
                %Verify if the node evaluated ("nodeval") is over the 
                %boundary. In this case we have to find the first element
                %and the last element in order maintain the counter
                %clockwise way (little more complicated).
                rowsinbedge = ...
                    bedge(any(logical(bedge(:,1:2) == nodeval),2),1:3);

                %The node evaluated is inside domain. In this case is
                %easer than the earlier situation once any element may be 
                %the first one. The cycle finishes when the last element 
                %before the first is found.
                if isempty(rowsinbedge)  
                    %Attribute the first element to "elemorder". 
                    %Considering "ielem" contain one of elements evaluated.
                    elemorder(1) = ielem;

                %The node evaluated is over the boundary. The first element
                %in "elemorder" cames from "bedge"
                else
                    %Initialize "elemcandidate". It stores the number of
                    %elements which can be the first element in "elemorder"
                    elemcandidate = rowsinbedge(:,3);
                    %"rowsvertices" receives just the number of vertices
                    rowsvertices = rowsinbedge(:,1:2);
                    %We need verify which other vertices belong to edge
                    %over the boundary.
                    rowsvertvec = reshape(rowsvertices',4,1);
                    %It is used insteady "setdiff".
                    vertoutnodeval = ...
                        rowsvertvec(logical(rowsvertvec ~= nodeval));
                    
                    %Get the corrdinate of "nodeval"
                    nodevalcoord = coord(nodeval,:);
                    %Get coordinate of "vertoutnodeval"
                    otherverexcoord = coord(vertoutnodeval,:);
                    %Evaluate which edge is the first:
                    %Get the vector (candidate to first edge). One vector
                    %by row
                    edgevec = ...
                        otherverexcoord - [nodevalcoord; nodevalcoord];
                    %Get a vector which points to out of domain.
                    pointout = [nodevalcoord; nodevalcoord] - ...
                        centelem(elemcandidate,:);
                    %Get a cross product between the "pointout" and 
                    %"edgevec"
                    %Candidate 1, for example:
                    crossprod = cross(pointout(1,:),edgevec(1,:));
                    %Evaluate the boolean option:
                    booleancross = (crossprod(3) > 0);
                    
                    %Attribute the other vertex (edge over the boundary).
                    vertexout = vertoutnodeval(1)*booleancross + ...
                        vertoutnodeval(2)*(1 - booleancross);
                    %Attribute the first element sharing the edge over the 
                    %boundary.
                    firstelem = elemcandidate(1)*booleancross + ...
                        elemcandidate(2)*(1 - booleancross);

                    %Attribute the first element to "elemorder(1)"
                    elemorder(1) = firstelem;
                end  %End of internal IF ("bedge" or "inedge")

                %Create a vector with number of "inedge" rows
                iinrow = 1:size(inedge,1);

                %Find the last element to be stored in the "elemorder"
                for iorder = 2:esurnqnt(nodeval + 1)
                    %Obtain the row or rows of "inedge" whose "nodeval" 
                    %belongs to first two columns of "inedge" and, to same 
                    %time, the aforefind element belongs to two last 
                    %columns of "inedge".
                    %Evaluate the vertices and elements in "inedge":
                    evalvertelem = all([any(inedge(:,1:2) == nodeval,2) ...
                        any(inedge(:,3:4) == elemorder(iorder-1),2)],2);
                    %Get the "inedge" row.
                    rowinedge = iinrow(evalvertelem)';

                    %Get the size of "rowinedge"
                    rowinedgesize = length(rowinedge);
                    %Initialize "rowpointer"
                    rowpointer = zeros(1,rowinedgesize);
                    %This loop points to entity of "rowinedge" which
                    %satisfy the necessary conditions. That row which
                    %satisfies receive the value "1", while that one 
                    %which does not satisfy receive "0"
                    for irow = 1:rowinedgesize
                        rowpointer(irow) = ...
                            (inedge(rowinedge(irow),2) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),3)) | ...
                            (inedge(rowinedge(irow),1) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),4));
                    end  %End of internal FOR

                    %"rowdef" is a definitive row which satisfy the
                    %necessary conditions aforesaid. That is the 
                    %position which have "1" insteady "0" 
                    rowdef = logical(rowpointer);   
                    %"rowinedge" receives, finaly, the unic row which 
                    %satisfy the conditions above.
                    rowinedge = rowinedge(rowdef);
                    
                    %A decision is made as function of "nodeval"'s
                    %value. The commands written below avoids "IF"  
                    %"nodeval" is lower than other node which 
                    %constitute the common edge among two elements. 
                    %In this case, a left element is attributed to 
                    %"elemorder" just if "nodeval" is lower than other 
                    %node   
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,3)*any((nodeval ~= ...
                        max(inedge(rowinedge,1:2))));
                    %In this case, a right element is attributed to 
                    %"elemorder" just if "nodeval" is major than other 
                    %element   
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,4)*any((nodeval == ...
                        max(inedge(rowinedge,1:2))));
                end  %End of FOR (swept the elements surrounding "nodeval")
                
                %"nodestore" receives the node evaluated. Thus, when this
                %will evaluated in the future that will be by passed,
                %because was already (see line 1241)
                nodestore(inodestore) = nodeval;
                %Add 1 to "inodestore"
                inodestore = inodestore + 1;

                %Put the number in its places in the "esurn1"
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "elem" row in the 
                %"esurn1" whose the position is termined by "strore" index
                
                %"retrocount" is a retroative counter.
                retrocount = esurnqnt(nodeval + 1);
                for iretro = store:-1:store - (esurnqnt(nodeval+1) - 1)
                    esurn1(iretro) = elemorder(retrocount);
                    %backward the counter
                    retrocount = retrocount - 1;
                end  %End of FOR

                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - ...
                    esurnqnt(nodeval + 1);
            end  %End of external IF
            %Turns on "pointvertdone" position
            pointvertdone(nodeval) = 1;
        end  %End of IF ("pointvertdone")
    end  %End of FOR (read each column of each row of "elem")
end  %End of FOR (Swept each row of "elem")

%It clears the vector "iinrow" and another
clear iinrow pointvertdone esurnqnt esurn2aux;


%Gives the information of "esurn" generated
disp('"esurn" was generated!');

%debugCon

%--------------------------------------------------------------------------
%"nsurn" vectors - Report the amount of NODES surronding each NODE.

%We will use the vectors "eserp1" and "esurn2" to swept the elements around
%each point. In the evaluated element we must find the near points to point
%who corresponds to position of "Esurn2" (the first position of this vector 
%corresponds to point 1 and consecutively)

%Initialization of variables
%Works as a counter to "nsurn1" vector
store = 0;
%"nsurn2" points to "nsurn1" and denote the begin and end of node 
%surrounding cycle of each node
nsurn2 = zeros(nnode + 1,1);

%Swept the vector "esurn2" to pass for each point (first for). Each one 
%index "i" is a position of aforesaid (abovementioned) vector.
%With the development of "i", points already visited will not any more  
    for i = 1:nnode
        %Swept from position "i" add 1 until position "i+1" (end of 
        %elements surrounding point cycle).
        %"iesurn" is a counter of thous elements surrounding each point

        %Initialization of variable
        %"nsurn1aux" is an auxiliar vector which avoid repeated node 
        %be stored in "nsurn1" (see below)
        nsurn1aux = 0;
        %"insurn" is a counter who will be used in "nsurn1aux" vector
        insurn = 0;
        %Swept all elements surround each node evaluated
        for iesurn = esurn2(i) + 1:esurn2(i + 1)
            %"ielem" will count which elements be in the cycle.
            %This parameter receives the elements 1,2 & 3 (whose positions 
            %are 7,8 and 9) - see basic example in Marcio's notebook
            ielem = esurn1(iesurn);   

            %Initializing "inode" (node counter in each element evaluated)
            inode = 1;
                %Verify in each element (use the "elem" matrix) which
                %nodes are linked with the one evaluated.
                %This happen while each node evaluated is different of zero 
                %and if that node does not belong to 5th column of "elem" 
                %matrix 

                while inode < 5 && elem(ielem,inode) > 0 
                    
                    %The decision command store the value from "elem" in 
                    %"nsurn1" if and only if the node evaluated does not be
                    %the node "i" (becouse we are searching the neighbour 
                    %nodes to node "i", not the own node "i"). Another
                    %condition is that the neighbour node already stored
                    %does not be stored again
                    %This decision command is exclusive to TRIANGLES
                    if elem(ielem,4) == 0 & elem(ielem,inode) ~= i & ...
                            elem(ielem,inode) ~= nsurn1aux 
                        %"store" stores the position in "Nsurn1" vector
                        store = store + 1;
                        %Add 1 to counter "insurn". When "i" changes its 
                        %value this parameter is nulled
                        insurn = insurn + 1;
                        %"Nsurn1" receives the node who attends to above
                        %conditions
                        nsurn1(store) = elem(ielem,inode);
                        %"Nsurn1aux" receives the nodes in order to those
                        %do not be repeated
                        nsurn1aux(insurn) = elem(ielem,inode);
                    %The next decisions command is exclusive to QUADRAN.
                    %The neighbour to node evaluated will be considerate if
                    %those are either forward or backward positioned with 
                    %regard to evaluated node.  
                    elseif elem(ielem,4) ~= 0 && elem(ielem,inode) ~= i
                        %A special case where the first node evaluated into 
                        %row "ielem" is different of value "i" but does not 
                        %is neighbour of node "i". This happen only when
                        %the node evaluated is the first and the node who
                        %corresponds to "i" is the third.
                        if inode == 1 && elem(ielem,3) == i
                            inode = inode + 1;
                        end  %End of first internal IF
                            %To avoid which node already accounted be 
                            %stored again this decision command just 
                            %storage nodes who does not be in Nsurn1aux 
                            %(auxiliar vector)
                            if elem(ielem,inode) ~= nsurn1aux
                                store = store + 1;
                                insurn = insurn + 1;
                                nsurn1(store) = elem(ielem,inode);
                                nsurn1aux(insurn) = elem(ielem,inode);
                                %Avoid which the neighbour node be stored, 
                                %once it does not is straightly linked with 
                                %node evaluated
                            end  %End of second internal IF
                            
                            %If the node evaluated was considered neighbour
                            % of "i" (and stored) the counter "inode"
                            %must jump one node because the next one
                            %certaly does not will neighbour 
                            inode = inode + 1;
                    end  %End of external IF
                    inode = inode + 1;  %Evaluate the next node
                
                 
                
                end  %End of WHILE
        
        end  %End of second FOR
        %Update "nsurn2". This variable receives the number of increments
        %that "store" had
        nsurn2(i + 1) = store;
    end  %End of first FOR

%Reorder "nsurn1" (counterclockwise). It returns a column vector.
nsurn1 = reordernsurn(esurn1,esurn2,nsurn1,nsurn2,bedge,coord,elem);

%It clears "nsurn1aux"
clear nsurn1aux;

%Gives the information of "nsurn" generated
disp('"nsurn" was generated!');

%% teste ordenamento

%debugCon

%% reoordering elements
%[centelem, elem, elemarea, esurn1, esurn2, inedge, bedge, elemloc,  coarseelem,wells] = reordering(centelem, elem, elemarea, esurn1, esurn2, inedge, bedge, elemloc,  coarseelem, npar,coord,wells);

%--------------------------------------------------------------------------
%Get the position of "bedge" and "inedge" rows.

%It has a structure similar to "nsurn1"
%rowposit = getrowposition(bedge,inedge,nsurn1,nsurn2);
%% remarcar aqui
% nsurn1 = [];
% nsurn2 = [];

%--------------------------------------------------------------------------
%Generate "esureface" (face neighbor element) and "esurefull" (vertices and 
%face) neighbor

%% DESABILITAR PARA TESTAR
% [esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,...
%     size(bedge,1),inedge,esurn1,esurn2,nsurn2,rowposit);

%[esureface1,esureface2,esurefull1,esurefull2] = [ 0,0,0,0]; 
 esureface1 = [];
 esureface2 = [];
 esurefull1 = [];
 esurefull2 = [];
%%
%Gives the information of "nsurn" generated
disp('"esure" were generated!');

%--------------------------------------------------------------------------
%porosity distribution - Create a vector with a map of porosity. The data
%could be obtained from a external file or data came from "Start.dat"

%Case the media has a constant porosity, this will attributed to each layer
%which the domain has.
if strcmp(keypor,'n')
    %Open the viscosity file
    readpor = ...
        fopen(sprintf('%s\\%s%s\\%s',char(keypath2),char(basename),...
            char(complemfile),char(porname)));
    %"kmap" reads the file put inside work folder
    getpordata = textscan(readpor,'%n%n%*[^\n]',...
        'delimiter',' ');
    %Attribute to "pormap" the porosity map
    pormap = cell2mat(getpordata);
        
    %Close the file
    fclose(readpor);
end  %End of IF
 
%--------------------------------------------------------------------------
%dens - Read of "start.dat" the amount of fluids considered and the
%respective densities of each one

%The density varies strongly with the pressure field
if strcmp(keydens,'y')
    %Open the viscosity file
    readens = ...
        fopen(sprintf('%s\\%s%s\\%s',char(keypath2),char(basename),...
            char(complemfile),char(densname)));
    %"kmap" reads the file puted inside work folder
    getdensdata = textscan(readens,'%n%n%*[^\n]',...
        'delimiter',' ');
    %Attribute to "kmap" the permeability map
    dens = cell2mat(getdensdata);
        
    %Close the file
    fclose(readens);
end  %End of IF

%--------------------------------------------------------------------------
%visc - Read of "start.dat" the amount of fluids considered and the
%respective viscosities of each one

%The viscosity varies strongly with the pressure field
if strcmp(keyvisc,'y')
    %Open the viscosity file
    readvisc = ...
        fopen(sprintf('%s\\%s%s\\%s',char(keypath2),char(basename),...
            char(complemfile),char(viscname)));
    %"kmap" reads the file puted inside work folder
    getviscdata = textscan(readvisc,'%n%n%*[^\n]',...
        'delimiter',' ');
    %Attribute to "kmap" the permeability map
    visc = cell2mat(getviscdata);
        
    %Close the file
    fclose(readvisc);
end  %End of IF

%--------------------------------------------------------------------------
%kmap - Create a vector with a map of permeability distribution whose data
%are obtained from "Start.dat" file

%"strcmp" compare the string readed with the letters "y", "r" of "n"

%[r] The media is hightly heterogeneous but its value will be generated
%from random command. To make this is necessary to enter with a kmax 
%value who will be read from "start.dat" 
if strcmp(keyperm,'r')
    %Initialize "kmap"
    kmap = zeros(size(elem,1),5);
    %"randomdist" has one permeability tensor (four components) to each
    %element with values variating between "kmin" and "kmax"
    randomdist = randi(klimits,[nelem,3]);
    %Organize "ramdomdist" by encreasing order (just into columns)
    randomdist = sort(randomdist,2);
    %Fill "kmap"
    %Create a counter from "1" untill the element's number
    ik = 1:size(elem,1);
    %Fill the first column of "kmap" and change the last column of 
    %"elem". Each element has a different permeability
    kmap(:,1) = ik;
    elem(:,5) = ik;
    %Fill the other columns, according elipsity condition
    kmap(:,2:5) = [randomdist(:,3) randomdist(:,1) randomdist(:,1) ...
        randomdist(:,2)];

%If "keyperm" receives [y] each element has a permeability tensor. 
%Its values will be read from a file like-table (excel or other) 
%whose name is read below:
elseif strcmp(keyperm,'y')
    %The permeability is obtaind by *.dat file. 
    if strcmp(hipermname,'0') == 0
        %"char()" converts the variable "hipermname" in a chain of char.
        %Open the permeability file
        readperm = ...
            fopen(sprintf('%s\\%s%s\\%s',char(keypath2),char(basename),...
            char(complemfile),char(hipermname)));
        %"kmap" reads the file put inside work folder
        getpermdata = textscan(readperm,'%n%n%n%n%n%*[^\n]');
        %Attribute to "kmap" the permeability map
        kmap = cell2mat(getpermdata);

        %Close the file
        fclose(readperm);
    end  %End of IF
    %The last column of "elem" is changed due highly range of heterogenity 
    %value. Create a counter from "1" untill the element's number.
    ik = 1:size(elem,1);
    %Change the last column of "elem" 
    elem(:,5) = ik;
end  %End of IF

%Message of preprocessor.
disp('physical properties were generated!');

%--------------------------------------------------------------------------
%wells - Gives details about well(s) set on the domain 

%Look to line 372. There the data referred to wells (written in the 
%"Start.dat") are read and used in this section.

%Initialize "wells" and "countwell". 
wells = 0;
countwell = 0;

%"j" is a counter of wells row.
j = 0;
%When the WELL is a POINT (Unstructured mesh with mesh adaptation)
%Find the nodes which contain the wells. These are related with rows of
%matrix "coord".
if numwell > 0 && any(well(:,9) == 0)
    %Swept the amount of wells
    for iwell = 1:numwell
        %"centwell" is a vector with the node which stores the center 
        %of well set. Surrounding these nodes the elements will have 
        %sourse terms.
            centwell(iwell) = find(coord(:,1) == well(iwell,2) & ...
            coord(:,2) == well(iwell,3));
        %Once the node was finded we can verify which elements to be
        %surrounding it. This information will be stored in the first
        %column of "wells matrix"
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,1) = esurn1(esurn2(centwell(iwell))+1:...
            esurn2(centwell(iwell)+1));
        %The second column will store the well's number. This order
        %follow the order of write in the "Start.dat". "iwell" gives
        %this number.
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,2) = iwell;
        %The third column stores the saturation flag. This may have 
        %values from [301 - 400]. Is important to remember that to 
        %productor well this parameter receives "0".             
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,3) = well(iwell,5);
        %The fourth column stores the saturation value. This will be "0"
        %for producer wells.
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,4) = well(iwell,6);
        %The fiveth column stores the pressure flag. That may be have 
        %values among [401 - 500] to injector well and [501 - 600] to 
        %productor well  
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,5) = well(iwell,7);
        %The sixth column stores the pressure value. That may receives the
        %flow rate value if the pressure flag is "0"
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,6) = well(iwell,8);
        %"countwell" receives a new value
        countwell = countwell + ...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell)));
    end  %End of FOR
end  %End of IF

%When the WELL is a POINT (STRUCTURED mesh with NO mesh adaptation)
if numwell > 0 && any(well(:,9) == 1)
    %Verify which ones row of "well" has sign related to point well in 
    %structured mesh (1)
    numstrwell = find(well(:,9) == 1);
    %Cont the number of wells
    countwell = countwell + 1;
    %Swept the wells alocated.
    for istrwell = 1:length(numstrwell)
        %Verify in "coord" the coordinate "x" of well center
        fstnodex = logical(coord(:,1) == well(numstrwell(istrwell),2)); 
        %Verify in "coord" the coordinate "y" of well center
        fstnodey = logical(coord(:,2) == well(numstrwell(istrwell),3)); 
        %Obtain the row of "coord" where the coordinate of "x" and "y" are
        %coincident. It is the number of node where the well is centered.
        wellnode = find((fstnodex + fstnodey) == 2);
        %Verify the elements which surrounding this node.
        %"poselem" guards the position in "esurn1" of elements which
        %surrounding the node  "wellnode"
        poselem = esurn2(wellnode) + 1:esurn2(wellnode + 1);

        %Fill the first column of "wells" (elements)
        wells(j + 1:j + length(poselem),1) = esurn1(poselem);
        %Fill the second column of "wells" (number of well)
        wells(j + 1:j + length(poselem),2) = countwell;
        %Fill the thrid column of "wells" (flag of saturation)
        wells(j + 1:j + length(poselem),3) = well(numstrwell(istrwell),5);
        %Fill the fourth column of "wells" (value of saturation)
        wells(j + 1:j + length(poselem),4) = well(numstrwell(istrwell),6);
        %Fill the fiveth column of "wells" (flag of pressure)
        wells(j + 1:j + length(poselem),5) = well(numstrwell(istrwell),7);
        %Fill the sixth column of "wells" (value of pressure)
        wells(j + 1:j + length(poselem),6) = well(numstrwell(istrwell),8);
        %Increment "j"
        j = j + length(poselem);
        %Increment "countwell" (when there is other well)
        countwell = countwell + 1;
    end  %End of FOR        
end  %End of IF

%BUCKLEY-LEVERETT Applications.
%When the WELL is a LINE and this LINE is over BOUNDARY
if numwell > 0 && any(well(:,9) == 2)
    %Catch the elements associated with edge defined.

    %1. Injector well:
    %Verify if there exists injector well
    pointinjec = logical(well(:,5) > 300 & well(:,9) == 2);
    %It exists:
    if any(pointinjec)
        %Count the number of wells
        countwell = countwell + 1;
        %Get the elements associated to injector well (in the face)
        %Verify "bcflag"
        getflag = logical(bcflag(:,1) > 200 & bcflag(:,1) < 300 & ...
            bcflag(:,2) > 0);
        %Define the flag for injector face
        injecflag = bcflag(getflag,1);
        %It points to "bedge" rows with "injecflag"
        pointflag = logical(bedge(:,5) == injecflag);
        %Get the elements for the injector well
        injecelem = bedge(pointflag,3);
        %Get saturation flag
        pointsatflag = logical(well(:,5) > 300 & well(:,5) < 400);
        getsatflag = unique(well(pointsatflag,5)); 
        getsatval = unique(well(pointsatflag,6));
        %Flag an value (pressure)
        getpressflag = unique(well(pointsatflag,7)); 
        getflowrateval = unique(well(pointsatflag,8));
    
        %Fill the first column of "wells" (elements)
        wells(j + 1:j + length(injecelem),1) = injecelem;
        %Fill the second column of "wells" (number of well)
        wells(j + 1:j + length(injecelem),2) = countwell;
        %Fill the thrid column of "wells" (flag of saturation)
        wells(j + 1:j + length(injecelem),3) = getsatflag;
        %Fill the fourth column of "wells" (value of saturation)
        wells(j + 1:j + length(injecelem),4) = getsatval;
        %Fill the fiveth column of "wells" (flag of pressure)
        wells(j + 1:j + length(injecelem),5) = getpressflag;
        %Fill the sixth column of "wells" (value of pressure)
        wells(j + 1:j + length(injecelem),6) = getflowrateval;
    
        %Increment "j"
        j = j + length(injecelem); 
    end  %End of IF (injector well)
    
    %2. Producer well:
    %Verify if there exists producer well
    pointproduc = logical(well(:,7) > 500 & well(:,9) == 2);
    %It exists:
    if any(pointproduc)
        %Count the number of wells
        countwell = countwell + 1;
        %Get the elements associated to Producer well (in the face)
        %Verify "bcflag"
        getflag = logical(bcflag(:,1) > 100 & bcflag(:,1) < 200);
        %Define the flag for producer face
        injecflag = bcflag(getflag,1);
        %It points to "bedge" rows with "injecflag"
        pointflag = logical(bedge(:,5) == injecflag);
        %Get the elements for the injector well
        producelem = bedge(pointflag,3);
        %Flag an value (pressure)
        pointpressflag = logical(well(:,7) > 500 & well(:,7) < 600); 
        getpressflag = unique(well(pointpressflag,7)); 
        getflowrateval = unique(well(pointpressflag,8));
        %Flag an value (saturation)
        getsatflag = unique(well(pointpressflag,5)); 
        getsatval = unique(well(pointpressflag,6));
    
        %Fill the first column of "wells" (elements)
        wells(j + 1:j + length(producelem),1) = producelem;
        %Fill the second column of "wells" (number of well)
        wells(j + 1:j + length(producelem),2) = countwell;
        %Fill the thrid column of "wells" (flag of saturation)
        wells(j + 1:j + length(producelem),3) = getsatflag;
        %Fill the fourth column of "wells" (value of saturation)
        wells(j + 1:j + length(producelem),4) = getsatval;
        %Fill the fiveth column of "wells" (flag of pressure)
        wells(j + 1:j + length(producelem),5) = getpressflag;
        %Fill the sixth column of "wells" (value of pressure)
        wells(j + 1:j + length(producelem),6) = getflowrateval;
    
        %Increment "j"
        j = j + length(producelem);
    end  %End of IF (producer well)
end  %End of IF

%When the WELL is a LINE 
%There is a source line INSIDE domain
if numwell > 0 && any(well(:,9) == 3) 
    %Define the well number
    wellnumber = countwell + 1;
    for i = 1:size(inboundedge,1)
        %Find the row of "inedge" which contain the two nodes which
        %constitute the line source
        indexedge = logical(inedge(:,1) == inboundedge(i,1) & ...
            inedge(:,2) == inboundedge(i,2));
        %Once the edge was finded we can verify which elements to be
        %surrounding it. This information will be stored in the first
        %column of "wells matrix"
        %Attribute the lft and right elements
        wells(countwell + i:countwell + i + 1,1) = ...
            [inedge(indexedge,3); inedge(indexedge,4)];
        %The second column will store the well's number. This order
        %follow the order of write in the "Start.dat". "iwell" gives
        %this number.
        wells(countwell + i:countwell + i + 1,2) = wellnumber;
        %The third column stores the saturation flag. This may have 
        %values from [301 - 400].              
        %Find a row of "well" which contain 1 in the seventh column.
        wellrow = logical(well(:,9) == 3);
        %Fill the third column of "wells" (flag of saturation)
        wells(countwell + i:countwell + i + 1,3) = well(wellrow(1),5);
        %Fill the fourth column of "wells" (value of saturation)
        wells(countwell + i:countwell + i + 1,4) = well(wellrow(1),6);
        %The fiveth column stores the pressure flag. That may be have 
        %values among [401 - 500] to injector well and [501 - 600] to 
        %productor well  
        wells(countwell + i:countwell + i + 1,5) = well(wellrow(1),7);
        %Fill the sixth column of "wells" (value of pressure)
        wells(countwell + i:countwell + i + 1,6) = well(wellrow(1),8);
        %Increment "countwell"
        countwell = countwell + 1;
    end  %End of FOR
end  %End of IF

%The domain has no wells
if numwell == 0
    wells = 0;
end  %End of IF

%Message of preprocessor.
disp('wells properties were generated!');
%pre allocating space for multiscale data 
%ativar aqui

%Jump a line
disp(' ');
%Final message
disp('>> All data were generated with success!!!');



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "getdatafile"
%--------------------------------------------------------------------------

function [head,benchkey,pathstore,numcase,complemfile,phasekey,xyrz,r0,...
    symaxe,keymsfv,coarseratio,pmethod,auxcvfactor,interptype,nonlinparam,...
    smethod,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,limiterflag,keyfile,keypath1,mshfile,keypath2,basename,keymesh,...
    numparx,transflabelx,numpary,transflabely,stdclflag,maxelemsize,...
    minelemsize,keyelem,keyalgor,keyrefine,nref,scalemesh,keyopen,keydens,...
    densname,dens,keyvisc,viscname,visc,satlimit,keyperm,klimits,...
    hipermname,numlayer,kmap,keypor,pormap,porname,keygravity,g,keycapil,...
    ncaplcorey,bcflag,numbcnode,bcnodepointer,courant,totaltime,numwell,...
    well,multiCC] = getdatafile
%It open the data file ("Start.dat")
readfile = fopen('start.dat');

%Write the head in this file. This command colect 8 lines of "Start.dat" 
%and copy it to file which will be created.
getdata = textscan(readfile,'%s',8,'delimiter','\n');
%Attribute the data to "head"
head = getdata{1};

%Ster to reading the parameters from the "Start.dat"
%//Do you want to run a benchmark case? [y] ==> "yes"; [n] ==> "no"  

%--------------------------------------------------------------------------
%Set Benchmark

%1. It reads the answer about to solve or not a benchmark case:
benchkey = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%It reads the path where the data are: 
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "pathstore"
pathstore = getdata{1};
pathstorechar = char(pathstore);

%Get the number of case:
%It uses "isstrprop" to identify who is number or digit.
binpath = isstrprop(pathstorechar,'digit');
%Verify if exist numbers among the last five digits in the patch given:
if (any(binpath(length(binpath) - 5:length(binpath))) == 0) && ...
        (strcmp(benchkey,'y') || strcmp(benchkey,'r'))
    %Send message to user:
    disp(' ');
    disp('---------------------------------------------------------------');
    disp('---------------------------------------------------------------');
    disp('>>> WARNING!');
    disp('---------------------------------------------------------------');
    disp('---------------------------------------------------------------');
    disp('You have set the option "Solve a Benchmark Case"!');
    disp('There is NO any number in the folder''s name that you have set.');
    disp('---------------------------------------------------------------');
    %The user attribute the number of case:
    numcase = input('Please, set below the number of the case that you want to run: \n>>> ');
    
%There is some number inside the path given
elseif (any(binpath(length(binpath) - 5:length(binpath)))) && ...
        (strcmp(benchkey,'y') || strcmp(benchkey,'r'))
    %Get the last digits of "binpath"
    fivedigit = binpath(length(binpath) - 5 + 1:length(binpath));
    %Get how many are "0"
    pointzero = find(fivedigit == 0);
    %Get how many are "1"
    pointone = find(fivedigit == 1);
    %Get the letters of "pathstore"
    getletter = pathstorechar(length(binpath) - 5 + 1:length(binpath));
    %Get only the numbers
    getonlynumb = str2double(getletter(fivedigit));
    %The number has fraction. It needs divide per 10.
    if pointzero(length(pointzero)) > pointone(1)
        %Attribute the number to "numcase"
        numcase = getonlynumb/10;
    %The number has NO fraction. It is integer.
    else
        %Attribute the number to "numcase"
        numcase = getonlynumb;    
    end  %End of IF
%We does not need of any benchmark number
else
    %Attribute null value to "numcase"
    numcase = 0;
end  %End of IF

%2. "complemfile" store the complement name or number of aforesaid file.
%It will be the name of the folder inside main folder
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "complemfile"
complemfile = getdata{1};

%--------------------------------------------------------------------------
%Solver Properties

%3. Type of Simulation. It chooses among Hyperbolic Equation [0], ONE-PHASE 
%problem [1] and TWO-PHASE problem [2]
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "phasekey"
phasekey = getdata{1};

%***OBS.:
%Add "30" to "numcase" when a two-phase flow case is resolved.
if phasekey > 1
    %Correct "numcase"
    numcase = numcase + 30;
end  %End of IF

%Coordinate System:
%4. Coordinate System (only 2D application)
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "xyrz"
xyrz = getdata{1};

%5. Read the initial radius "r0" and the symmetry axe "symaxe" 
%Get radius
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "r0"
r0 = getdata{1};
%Get symmetry axe
getdata = textscan(readfile,'%f %f %f',1,'CommentStyle',{'//'});
%Attribute the data to "symaxe"
symaxe = cell2mat(getdata);
% TO BE MODIFIED
%ghostelem Properties:
%6. Do you want to run Multi-Scale Finite Volume Method? It reads the
%answer
getdata = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Attribute the data to "keymsfv"
keymsfv = getdata{1};
%7. Get Coarsening Ratio
getdata = textscan(readfile,'%u %u',1,'CommentStyle',{'//'});
%Attribute the data to "coarseratio"
coarseratio = cell2mat(getdata);

%Get multiscale localization condition
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
multiCC = cell2mat(getdata);


%Type of Solver for the Eliptic Equation (Pressure Equation)
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "pmethod"
pmethod = getdata{1};

%The auxiliary control volume factor.
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "auxcvfactor"
auxcvfactor = getdata{1};

%9. Type of Nodal Interpolation (for pressure Equation)
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "interptype"
interptype = getdata{1};

%10. Tolerance for the non-linear system solution
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "nltol"
nltol = getdata{1};

%11. Maximum Number of Iteraction allowed for the non-linear system solut.
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "maxiter"
maxiter = getdata{1};

%Attribute "nltol" and "maxiter" to "nonlinparam"
nonlinparam = [nltol maxiter];

%12. Type of Solver for the Hyperbolic Equation (Saturation Equation)
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "smethod"
smethod = getdata{1};

%13. Truly MultiDimensional Parameters (TMU, SMU or SMU2)
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "multdopt"
multdopt(1) = getdata{1};

%14. Where the multidimensionality occurs
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "multdopt"
multdopt(2) = getdata{1};

%15. If you chose [goef], write below the value for the weight "w" and
%the parameter "ni". It must be written as a row vector ("w" "ni"):
getdata = textscan(readfile,'%f %f',1,'CommentStyle',{'//'});
%Attribute the data to "goefopt"
goefreeopt = cell2mat(getdata); 

%Orders for the Approximations (space and time):
%16. The Order for the Hyperbolic Operator approximation
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "order"
order = getdata{1};

%17. The type of approximation for the time term:
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "timeorder"
timeorder = getdata{1};

%Strategies for the Variable Reconstruction
%18. Choose which type of recovery strategy do you want:
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "recovtype"
recovtype = getdata{1};

%Define the type of vicinity for the Least Squares reconstruction:
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "lsneightype"
lsneightype = getdata{1};

%19. Set a value for the inverse distance exponent.
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "lsexp"
lsexp = getdata{1};

%----------------------------------------
%Limiting Strategies and Greek Correction

%20. Edge-based Limiter:
getdata = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Attribute the data to "key_eblim"
key_eblim = getdata{1};
%According "key_eblim" turn on the "eblimkey"
%There is edge based limiter
if strcmp(key_eblim,'y') && order > 1 && strcmp(recovtype,'ggg') == 0
    %1.1. Define if there is a edge-based limiter:
    eblimkey = 'on';
%There is NO edge based limiter
else
    eblimkey = 'off';
end  %End of IF

%21. choose among the Edge-Based limiters available
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "eblimtype"
eblimtype = getdata{1};

%22. Gradient Interpoltion parameter "k" (Blazek, pp. 167). 
%It can receives "1", "0", "-1", "1/3"
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "kint"
kint = getdata{1};

%23. Patch-based Limiter
getdata = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Attribute the data to "key_cvblim"
key_cvblim = getdata{1};
%According "key_cvblim" turn on the "cvblimkey"
%There is edge based limiter
if strcmp(key_cvblim,'y') && order > 1 && strcmp(recovtype,'ggg') == 0
    %1.1. Define if there is a edge-based limiter:
    cvblimkey = 'on';
%There is NO edge based limiter
else
    cvblimkey = 'off';
end  %End of IF

%Choose among the Patch-Based limiters below:
%24. choose among the Edge-Based limiters available
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "cvblimtype"
cvblimtype = getdata{1};

%25. If you chose [wf], write below the additional parameter ("delta"):
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "delta"
delta = getdata{1};

%-----------------
%Greek Correction:

%26. Turn the Greek Correction on:
getdata = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Attribute the data to "keygc"
keygc = getdata{1};
%According "keygc" turn on the "gckey"
%There is Greek Correction
if strcmp(keygc,'y') && order > 1 && strcmp(recovtype,'ggg') == 0
    %1.1. Define if there is Greek Correction:
    gckey = 'on';
%There is NO Greek Correction
else
    gckey = 'off';
end  %End of IF

%-----------------
%"MOOD" Definition

%27. Do you want to turn the MOOD on?
getdata = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Attribute the data to "moodkey"
moodkey = getdata{1};
%Is MOOD turned on?
if strcmp(moodkey,'y') && order > 1 && strcmp(recovtype,'ggg') == 0
    %Define if MOOD is turned on:
    mood = 'on';
%There is NO Greek Correction
else
    mood = 'off';
end  %End of IF

%28. Type of the MOOD strategy
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "moodtype"
moodtype = getdata{1};

%Construct "limiterflag" according options above:
limiterflag = {eblimkey char(eblimtype) cvblimkey char(cvblimtype) delta ...
    gckey kint mood moodtype};

%--------------------------------------------------------------------------
%Geometry and Mesh Properties

%The parameter "keyfile" is a string type which defines if the user has or 
%not a *.msh file
keyfile = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%When the mesh exists a path to the mesh is read.
%"keypath1" will store this path (There is a mesh file).
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "keypath1"
keypath1 = getdata{1};

%"meshfile" is a string parameter which receives a name of *.msh file. 
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "mshfile"
mshfile = getdata{1};

%When the mesh does not exists other path is read
%"keypath2" will store this path (There is no mesh file).
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "keypath2"
keypath2 = getdata{1};

%Get the BASE name ("basename") and a label for complement ("complemfile"). 
%The combination of these will gives the folder name to files which will be 
%created.
%"basename" store the name attributed to *.geo file.
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "basename"
basename = getdata{1};         

%--------------------------------------------------------------------------
%General Options (Configure gmsh)

%----------------
%Mesh Parameters:

%"keymesh" is a parameter which set either structured or unstructured mesh 
keymesh = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%------------------------------
%Read the STRUCTURED MESH data:

%"numparx" is the number of parameters to be set in direction "x". This 
%parameters are read including when structured mesh will not used due 
%reading of "Start.dat".
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numparx"
numparx = getdata{1};

%"transflabelx" catches the parameter reported to "x" direction in order to 
%building a structured mesh in a regular domain 
getdata = textscan(readfile,'%f %f %f %f',numparx,'CommentStyle',{'//'});
%Attribute the data to "transflabelx" 
transflabelx = cell2mat(getdata);

%It adds "1" to first position of "transflabelx" because "gmsh" works with 
%the number of nodes (5, for example), but the user put the number of cntrl
%volumes (4, for example).
transflabelx(:,1) = transflabelx(:,1) + 1;

%"numpary" is the number of parameters to be set in direction "x". This 
%parameters are read including when structured mesh will not used due 
%reading of "Start.dat".
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numpary"
numpary = getdata{1};

%"transflabely" catches the parameter reported to "y" direction in order to 
%building a structured mesh in a regular domain 
getdata = textscan(readfile,'%f %f %f %f',numpary,'CommentStyle',{'//'});
%Attribute the data to "transflabelx" 
transflabely = cell2mat(getdata);
%It adds "1" to first position of "transflabely" because "gmsh" works with 
%the number of nodes (5, for example), but the user put the number of cntrl
%volumes (4, for example).
transflabely(:,1) = transflabely(:,1) + 1;

%--------------------------------
%Read the UNSTRUCTURED MESH data:

%"stdclflag" receives the half size mesh and rewrite it in the new "*.geo"
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "stdclflag"
stdclflag = getdata{1};

%MAX ELEMENT SIZE ("maxelemsize"):
%The defalt value recomended in user manual is 1e22.
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "maxelemsize"
maxelemsize = getdata{1};

%MIN ELEMENT SIZE ("minelemsize"):
%Its default value is "0", but the user can put any value.
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "maxelemsize"
minelemsize = getdata{1};

%-----------------------------------------------
%Type of ELEMENT (Triangle, Quadrangle or Both):

%"keyelem" verify if the element used in meshed domain will be a trinagle, 
%quadrangle or both. Is important attempt which in structured mesh we can 
%not use "both" option.
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "keyelem". 
keyelem = getdata{1};

%-----------------------------------------
%Mesh Algorithm and Refinement Parameters:

%"keyalgor" report the type of algorithm used to build the mesh. 
%You can choose the following options:
%[1] ==> "MeshAdapt"; [5] ==> "Delaunay"; [6] ==> "Frontal"
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "keyalgor"
keyalgor = getdata{1};

%"keyrefine" report the method used to refine the mesh already built. 
%You can choose [0] ==> "No split"; [1] ==> "Automatic
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "keyrefine"
keyrefine = getdata{1};

%Uniform REFINEMENT Parameter:
%"nref" obtain the number that define the quantity of refinement process to 
%use it onto refinement loop
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "nref"
nref = getdata{1};

%"scalemesh" apply a scale in mesh size (constant factor in all mesh).
%Naturaly, the default value is 1.
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "scalemesh"
scalemesh = getdata{1};

%--------------------------------
%Opening of gmsh's user interface

%The parameter "keyopen" is of string type which defines if the gmsh will 
%be opened to finish of mesh generation process.   
keyopen = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%--------------------------------------------------------------------------
%PHYSICAL PARAMETERS (Density):

%"keydens" is a key which define if the density value cames of *.dat file
%or will be constant
keydens = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%"densname" receives the file's name which contain the density map when 
%there is a strongly influence of pressure over the density.
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "densname"
densname = getdata{1};

%"ndens" is the number of fluids flowing in the domain. Each fluid has a 
%Density value (first "water", secondely "oil")
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "ndens"
ndens = getdata{1};

%Read the Density values:
getdata = textscan(readfile,'%f',ndens,'CommentStyle',{'//'});
%Attribute the data to "dens"
dens(1:ndens) = cell2mat(getdata);

%--------------------------------------------------------------------------
%PHYSICAL PARAMETERS (Viscosity):

%"keyvisc" is a key which define if the viscosity value cames of *.dat file
%or will be constant
keyvisc = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%"viscname" receives the file's name which contain the density map when 
%there is a strongly influence of pressure over the density.
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "viscname"
viscname = getdata{1};

%"nvisc" is the number of fluids flowing in the domain. Each fluid has a 
%Density value (first "water", secondely "oil")
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "nvisc"
nvisc = getdata{1};

%Read the Viscosity values:
getdata = textscan(readfile,'%f',nvisc,'CommentStyle',{'//'});
%Attribute the data to "visc"
visc(1:nvisc) = cell2mat(getdata);

%--------------------------------------------------------------------------
%PHYSICAL PARAMETERS (Residual Saturation):

%The first value read is referred to unredutible water saturation. The
%second one is referred to residual oil saturation.
getdata = textscan(readfile,'%f',2,'CommentStyle',{'//'});
%Attribute the data to "satlimit"
satlimit(1:2) = cell2mat(getdata);

%-----------------------------------
%PHYSICAL PARAMETERS (Permeability):

%There is three ways to obtain the permeability distribution:
%The "keyperm" difine it.
keyperm = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%Read the value of "klimits" (max and min). When the permeability is 
%randomic 
getdata = textscan(readfile,'%f',2,'CommentStyle',{'//'});
%Attribute the data to "klimits"
klimits(1:2) = cell2mat(getdata);

%"hipermname" receives the file's name which contain the permeability map 
%When a hightly heterogeneous media is choosen.
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "hipermname"
hipermname = getdata{1};

%"numlayer" read the amount of layers considered. 
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numlayer"
numlayer = getdata{1};

%"kmap" is vector who receives the tensor componentes in each layer        
getdata = textscan(readfile,'%f %f','CommentStyle',{'//'});
%Fill an auxiliary matrix
auxmat = cell2mat(getdata);
%Initialize "kmap" and other parameters:
kmap = zeros(numlayer,5);
c = 1;
%Swept all layers
for i = 1:3:3*numlayer
    %Attribute the tensor number to "kmap"
    kmap(c,1) = auxmat(i,1);
    %Attribute the tensor to a auxiliary vector
    auxvec(1:4) = auxmat(i + 1:i + 2,1:2)';
    %Finaly, "kmap" receives the components.
    kmap(c,2:5) = auxvec;
    %Increment the auxiliary counter "c"
    c = c + 1;
end  %End of IF

%--------------------------------------------------------------------------
%PHYSICAL PARAMETERS (Porosity):

%"keypor" defines if the porosity map cames of owner "Start.dat" or of a
%external file.
keypor = textscan(readfile,'%c',1,'CommentStyle',{'//'});

%"poramount" read the amount of diferent porous value which must attributed 
%to each material
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "poramount"
poramount = getdata{1};

%"pormap" is a vector which receives the porosity values
getdata = textscan(readfile,'%f',poramount,'CommentStyle',{'//'});
%Attribute the data to "pormap"
pormap = cell2mat(getdata);

%The porosity map must came of a external file. The name of file must read 
%from "Start.dat".
%"porname" receives the file's name which contain the permeab. map to 
%hightly heterogeneous media
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'},'whitespace','%');
%Attribute the data to "porname"
porname = getdata{1};

%--------------------------------------------------------------------------
%PHYSICAL PARAMETERS (Gravity and Capilarity):

%Verify if or no has gravity effects
keygravity = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%"g" is a vector which receives the gravity vector
getdata = textscan(readfile,'%f %f %f',poramount,'CommentStyle',{'//'});
%Attribute the data to "g"
g = cell2mat(getdata)';

%Verify if or no has capilarity effects
keycapil = textscan(readfile,'%c',1,'CommentStyle',{'//'});
%Get a exponent of the Corey Model
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "ncaplcorey"
ncaplcorey = getdata{1};

%--------------------------------------------------------------------------
%BOUNDARY CONDITIONS:

%"numbc" is the number of boundary conditions statment
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numbc"
numbc = getdata{1};

%"bcflag" reads the parameters related to boundary condition. 
%This is compared with *.geo file and the label of each edge may 
%change according the type of boundary condition seted.
getdata = textscan(readfile,'%f %f',numbc,'CommentStyle',{'//'});
%Attribute the data to "bcflag"
bcflag = cell2mat(getdata);

%"numbcnode" reads the amount of boundary condition applyed only in
%the nodes.
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numbcnode"
numbcnode = getdata{1};

%"bcnodepointer" reads the parameters related to boundary condition applyed 
%in each node. This will used to rename the node physical groups.
getdata = textscan(readfile,'%u',numbcnode,'CommentStyle',{'//'});
%Attribute the data to "bcnodepointer"
bcnodepointer = cell2mat(getdata);

%--------------------------------------------------------------------------
%CONVERGENCE CRITERIA:

%The Courant number will be used as a stability parameter of the hyperbolic 
%saturation equation. 
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "courant"
courant = getdata{1};

%"timedimkey" receives "s" to dimentional unit or "vpi" to adment. unit. 
getdata = textscan(readfile,'%s',1,'CommentStyle',{'//'});
%Attribute the data to "timedimkey"
timedimkey = getdata{1};

%Initialize "totaltime"
totaltime = zeros(2,1);
%When a dimentional unit is used
totaltime(1) = totaltime(1) + (strcmp(timedimkey,'pvi'));

%This parameter denots the time of simulation. Will be used in the unstead 
%state simulation. 
getdata = textscan(readfile,'%f',1,'CommentStyle',{'//'});
%Attribute the data to "totaltime"
totaltime(2) = getdata{1};

%--------------------------------------------------------------------------
%WELLS:

%"numwell" reads the number of wells seted 
getdata = textscan(readfile,'%u',1,'CommentStyle',{'//'});
%Attribute the data to "numwell"
numwell = getdata{1};

%Fill the matrix "well". This contains in its columns: the number
%of well (1, 2, ..., n); the cartesian position of well (x,y); the
%influence ray (ray of well hole); finaly, the information if it is 
%injector (0) or productor (1). 
if numwell > 0
    getdata = textscan(readfile,'%f64%f64%f64%f64%f64%f64%f64%f64%f64',...
        numwell,'CommentStyle',{'//'});
    %Attribute the data to "well"
    well = cell2mat(getdata);
%There is no any well
else
    well = 0;
end  %End of IF

%Close the data file ("Start.dat")
fclose(readfile);

%--------------------------------------------------------------------------
%FUNCTION "getcoord"
%--------------------------------------------------------------------------

function [coord,nnode] = getcoord(verifymshfile)

%Open the *.msh file
readmsh = fopen(verifymshfile);
%"nnode" is the number of nodes in the discrete domain
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',4);
%Attribute the data to "nnode"
nnode = getmshdata{1};

%"coord" is a matrix which contain the coordinate of each point of domain
getmshdata = textscan(readmsh,'%*u %f64 %f64 %f64',nnode,'HeaderLines',1);
%Fill the matrix "coord" with the x, y and z coordinates.
coord = cell2mat(getmshdata);

%Clear "getmshdata"
clear getmshdata;

%Close the *.msh file
fclose(readmsh);

%--------------------------------------------------------------------------
%FUNCTION "getelem"
%--------------------------------------------------------------------------

function [elem,nbe,nelem,nodelim,intnode,flaglim,inboundedge] = ...
    getelem(verifymshfile,nnode,numwell,well)
%Initialize "inboundedge"
inboundedge = 0;
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

%There is an INTERNAL LINE
%debug

if numwell > 0 && any(well(:,9) == 3) || numwell > 0 && any(well(:,9) == 0)
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n%n%*n%n%*n%n%n%*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Gives the contribution to "entitype"
    entitype = auxmat(:,1);

    %Get "inboundedge" from "auxmat"
    inboundedge = ...
        auxmat(logical(auxmat(:,1) ~= 15 & auxmat(:,2) > 1000),3:4);
    %Orient in encreasing order the vertices which define each boundary
    %edge (ex.: 5, 1 ==> 1, 5).
    inboundedge = sort(inboundedge,2);

    %"nbe" verifies how many entities are edges which constitute the 
    %boundary (scalar)
    nbeaux = sum(entitype == 1);
    nbe = nbeaux - size(inboundedge,1);
    %"nodelim". These nodes constitute the limits of geometry. 
    %Get the amount of internal node (scalar)
    intnode = sum(logical(auxmat(:,1) == 15 & auxmat(:,2) > 1000));
    %We exclude any node inside the domain.    
    nodelim = sum(entitype == 15) - intnode;
%There is no well
else
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Fill "entitype"
    entitype = auxmat(:,1);
    %"nbe" verifies how many entities are edges which constitute the 
    %boundary
    nbe = sum(entitype == 1);
    nbeaux = nbe;
    %"nodelimit". These nodes constitute the limits of geometry
    nodelim = sum(entitype == 15);
    intnode = 0;
end  %End of IF

%Define the amount of each entity:
%Verify how many "entities" are nodes. This value is attributed to 
%"ntri" verifies how many triangles there are in the domain
ntri = sum(entitype == 2); 
%"nquad" verifies how many quadrangles there are in the domain
nquad = sum(entitype == 3);
%"nelem" is the number of elements in the domain (triangles and/or 
%quadrangles)
nelem = ntri + nquad;

%Get "flaglim" from "auxmat". It is a vector with the flag associated to
%external nodes.
flaglim = auxmat(logical(auxmat(:,1) == 15 & auxmat(:,2) < 1000),2);

%Close the *.msh file
fclose(readmsh);

%---------------------------
%Define the matrix ("elem"):

%Open again the *.msh file
readmsh = fopen(verifymshfile);

%"elem" is a matrix which contain the nodes which constitute it
%Create and initialize the parameter "elem"
elem = zeros(nelem,5);  %nelem rows and 5 columns

%Fill the matrix "elem" with 4 nodes for a quadrangle element and 3 nodes 
%for a triangle element. In this case the fourth column receives null
%value. The fiveth column receive material propertie flag.
%There are tringles:
if ntri > 0  
    %Triangle
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u',ntri,...
        'HeaderLines',8 + nnode + nodelim + intnode + nbeaux);
    elemaux = cell2mat(getmshdata);
    %Attribute to "elem" (fifth column)
    elem(1:ntri,5) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(1:ntri,1:3) = elemaux(:,2:4); 

    %Close the *.msh file
    fclose(readmsh);
end  %End of first IF

%There are quadrangles:
if nquad > 0  
    %Open again the *.msh file
    readmsh = fopen(verifymshfile);
    
    %Quadrangles
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u %u',nquad,...
        'HeaderLines',8 + nnode + nodelim + intnode + nbeaux + ntri);
    elemaux = cell2mat(getmshdata);
    %Attribute to "elem" (fifth column)
    elem(ntri + 1:ntri + nquad,5) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(ntri + 1:ntri + nquad,1:4) = elemaux(:,2:5); 

    %Close the *.msh file
    fclose(readmsh);
end  %End of second IF

%Clear matrices
clear auxmat getmshdata elemaux; 

%It avoids that the fifth column of "elem" get null values
if any(logical(elem(:,5) == 0))
    %Fill the fifth column of "elem" with "1". Thus the material can be
    %attributed.
    elem(:,5) = 1;
end  %End of IF

%--------------------------------------------------------------------------
%FUNCTION "getinfoedge"
%--------------------------------------------------------------------------

function [bedge,inedge,klb] = getinfoedge(coord,elem,nnode,nbe,nelem,...
    nodelim,intnode,flaglim,verifymshfile,bcflag)

%"bedge" and "inedge" matrix - Contain respectively the node who constitute
%the boundary edges and iner edges, besides the element which contain such 
%nodes. In case of internal nodes "Inedge" offers the node to both left and 
%right side to edge evaluated 

%The number of external edge (where the boundary conditions are applieds) 
%may be obtained throughout *.msh generated by gmsh.
    
%"bedge" is a matrix which contain the two nodes which constitute each 
%external edge, the element which edge belongs and a geometric flag

%Create and initialize "bedgeaux" matrix. Its order is (nbe x 5). 
%"nbe" is the number of boundary edges. The first two columns are the nodes 
%which constitute the edge. The third node is the element which the edge 
%belongs (element to left). The last columns is the boundary condition's 
%type (Dirichlet or Neumann) flag. 
bedgeaux = zeros(nbe,5);

%Create and initialize "inedge". In principle the number of rows is not
%known. To long of fulfill of that we will know its row amount.
%However its column number is known. The two firsts values are the nodes
%who defines the edge. The two last values are respectively the left and
%right elements with regard to edge evaluated
inedge = zeros(1,4);

%Fill the matrix "bedge" with a boundary condition flag and the two nodes 
%which constitute it

%Open the *.msh file
readbound = fopen(verifymshfile);

%In principle "bedgeaux" read in its fourth column the geometry flag from 
%*.msh file generated by gmsh. After that "bedge" receives this info.
getboundata = textscan(readbound,'%*u %*u %*u %u %*u %u %u',nbe,...
    'HeaderLines',8 + nnode + nodelim + intnode);
%Attribute the data to "bedgeaux"
bedgeaux(:,5) = getboundata{1};
bedgeaux(:,1:2) = [getboundata{2} getboundata{3}];

%Close file
fclose(readbound);

%Fills "bedge" according to value in its fifth column:
%It points the correct row in "bedgeaux"
pntcorrectrow = ismember(bedgeaux(:,5),bcflag(:,1));
%Fill "bedge"
bedge = bedgeaux(pntcorrectrow,:);

klb = 0;
%Verify how many rows in "bedge" has known Neumann Boundary condition.
%Get the flag that feature a Neumann NON-Null boundary
pointflag = logical(bcflag(:,1) > 200 & bcflag(:,2) ~= 0);
%Define the "flagval"
flagval = bcflag(pointflag,1);
%For each non-null flag, we verify the row number and store this in "klb"
%(Known Line Boundary)
if any(pointflag)
    bedgnumrow = 1:size(bedge,1);
    %Swept the non-null bc.
    for i = 1:length(flagval)
        %Point to "bedge" row
        boundpoint = logical(bedge(:,5) == flagval(i));
        klb = union(klb,bedgnumrow(boundpoint));
    end  %End of FOR
    klb = klb(2:length(klb));
end  %End of IF
        
%Fill the fourth column, which receive the boundary code related to first 
%point of each element. The second point is always the first one of behind 
%row. Thus, the fourth column receives the values of the fiveth column.
bedge(:,4) = bedge(:,5);

%The nodes which was marked with a different bounday condition will be
%ajusted below.
%"nodelim" is the amount of nodes over the boundary
%"flaglim" stores the flags associated to the external nodes
for ipoint = 1:nodelim
    %"rowmod" are the rows which will be modified
    rowmod = logical(bedge(:,1) == ipoint);
    bedge(rowmod,4) = flaglim(ipoint);
end  %End of FOR

%Complete "bedge" with number of elements and fill "inedge".

%"inaux" is a counter of "inedge" rows because we do not know how many rows
%we will produce in this matrix
inaux = 1;

%Swept all elements of "elem" (rows)
for ielem = 1:nelem
    %Get the amount of edges in each element.
    amountedges = sum(elem(ielem,1:4) ~= 0);        
    %Swept all nodes of each element (columns)
    for inode = 1:amountedges
        %"next" is a variable which search the neighbour node in the
        %"elem" matrix. When inode is 4 (quadrangle mesh), "next" 
        %receives 1 because the neighbour of node 4 is the node 1. 
        %In any other situation the neighbour node is simply inode + 1 
        next = (inode == amountedges) + (inode + 1)*(inode < amountedges);
        %The same logical for "behaind"
        behaind = (inode == 1)*amountedges + (inode > 1)*(inode - 1);
            
        %------------------------------------------------------------------
        %>> When the edge belongs to boundary <<
            
        %Obtain the "bedge"'s row which contain "inode" and "next"
        i = all(elem(ielem,inode) == bedge(:,1:2) | elem(ielem,next) == ...
            bedge(:,1:2),2);

        %Verify if there is a common row between "findbedgerow1" 
        %and "findbedgerow2"
        if any(i)
            %If there is coincidence between nodes which are
            %conecteds to node ielem,inode and the nodes which
            %belong to boundary (obtained by *.msh file), the third
            %column of "bedge" receives the index which corresponds
            %to element in analisys.
            bedge(i,3) = ielem;
                
        %------------------------------------------------------------------
        %>> When the edge does not belongs to boundary <<
            
        %The couple of nodes which do not belong to "bedge" probably
        %belongs to "inedge". This ELSE attribute the couple of nodes to
        %"inedge". Before that, some conditions must be satisfy.
        else
            %Define the "inode" coordinate
            inodecoord = coord(elem(ielem,inode),:);
            %Get two vectors
            vecnext = coord(elem(ielem,next),:) - inodecoord;
            vecbehind =  coord(elem(ielem,behaind),:) - inodecoord;
            
            %Define the angle between the edges evaluated.
            angbtweenedges = acosd(dot(vecnext,vecbehind)/...
                (norm(vecnext)*norm(vecbehind)));
            %Get an angle condition. It is used because some edges form
            %angles near to 180 degrees. In order avoid ill interpretation
            %by "preprocessor", when the angle is near to 180, the way of
            %element is tuned counterclockwise.
            angcond = ((angbtweenedges/180) <= 0.99);

            %Get the cross product.
            pointsign = cross(vecnext,vecbehind);
            %Attribute just the cross "k" component.
            pointsign = pointsign(3)*angcond;
                            
            %Verify if couple of nodes already were contabilized. Case
            %yes, the element associate is a right element.
             
            %Obtain the "inedge"'s row which contain "inode" and "next"
            j = all(elem(ielem,inode) == inedge(:,1:2) | ...
                elem(ielem,next) == inedge(:,1:2),2);
            
            %Verify if the couple of nodes already were evaluated AND
            %if the node associated to "inode" is bigger than the node
            %associated to "next". In this case the element "ielem"
            %must be an element to right.
            if any(j) && elem(ielem,inode) > elem(ielem,next)
                %Define the column
                column = 4*(pointsign >= 0) + 3*(pointsign < 0); 
                %Attribute the element to column according to sign of
                %the "pointsign"
                inedge(j,column) = ielem;
            %The couple of nodes were evaluated but the node associated 
            %to "inode" is lower than the node associated to "next". In 
            %this case the element "ielem" is a left element.
            elseif any(j) && elem(ielem,inode) < elem(ielem,next)
                %Define the column
                column = 3*(pointsign >= 0) + 4*(pointsign < 0); 
                %Attribute the element to column according to sign of
                %the "pointsign"
                inedge(j,column) = ielem;
            %In case of the couple does not be already stored, Although
            %the node associated with "inode" is major than node
            %associated with "next". The "ielem" go to fourth column
            %too.
            elseif any(j) == 0 && elem(ielem,inode) > elem(ielem,next)      
                %Attribute the lower node to first column and the 
                %bigger node to second column
                inedge(inaux,1:2) = sort([elem(ielem,inode) ...
                    elem(ielem,next)]);
                %Define the column
                column = 4*(pointsign >= 0) + 3*(pointsign < 0); 
                %Attribute the element to column according to sign of
                %the "pointsign"
                inedge(inaux,column) = ielem;
                %Increment "inaux"
                inaux = inaux + 1;
            %The couple of nodes were never evaluated but the node
            %associated to "inode" is minor than node associated to 
            %"next". The element "ielem" go to third col. (left elem.)
            elseif  any(j) == 0 && elem(ielem,inode) < elem(ielem,next)
                %Attribute the lower node to first column and the 
                %bigger node to second column
                inedge(inaux,1:2) = sort([elem(ielem,inode) ...
                    elem(ielem,next)]);
                %Define the column
                column = 3*(pointsign >= 0) + 4*(pointsign < 0); 
                %Attribute the element to column according to sign of
                %the "pointsign"
                inedge(inaux,column) = ielem;
                %Increment "inaux"
                inaux = inaux + 1;
            end  %End of IF (fill "inedge")
        end  %End of IF ("bedge" and "inedge")
    end  %End of FOR (swept the vertices in each element)
end  %End of FOR (build "bedge" and "inedge")

%Clear matrices
clear getboundata bedgeaux;

%--------------------------------------------------------------------------
%Function "calcnormals"
%--------------------------------------------------------------------------

function [normals] = calcnormals(coord,centelem,bedge,inedge)
%Initialize "normals"
normals = zeros(size(bedge,1) + size(inedge,1),3);

%Initialize "bedgesize"
bedgesize = size(bedge,1);
%Initialize "iflowb"
iflowb = 1:bedgesize;
%Fill "normals" (in "bedge")
normals(iflowb,1:2) = ...
    fliplr((coord(bedge(iflowb,2),1:2) - coord(bedge(iflowb,1),1:2)));
normals(iflowb,2) = -normals(iflowb,2);
%Correct (if necessary) the normal orientation

%Calculate a pointer (centroid to out of control volume)
pointer = coord(bedge(iflowb,1),1:2) - centelem(bedge(iflowb,3),1:2);
%Correct (if necessary) the normal orientation
changesign = (pointer(:,1).*normals(iflowb,1) + ...
    (pointer(:,2).*normals(iflowb,2))) < 0;
%Correct the sign of "normals"
normals(changesign,1:2) = -normals(changesign,1:2);

%Initialize "iflowb"
iflowin = 1:size(inedge,1);
%Fill "normals" (in "inedge")
normals(bedgesize + iflowin,1:2) = ...
    fliplr((coord(inedge(iflowin,2),1:2) - coord(inedge(iflowin,1),1:2)));
normals(bedgesize + iflowin,2) = -normals(bedgesize + iflowin,2);

%Clear auxiliary vectors
clear iflowb iflowin pointer changesign;

%--------------------------------------------------------------------------
%Function "calcnormals_polarsys"
%--------------------------------------------------------------------------

function [normals] = calcnormals_polarsys(coord,centelem,bedge,inedge)
%Define "bedgesize"
bedgesize = size(bedge,1);
%Get the normals without the "r0" contribution
normals = calcnormals(coord,centelem,bedge,inedge);
%Attribute the intervention of polar coordinate
%Swept boundary edges
for iflowb = 1:size(bedge,1)
    %Get midedge coordinate
    midedgecoord = mean(coord(bedge(iflowb,1:2),:));
    normals(iflowb,:) = midedgecoord(1)*normals(iflowb,:);
end  %End of FOR

%Swept internal edges
for iflowin = 1:size(inedge,1)
    %Get midedge coordinate
    midedgecoord = mean(coord(inedge(iflowin,1:2),:));
    normals(iflowin + bedgesize,:) = ...
        midedgecoord(1)*normals(iflowin + bedgesize,:);
end  %End of FOR

%--------------------------------------------------------------------------
%FUNCTION "getsurnode"
%--------------------------------------------------------------------------

function [esurn,nsurn] = getsurnode(inode,esurn1,esurn2,nsurn1,nsurn2)
%"countelem" counts how many elements there are surrounding each node 
%evaluated.
countelem = esurn2(inode + 1) - esurn2(inode);
%"esurn" calculates what are the elements surrounding each node evaluated.
ie = 1:countelem;
esurn(ie) = esurn1(esurn2(inode) + ie);

%"countnode" counts how many nodes there are surrounding each node
%evaluated.
countnode = nsurn2(inode + 1) - nsurn2(inode);
%"nsurn" calculates what are the nodes surrounding each node 
%evaluated.
in = 1:countnode;
nsurn(in) = nsurn1(nsurn2(inode) + in);

%--------------------------------------------------------------------------
%FUNCTION "shiftchoosen"
%--------------------------------------------------------------------------

function [vecout] = shiftchoosen(vecin,numelem,letter)
%Initialize "vecout". It is the vector reordered according vector element
%choosen.
vecout(1:length(vecin),1) = vecin;

%It finds the position of element choosen in vector "vecin"
if strcmp(letter,'val')
    %Define "i" (auxiliary conuter 1:length(vecin))
    i = 1:length(vecin);
    %Poits its position
    pointpos = i(logical(vecin == numelem));
%"pointpos" receives the position
else
    pointpos = numelem;
end  %End of IF

%It uses "circshift" to reorder the vector
if pointpos > 1
    %From "pointpos" on
    auxvecout = vecout(pointpos:length(vecout));
    %From "pointpos" on
    auxvecout(length(auxvecout) + 1:length(vecout)) = ...
        vecout(1:pointpos - 1); 
    %Update "vecout"
    vecout = auxvecout;
end  %End of IF

%--------------------------------------------------------------------------
%FUNCTION "getcentelem"
%--------------------------------------------------------------------------

%This function calculate the coordinate of center's element evaluated 
%(xunkn, yunkn, zunkn). A vector "coordunkn" (3x1) is returned.
%"coordunkn(1)" does mean the coordinate x of unknown point.
%The inflow data is "elemnode" which is a vector with all nodes which
%constitute the element evaluated.
%That is the coordinate of each node who constitute the element evaluated.
%"elemnnode(4)" is the value of fourth column of "elem" matrix. If the elem 
%is a triangle "elemnode(4)" is 0 and the statment (elemnode(4) > 0) is 0. 
%If the elem is a quadrangle the statment is 1 and is added to 3. 
function [centelem] = getcentelem(coord,elem)                         
%Initialize the matrix "centelem"
centelem = zeros(size(elem,1),3);

%Poits to the triangle element 
pointri = elem(:,4) == 0;
%Poits to the quadrangle element 
pointquad = elem(:,4) ~= 0;
%Fill "centelem" (triangles)
centelem(pointri,1:2) = coord(elem(pointri,1),1:2);
centelem(pointri,1:2) = centelem(pointri,1:2) + ...
    coord(elem(pointri,2),1:2);
centelem(pointri,1:2) = (centelem(pointri,1:2) + ...
    coord(elem(pointri,3),1:2))./3;
%Fill "centelem" (quadrangles)
centelem(pointquad,1:2) = coord(elem(pointquad,1),1:2);
centelem(pointquad,1:2) = centelem(pointquad,1:2) + ...
    coord(elem(pointquad,2),1:2);
centelem(pointquad,1:2) = centelem(pointquad,1:2) + ...
    coord(elem(pointquad,3),1:2);
centelem(pointquad,1:2) = (centelem(pointquad,1:2) + ...
    coord(elem(pointquad,4),1:2))./4;

%Clear auxiliary vectors
clear pointri pointquad;

%--------------------------------------------------------------------------
%Function "calcelemarea"
%--------------------------------------------------------------------------

%This function calculate the area of each element. The element is recognaze
%due parameter "numelem"
function [elemarea] = calcelemarea(coord,elem)
%Initialize the vector
elemarea = zeros(size(elem,1),1);
%Swept all elements of domain
for ielem = 1:size(elem,1)
    %Calculate the element's area to quadrangular element
    if elem(ielem,4) > 0
        %"atriang1" calculates the area of first internal triangle
        atriang1 = 0.5*norm(cross((coord(elem(ielem,2),:) - ...
            coord(elem(ielem,1),:)),(coord(elem(ielem,2),:) - ...
            coord(elem(ielem,3),:))));
        %"atriang2" calculates the area of second internal triangle
        atriang2 = 0.5*norm(cross((coord(elem(ielem,4),:) - ...
            coord(elem(ielem,3),:)),(coord(elem(ielem,4),:) - ...
            coord(elem(ielem,1),:))));
        %The element's area is constituted by sum of two calc. areas above. 
        %This case worth whenever the element is a quadrangle 
        elemarea(ielem) = atriang1 + atriang2;
    %Calculate the element's area to triangular element
    elseif elem(ielem,4) == 0
        %"elemarea" is calculated in the one step
        elemarea(ielem) = ...
            0.5*norm(cross((coord(elem(ielem,2),:) - ...
            coord(elem(ielem,1),:)),(coord(elem(ielem,2),:) - ...
            coord(elem(ielem,3),:))));
    end  %End of IF
end  %End of FOR

%--------------------------------------------------------------------------
%Function "calcelemarea_polarsys"
%--------------------------------------------------------------------------

%This function calculate the area of each element. The element is recognaze
%due parameter "numelem"
function [elemarea] = calcelemarea_polarsys(coord,elem,centelemr)
%Get "elemarea"
elemarea = calcelemarea(coord,elem);
%Attribute the polar contribution ("r0")
elemarea = centelemr.*elemarea;

%--------------------------------------------------------------------------
%Function "getesure"
%--------------------------------------------------------------------------

function [esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,...
    bedgesize,inedge,esurn1,esurn2,nsurn2,rowposit)
%Initialize the pointers "esureface2" and "esurefull2"
esureface2 = 0;
esurefull2 = 0;
%Initialize counters
mface = 0;
mfull = 0;
%Sewpt all control volumes
for ielem = 1:size(elem,1)
    %Initialize "neighborelemface" and "neighborelemfull"
    neighborelemface = 0;
    neighborelemfull = 0;
    %Get the amount of edge in each control volume
    amountedge = sum(logical(elem(ielem,1:4) ~= 0));
    %Initialize the auxiliary counter
    m = 1;
    cface = 1;
    cfull = 1;
    boundflagface = 0;
    boundflagfull = 0;
    %Define "elemvertices", "vertex" and "next"
    elemvertices = elem(ielem,1:amountedge);
    %Swept the internal half-edges
    for j = 1:amountedge + 1
        %In each face, define the first and the second ("next") vertices.
        vertex = elemvertices(1);
        next = elemvertices(2);
        %Get the row number for the edge evaluated
        %Get the position for each vertex:
        survtx1 = rowposit(nsurn2(vertex) + 1:nsurn2(vertex + 1));
        survtx2 = rowposit(nsurn2(next) + 1:nsurn2(next + 1));
        %Verify which position matches between the two vertices:
        %Initialize "nsurnptr"
        nsurnptr = zeros(length(survtx1),1);
        %Swept "nsurn"
        for isur = 1:length(survtx1)
            pointinters = logical(survtx1(isur) == survtx2);
            nsurnptr(isur) = sum(pointinters);
        end  %End of FOR
        %Get the position
        pointrow = survtx1(logical(nsurnptr));
        
        %Verify if "pointrow" belongs to "inedge"
        if pointrow > bedgesize
            %Get the neighbor elements:
            %It stores initialy neighbor of faces but put the neighbor of
            possiblelem = 1:2;
            pointactualelem = logical(inedge(pointrow - bedgesize,3:4) ~= ...
                ielem);
            %vertices (if "m" == 2, see below)
            neighborelemfull(cfull) = inedge(pointrow - bedgesize,2 + ...
                possiblelem(pointactualelem));
            %It stores only neighbor of face.
            neighborelemface(cface) = inedge(pointrow - bedgesize,2 + ...
                possiblelem(pointactualelem));
            %when "m" is equal to 2, verify if exists vertex neighboring
            if m == 2
                %Initialize "neighboraux"
                neighboraux = 0;
                %Get the elements surrounding the first vertex of second
                %edge evaluated.
                esurn = esurn1(esurn2(vertex) + 1:esurn2(vertex + 1));
                %Reorder the elements position in "esurn"
                [newesurn] = shiftchoosen(esurn,...
                    neighborelemfull(length(neighborelemfull) - 1),'val');
                %Evaluate if there is vertex neighbor
                elemvector = [ielem ...
                    neighborelemfull(length(neighborelemfull) - 1:...
                    length(neighborelemfull))];

                %Verify "setdiff"
                %Initialize "nsurnptr"
                nsurnptr = zeros(length(newesurn),1);
                %Swept "nsurn"
                for idiff = 1:length(newesurn)
                    pointinters = logical(newesurn(idiff) == elemvector);
                    nsurnptr(idiff) = sum(pointinters);
                end  %End of FOR
                %Get the position
                extraelem = newesurn(logical(nsurnptr == 0));

                %There are extra elements
                if any(extraelem)
                    neighboraux(1) = ...
                        neighborelemfull(length(neighborelemfull) - 1);
                    neighboraux(2:length(extraelem) + 1) = extraelem;
                    neighboraux(length(extraelem) + 2) = ...
                        neighborelemfull(length(neighborelemfull));
                    %Attribute to "neighborelem" the new neighbor
                    neighborelemfull(cfull - 1:cfull + ...
                        length(neighboraux) - 2) = neighboraux;
                    %Redifine "m" and "c"
                    m = 1;
                    cfull = cfull + length(extraelem);
                end  %End of IF
            end  %End of IF
            %Update "m" and "c"
            m = m + 1;
            cface = cface + 1;
            cfull = cfull + 1;
        %The edge belongs to "bedge"
        else
            %Decrement "m" in order it does not reach 2.
            m = 1;
            boundflagface = cface;
            boundflagfull = cfull;
        end  %End of IF
        
        %It shifts the vertices position in each control volume:
        %The control volume is a Triangle
        if amountedge == 3
            elemvertices = circshift(elemvertices,1:3);
        %The control volume is a Quadrangle
        else
            elemvertices = circshift(elemvertices,2:4);
        end  %End of IF
    end  %End of internal FOR (each control volume)
    %Ensure there is no repeated element
    neighborelemface = unique(neighborelemface,'stable');
    neighborelemfull = unique(neighborelemfull,'stable');

    %Reorder adequadely "neighborelemface"
    if boundflagface ~= 0 && boundflagface <= length(neighborelemface)
        neighborelemface = shiftchoosen(neighborelemface,...
            neighborelemface(boundflagface),'val');
    end  %End of IF ("boundflagface")
    %Reorder adequadely "neighborelemfull"
    if boundflagfull ~= 0 && boundflagfull <= length(neighborelemfull)
        neighborelemfull = shiftchoosen(neighborelemfull,...
            neighborelemfull(boundflagfull),'val');
    end  %End of IF ("boundflagface")
        
    %Attribute the contribution from "esureface" and "esurefull" to vectors 
    %"esureface1" and "esurefull1"
    esureface1(mface + 1:mface + length(neighborelemface)) = ...
        neighborelemface;
    esurefull1(mfull + 1:mfull + length(neighborelemfull)) = ...
        neighborelemfull;
    
    %Attribute the contribution from "esureface" and "esurefull" to vectors 
    %"esureface2" and "esurefull2"
    esureface2 = vertcat(esureface2,length(esureface1));
    esurefull2 = vertcat(esurefull2,length(esurefull1));

    %Increment counters "mface" and "mfull"
    mface = length(esureface1);
    mfull = length(esurefull1);
end  %End of FOR (Swept all elements)

%Turn "esureface1" and "esurefull1" vertical vector
esureface1 = esureface1';
esurefull1 = esurefull1';

%--------------------------------------------------------------------------
%FUNCTION "getrowposition"
%--------------------------------------------------------------------------

function [rowposit] = getrowposition(bedge,inedge,nsurn1,nsurn2)
%Initialize "rowposit"
rowposit = zeros(length(nsurn1),1);
%Initialize "bedgesize"
bedgesize = size(bedge,1);

%Swept the boundary edges
for i = 1:bedgesize
    %Get the vertices:
    vertices = bedge(i,1:2);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(1)) + 1:nsurn2(vertices(1) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    nsurnpos = nsurnsize(nsurn == vertices(2));
    %Fill "rowposit"
    rowposit(nsurn2(vertices(1)) + nsurnpos) = i;

    %----------------------
    %Evaluate the vertex 2:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(2)) + 1:nsurn2(vertices(2) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    nsurnpos = nsurnsize(nsurn == vertices(1));
    %Fill "rowposit"
    rowposit(nsurn2(vertices(2)) + nsurnpos) = i;
end  %End of FOR ("bedge")

%Swept the internal edges
for i = 1:size(inedge,1)
    %Get the vertices:
    vertices = inedge(i,1:2);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(1)) + 1:nsurn2(vertices(1) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    nsurnpos = nsurnsize(nsurn == vertices(2));
    %Fill "rowposit"
    rowposit(nsurn2(vertices(1)) + nsurnpos) = i + bedgesize;

    %----------------------
    %Evaluate the vertex 2:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(2)) + 1:nsurn2(vertices(2) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    nsurnpos = nsurnsize(nsurn == vertices(1));
    %Fill "rowposit"
    rowposit(nsurn2(vertices(2)) + nsurnpos) = i + bedgesize;
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%FUNCTION "reordernsurn"
%--------------------------------------------------------------------------

function [nsurn1] = reordernsurn(esurn1,esurn2,nsurn1,nsurn2,bedge,coord,...
    elem)
%Put "nsurn1" in a vector column
nsurn1 = nsurn1'; 

%--------------------------------------------------------------------------
%Fernando's modification

% a aloca��o do vetor serve para identificar todos os nos
% do contorno e damos flag os demais que est�o no interior � zero
is_bound = zeros(size(coord,1),1);
for i = 1:size(bedge,1),
   is_bound(bedge(i,1)) = 1;
   % procura se tem algum no se repitidno no coluna 1 do bedge
   % quando passa isso que dizer que ordena��o do bedge esta 
   % sentido horario e/ou antihorario
   auxc = find(bedge(:,1) == bedge(i,1));
   if ~isempty(auxc) && length(auxc) > 1
       is_bound(bedge(auxc(1),1)) = 1;
   end
   % procura se o n� esta se repitindo na coluna 2 bedge
  
    auxc1 = find(bedge(:,2) == bedge(i,2));
    if ~isempty(auxc1) && length(auxc1) > 1
        is_bound(bedge(auxc1(1),2))=1;  
    end
end  %End of FOR

%--------------------------------------------------------------------------

for i=1:size(coord,1),
    if is_bound(i) == 0
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if elem(prim,4)==0
            bp=3;
        else
            bp=4;
        end
        if elem(ult,4)==0
            bu=3;
        else
            bu=4;
        end
        for j=1:bp,
            for u=1:bu,
                if ((elem(ult,u)==elem(prim,j))&&(elem(prim,j)~=i))
                    nnsn(1)=elem(prim,j);
                end
            end
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if elem(atual,4)==0
                bp=3;
            else
                bp=4;
            end
            if elem(ant,4)==0
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end     
    else 
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if elem(prim,4)==0
            bp=3;
        else
            bp=4;
        end
        if elem(ult,4)==0
            bu=3;
        else
            bu=4;
        end
        for j=1:size(bedge,1)
            for p=1:bp,
                if (i==elem(prim,p))
                    if p==bp
                        n1=elem(prim,1);
                    else
                        n1=elem(prim,p+1);
                    end
                end
            end
            nnsn(1)=n1;            
            for u=1:bu,
                if (i==elem(ult,u))
                    if u==1
                        n1=elem(ult,bu);
                    else
                        n1=elem(ult,u-1);
                    end
                end
            end
            nnsn((nsurn2(i+1)-nsurn2(i)))=n1;
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if elem(atual,4)==0
                bp=3;
            else
                bp=4;
            end
            if elem(ant,4)==0
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end
    end            
end






