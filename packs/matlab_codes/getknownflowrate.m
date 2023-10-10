function [pressure,flowrate,flowresult] = getknownflowrate
%Define Global parameters:
global elem bedge inedge normals;

%Initialize "sidelength"
sidelength = 0.25;
%Initialize "v"
v = [1 0];
%Initialize "bedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "pressure", "flowrate" and "flowresult"
pressure = ones(size(elem,1),1);
flowresult = zeros(size(elem,1),1);
flowrate = zeros(bedgesize + inedgesize,1);

%Swept "bedge"
for i = 1:bedgesize
    %Get "leftelem"
    leftelem = bedge(i,3);
    %Get the flowrate
    flowrate(i) = dot(v,normals(i,1:2))/sidelength;
    %Attribute "flowrate" to "flowresult"
    flowresult(leftelem) = flowresult(leftelem) + flowrate(i);  
end  %End of FOR

%Swept "inedge"
for i = 1:inedgesize
    %Get "leftelem" and "rightelem"
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    %Get the flowrate
    flowrate(bedgesize + i) = dot(v,normals(bedgesize + i,1:2))/sidelength;

    %Attribute "flowrate" to "flowresult"
    flowresult(leftelem) = flowresult(leftelem) + flowrate(bedgesize + i);  
    flowresult(rightelem) = flowresult(rightelem) - flowrate(bedgesize + i);  
end  %End of FOR

 