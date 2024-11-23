%função que calcula o somatório das vazões em cada elemento

function q =calflowratetotal(influx)
% somatório das vazões em cada elemento
global elem inedge bedge coord
q=sparse(size(elem,1),1);
bflux=influx(1:size(bedge,1));
% faces da fronteira
for i=1:size(bedge,1)
    lef=bedge(i,3);
    q(lef)=q(lef)+bflux(i);
end
% faces internas
for i=1:size(inedge,1)
    lef=inedge(i,3);
    rel=inedge(i,4);
    q(lef)=q(lef)+influx(i+size(bedge,1));
    q(rel)=q(rel)-influx(i+size(bedge,1));
end
end
