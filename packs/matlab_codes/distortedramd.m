function [auxcoord]=distortedramd
global coord bedge

auxcoord=zeros(size(coord,1),4);

a=-1;
b=1;
% para malha quadilateral da SPE usei limites 10 e -10
% com alpha 0.8 e h=1/220
ex = a + (b-a)*rand(size(coord,1),1);
ey = a + (b-a)*rand(size(coord,1),1);
% quando mas proximo a 1 eh mas ditorcido
alpha=0.35;
h=1000/25;

for icoord=1:size(coord,1)
    x=coord(icoord,1);
    
    y=coord(icoord,2);
    if icoord>size(bedge,1) %&& abs(coord(icoord,1)-0.5)>1e-10
      % if single(y)~=0.5       
%        if x~=0.5 && y~=0.5
            x=coord(icoord,1)+h*alpha*ex(icoord,1);
            
            y=coord(icoord,2)+h*alpha*ey(icoord,1);
            
            auxcoord(icoord,1)=icoord;
            auxcoord(icoord,2)=x;
            auxcoord(icoord,3)=y;
%          else
%              auxcoord(icoord,1)=icoord;
%              auxcoord(icoord,2)=x;
%              auxcoord(icoord,3)=y;
%         end
    else

        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;
        
    end
end

end
