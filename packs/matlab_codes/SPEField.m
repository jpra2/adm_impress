

load('spe_perm.dat')
spe_perm = reshape(spe_perm',1,[]);
%60 x 220 x 85 cells (1.122x106 cells).
startPoint = (layer-1)*60*220 +1;
endPoint =  startPoint + 60*220 ;
% 
kxx = spe_perm(startPoint:endPoint-1)';
clear spe_perm

if triangular == 1
    flagS = logical(mod(1:size(elem,1),2)); 
    elem(flagS,end) = 1:(size(elem,1)/2);
    elem(~flagS,end) = 1:(size(elem,1)/2);
else
    elem(:,end) = 1:size(elem,1);

end

kym = zeros(size(kxx));
kym = reshape(kxx, 60,220);

%kym = rot90(kym,2);
kym = flip(kym,1);
if Rotacionado == 1
    kym = kym';    
end

if triangular == 1
    kxx = reshape(kym,1,[]); 
    kmap(1:(size(elem,1)/2) ,1) = 1:(size(elem,1)/2); 
    kmap(1:(size(elem,1)/2),2) = kxx(1:(meshX*meshY));
    kmap(1:(size(elem,1)/2),5) =  kxx(1:(meshX*meshY));
    
    
else
    kxx = reshape(kym,1,[]); 
    kmap(1:size(elem,1),1) = 1:size(elem,1); 
    kmap(1:size(elem,1),2) = kxx(1:(meshX*meshY));
    kmap(1:size(elem,1),5) =  kxx(1:(meshX*meshY));

end

%postprocessorName(kmap(elem(:,end),2),elem(:,end),superFolder, 'PermFIELDSPE');