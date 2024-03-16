%function that calculates the flux density on boundary edges - Eq (32)

function VN_B = bflux(bedge,P,BI,BJ)

VN_B = zeros(size(bedge,1),1);

for i = 1:size(bedge,1)
    Iele = bedge(i,3);
    VN_B(i) = BI(i) * P(Iele) + BJ(i);
end