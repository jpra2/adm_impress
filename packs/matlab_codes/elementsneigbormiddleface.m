function [element1]=elementsneigbormiddleface
global coord esurn1 esurn2 centelem
j=1;
auxnec1=0;
for i=1:size(coord,1)
    x=coord(i,1);
    y=coord(i,2);
    if (x==0.5 && y==0.5)
    %if ((0.25<x || x==0.25) &  (x<0.75 || x==0.75)) && y==0.5
        nec1=esurn2(i+1)-esurn2(i);
        p=1;
        for k=1:nec1
            
            if 0.25<centelem(esurn1(esurn2(i)+k),1) && centelem(esurn1(esurn2(i)+k),1)<0.75
                element1(p+auxnec1,1)=esurn1(esurn2(i)+k);
                p=p+1;
            end
        end
        auxnec1=auxnec1+p-1;
        j=j+1;
    end
end
k=1;
n=length(element1);
while k<=n
    j=1;
    while j<=n
        if k~=j
            if element1(k)==element1(j)
                element1(j)=[];
                n=length(element1);
            end
        end
        j=j+1;
    end
    k=k+1;
end

end