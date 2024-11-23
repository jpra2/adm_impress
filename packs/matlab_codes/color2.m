clear;
clc;
% This code solves the node coloring problem.
% The data is taken from the coursera online course on discrete
% optimization.https://www.coursera.org/learn/discrete-optimization#about
%The python cpdess I developed for that course are reproduced
% in MATLAB.The data is in python index format.[counts from 0]
% Part 1 data formatting.
%-------------------
fname='gc_20_1';
%fname='gc_50_3';

data=load(fname);
%The data matris is  (m+1)X2 matrix;
% m the number of edges(connecting elements).
%The first row is the number of  nodes(vertices,points) n and m.
n=data(1,1);m=data(1,2);
data(1,:)=[];
data=data+1;% python index starts with zero;
%-------------------------------

%Part 2 Integet programming model formulation

k=round(n/2);% The number of maximum colors(Groups)
k=7;
N=n*k+k;% The number of binary variables in this formulation.

f=[zeros(n*k,1);ones(k,1)];% The linear objective vector.
%The objective is to minimize the number of colors used.%
%i.e min sigma(y_k)
%2a Equality constraint sigma(xi_k)=1;
r1=n*(0:k-1);

Aeq=[repmat(eye(n),1,k) zeros(n,k)];
beq=ones(n,1);
% % Part 2b Inequality constraints  x_i_k+x_j_k=<1.
% 
for i=1:m
Am1(i,data(i,:))=1;
end
 Ay1=zeros(k*m,k);
 for i=1:k
Ay1(1+m*(i-1):i*m,i)=1;
end

 Aieq=[kron(eye(k) ,Am1)  -Ay1];
bieq=zeros(k*m,1);
 

intcon=1:N;
lb=zeros(N,1);
ub=lb+1;
%intlinprog.options.MaxTime = 72000 ;
 [x,F] = intlinprog(f,intcon,Aieq,bieq,Aeq,beq,lb,ub);
 x1=x(1:n*k);
 y=x(k*n+1:end);
x1=reshape(x1,n,k);
 x1=x1(:,y>.5);
 %X:The colors of the nodes.
 X=x1*(1:F)'
 %F: The number of colors used.
 F
% 
% %% eq
% 
% % 
% % 
