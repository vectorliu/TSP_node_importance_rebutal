%Using an simple example to see the divergence
%Topology is just the simple example used in the TSP paper
%According to our results, the results of Bertrand seems skeptical

clc;clear all;close all;

%%%Topology
nodes = 10; 
delta = 0.1;
one = ones(nodes,1);

adja = [0 1 0 0 1 0 0 1 1 1;1 0 0 0 0 0 0 0 1 0;0 0 0 1 0 1 1 1 1 0;...
    0 0 1 0 0 1 0 0 0 0;1 0 0 0 0 0 1 1 1 1;0 0 1 1 0 0 1 1 0 0;...
    0 0 1 0 1 1 0 1 0 1;1 0 1 0 1 1 1 0 1 1;1 1 1 0 1 0 0 1 0 1;1 0 0 0 1 0 1 1 1 0];
laplacian = diag(sum(adja,2)) - adja;
matrixC = eye(nodes) - delta*laplacian;

%%%Initial vector
initial = zeros(nodes,1);
initial(1) = 1;
initial = laplacian*initial;

%%%True Fiedler vector
[v,d] = eig(laplacian);
Fiedler = v(:,2);
All = ones(nodes,1)/sqrt(nodes);

iter = 1e3;
FiedlerE = zeros(nodes, iter);
Error = zeros(1, iter);
FiedlerE(:,1) = initial;
FiedlerE(:,1) = FiedlerE(:,1)/norm(FiedlerE(:,1), 2);
Error(1) = (norm(FiedlerE(:,1) - Fiedler, 2))^2/nodes;
ErrorOne = zeros(1,iter);
ErrorOne(1) = (norm(FiedlerE(:,1) - All, 2))^2/nodes;

for i=1:iter-1
    temp = matrixC*FiedlerE(:,i);
    FiedlerE(:,i+1) = temp/norm(temp, 2);
    
    Error(i+1) = (norm(FiedlerE(:,i+1) - Fiedler, 2))^2/nodes;
    ErrorOne(i+1) = (norm(FiedlerE(:,i+1) - All, 2))^2/nodes;
end

figure;
plot(1:iter, Error(:), 'r-.', 1:iter, ErrorOne(:), 'b--', 'linewidth', 2);
xlabel('iteration');
ylabel('mean-squared estimation error');
