function [laplacian, fiedler, delta, lambda] = create_ER_graph(nodes, p)
%%%function to create the random ER graph
%%%input:number of nodes, connection probability
%%%output:laplacian matrix, fiedler vector, and optimal delta = 1/\lambda_n
epsilon = 1e-5;
while 1
    adja = zeros(nodes, nodes);
    for i=1:1:nodes
        for j=i+1:1:nodes
            if rand(1,1) <= p
                adja(i,j) = 1;
                adja(j,i) = 1;
            end
        end
    end
    laplacian = diag(sum(adja,2)) - adja;
    [v,d] = eig(laplacian);
    fiedler = v(:,2);
    %delta = (0.5)./max(sum(adja,2));
    delta = 1./d(end,end);
    lambda = d(2,2);
    if abs(d(2,2)) > epsilon%To judge whether connected network
        break;
    end
end
end