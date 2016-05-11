%Use random ER topology to calculate the Fiedler vector, to evaluate the divergence
%Created by Hao Liu, May 10, 2016.
%Generate the initial vector randomly, and then apply x(0)=Lx(0)
%Two random factors: random graph (10~20 nodes), random initial vector

%%%Main function
function Fiedler_divergence_random()
clc;clear all;close all;
%%%Basic parameters definition
k = 1000;%Times of Monte-Carlo simulations
n = 9 + randi(11, k, 1);%Uniformly ranodm network size, 10~20
p = 0.4;%Connect probability of ER graph
iter = 1000;%iterations in each monte-carlo simulation

%%%Generate the random graph
matrixL = cell(k,1);
for i=1:k
    [matrixL{i}.laplacian, matrixL{i}.fiedler, matrixL{i}.delta] = create_ER_graph(n(i), p);
end

%%%Main loop to calculate the mean-squared error
error_to_fiedler = zeros(k, iter);
error_to_one = zeros(k,iter);
for i=1:k
    %%%True values setting
    L = matrixL{i}.laplacian;
    fiedler = matrixL{i}.fiedler;
    delta = matrixL{i}.delta;
    All = ones(n(i),1)/sqrt(n(i));
    %%%Initial vector
    initial = randn(n(i),1);
    initial = L*initial;
    initial = normalization(initial);
    
    error_to_fiedler(i,1) = mean_squared_error(fiedler, initial);
    error_to_one(i,1) = mean_squared_error(All, initial);
    %%%matrix C
    %C = eye(n(i)) - one*one'/n(i) - delta*L;
    C = eye(n(i)) - delta*L;
    
    for j=1:iter-1
        temp = C*initial;
        initial = temp/norm(temp, 2);
        %%%To control the precision manually
        %initial = roundn(initial, -1);
    
        error_to_fiedler(i, j+1) = mean_squared_error(fiedler, initial);
        error_to_one(i, j+1) = mean_squared_error(All, initial);
    end    
end
%%%plot the error
figure;
quan = quantile(error_to_fiedler(1:k,:),[.25 .5 .75],1);
semilogy(1:iter, quan(1,:), 'r-', 1:iter, quan(2,:), 'b--', 1:iter, quan(3,:),'g-.', 1:iter, error_to_fiedler(1,:),'k-', 'linewidth', 2);
xlabel('iteration');
ylabel('mean-squared estimation error');
legend('25% quantile', 'Median', '75% quantile', 'result of single run');
end

function vector = normalization(base)
%%%normalization
vector = base/norm(base, 2);
end

function error = mean_squared_error(base, guess)
%%%calculate meansquare error
error1 = (norm(base - guess, 2))^2/size(base, 1);
error2 = (norm(base + guess, 2))^2/size(base, 1);
error = min(error1, error2);
end