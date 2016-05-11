%Use random ER topology to calculate the Fiedler vector, and to evaluate
%the performance of our 5 methods: ME, ZE, SE, DZE, DSE

%Manually control the precision
%Perfromance evaluated by NRD1, i.e., the distorted order with respect to
%the most important node
%To sum different Monte-Carlo methods, we fix the network scale to be n=100

function Performance_precision()
clc;clear all;close all;
%%%Basic parameters definition
k = 500;%Times of Monte-Carlo simulations
n = 100;%Fixed network size, 100
%p = 0.5;%Connect probability of ER graph
precision = 8:-1:1;

%%%Generate the random graph
matrixL = cell(k,1);
for i=1:k
    p = 0.3 + 0.5*rand(1);%0.3~0.8
    [matrixL{i}.laplacian, matrixL{i}.fiedler, matrixL{i}.delta, matrixL{i}.lambda] = create_ER_graph(n, p);
end
%%%Main results save
MEC = zeros(size(precision, 2), k);
ZEC = zeros(size(precision, 2), k);
SEC = zeros(size(precision, 2), k);
DZEC = zeros(size(precision, 2), k);
DSEC = zeros(size(precision, 2), k);

tic
%%%Main loop to evaluate performance
for i=1:size(precision,2)
    pre = precision(i);
    precision(i) = 10^(-pre);
    for j=1:k
        fiedler = matrixL{j}.fiedler;
        fiedler = normalization(fiedler);
        fiedler = roundn(fiedler, -pre);
        delta = matrixL{j}.delta;
        L = matrixL{j}.laplacian;
        lambda = matrixL{j}.lambda;
        
        %%%calculate the theoretical vector and their order
        [theoretical, ME, ZE, SE, DZE, DSE] = calculate_approximated_descent(L, fiedler, delta, lambda);
        
        %%%calculate the NRD1 criteria
        MEC(i, j) = compare_order(ME, theoretical, 1);
        ZEC(i, j) = compare_order(ZE, theoretical, 1);
        SEC(i, j) = compare_order(SE, theoretical, 1);
        DZEC(i, j) = compare_order(DZE, theoretical, 1);
        DSEC(i, j) = compare_order(DSE, theoretical, 1);
    end
end
toc

%%%Plot the results to see the impact of precision
figure;%Mean of 1000 iterations
semilogx(precision, mean(MEC, 2), 'r*-', precision, mean(ZEC, 2), 'gs-.', precision, mean(SEC, 2), 'bo-', ...
    precision, mean(DZEC, 2), 'gd-.', precision, mean(DSEC, 2), 'bp-.', 'linewidth', 1.5, 'markersize', 8);
legend('ME(DME)','ZE','SE','DZE','DSE');
xlabel('Precision of Fiedler vector');
xlim([10^(-9),10^(0)]);
ylabel('NRD1 of different methods');
ylim([-1,16]);

figure;%Plot DSE, with median and quantile
quan = quantile(DSEC(:,1:k), [.25 .5 .75], 2);
semilogx(precision, quan(:,1), 'r*-', precision, quan(:,2), 'gs-.', precision, quan(:,3), 'bo--', 'linewidth', 1.5, 'markersize', 8);
legend('25% quantile', 'median', '75% quantile');
xlabel('Precision of Fiedler vector');
xlim([10^(-9),10^(0)]);
ylabel('NRD1 of DSE method');
ylim([-0.5, 2.5]);

save('performance_precision.mat', 'ME', 'ZE', 'SE', 'DZE', 'DSE', 'matrixL', 'MEC', 'ZEC', 'SEC','DZEC', 'DSEC', 'theoretical', 'k', 'precision');

end

function Delta = calculate_Delta(L, i)
%%%Calculate the descent of \Delta_i
temp = L;
temp(i,:) = zeros(1, size(L,1));
temp(:,i) = zeros(size(L,1), 1);
temp = temp - diag(sum(temp,1));
Delta = L -temp;
end

function [theoretical, ME, ZE, SE, DZE, DSE] = calculate_approximated_descent(L, fiedler, delta, lambda)
%%%calculate the approximated descent
theoretical = zeros(size(fiedler, 1), 1);
ME = zeros(size(fiedler, 1), 1);ZE = zeros(size(fiedler, 1), 1);
SE = zeros(size(fiedler, 1), 1);DZE = zeros(size(fiedler, 1), 1);
DSE = zeros(size(fiedler, 1), 1);
identity = eye(size(fiedler, 1));
one = ones(size(fiedler, 1), 1);

for i=1:size(fiedler, 1)
    %calculate matrix \Delta_i first
    Delta = calculate_Delta(L,i);
    %calculate the theoretical descent
    eigen = eig(L - Delta);
    algebra = eigen(3);
    theoretical(i) = -algebra;
    %calculate criteria
    ME(i) = fiedler'*Delta*fiedler;
    ZE(i) = fiedler'*Delta*(fiedler - fiedler(i)*identity(:,i))./(fiedler'*(fiedler - fiedler(i)*identity(:,i)));
    SE(i) = fiedler'*Delta*(identity - delta*L + delta*Delta - ...
        one*one'./size(fiedler, 1) - identity(:,i)* identity(:,i)')*(fiedler - fiedler(i)*identity(:,i))./(fiedler'*(identity - delta*L + delta*Delta...
        - one*one'./size(fiedler, 1) - identity(:,i)* identity(:,i)')*(fiedler - fiedler(i)*identity(:,i)));
    DZE(i) = fiedler'*Delta*(fiedler - fiedler(i)*identity(:,i));
    DSE(i) = fiedler'*Delta*(identity - delta*L + delta*Delta -...
        one*one'./size(fiedler, 1) - identity(:,i)* identity(:,i)')*(fiedler - fiedler(i)*identity(:,i))./(1 - delta.*lambda);
end
end

function disorder = compare_order(approximated, base, number)
%%%calculate the distortion with respect to the most critical node
disorder = 0;
[~, indices1] = sort(approximated, 'descend');%The largerer, the more important
[~, indices2] = sort(base, 'descend');%The larger, the more important
for i=1:number
    for j=1:length(approximated)
        if (indices1(j) == indices2(i)) && (j-i)>0
            disorder = disorder + j - i;
        end
    end
end
end

function vector = normalization(base)
%%%normalization
vector = base/norm(base, 2);
end
