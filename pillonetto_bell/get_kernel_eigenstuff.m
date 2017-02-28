function [eivs,eifs]=get_kernel_eigenstuff(P,T,string)
    eivs = cell(1,P); %functions of q
    eifs = cell(1,P); %functions of t
    for j=1:P
        [eivs{j},eifs{j}] = feval(string,j,T);
    end
end