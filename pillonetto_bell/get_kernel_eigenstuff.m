function [eivs,eifs]=get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile)
    eivs = zeros(1,P); %numbers
    eifs = cell(1,P); %functions of t
    for j=1:P
        [eivs(j),eifs{j}] = feval(rkhs_eigenfile,j,theta,T);
    end
end