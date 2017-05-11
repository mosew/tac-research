function kernel = computeKernel(s)
    load('operators.mat','ANhat','BNhat','tau');
    kernel = mpower(ANhat,fix(s/tau)-1) * BNhat;
end