function M_1 = Build_M_1(n)
      
D = [(1/(3*n))*[1,2*ones(1,n-1),1]];
   
S = [(1/(6*n))*[ones(1,n)]];
   
M_1 = diag(D) + diag(S,1) + diag(S,-1);

end
