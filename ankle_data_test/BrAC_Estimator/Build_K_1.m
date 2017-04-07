function K_1 = Build_K_1(n)
      
D = [n*[1,2*ones(1,n-1),1]];
   
S = [-n*[ones(1,n)]];
   
K_1 = diag(D) + diag(S,1) + diag(S,-1);

end
