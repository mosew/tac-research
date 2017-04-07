function Q = Build_Q_K(n)
      
D = [n*[1,2*ones(1,n-1),1]];
   
S = [-n*[ones(1,n)]];
   
Q = diag(D) + diag(S,1) + diag(S,-1);

end

