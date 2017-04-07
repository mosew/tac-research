function Q = Build_Q_M(n)
      
D = [(1/n)*[(1/3),(2/3)*ones(1,n-1),(1/3)]];
   
S = [(1/n)*[(1/6)*ones(1,n)]];
   
Q = diag(D) + diag(S,1) + diag(S,-1);

end

