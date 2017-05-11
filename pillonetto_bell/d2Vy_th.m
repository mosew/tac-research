function d2Vy_th = d2Vy_th(theta,P,T,n,eivs,Lmatrix,dLmatrix,d2Lmatrix)

    % OUTPUT:
    % n x n x nTheta x nTheta

    eivs = evaluate_eivs(eivs,theta);
    
    nTheta=length(theta);
    
    deivs = deiv(theta,P,T);
%     d2Lmatrix_ = d2Lmatrix(theta,P,n,Lmatrix);
    d2eiv_ = d2eiv(theta,P,T);
    
    d2Vy_th = zeros(n,n,nTheta,nTheta);
    
    for i = 1:n
        for k = 1:n
            for s = 1:nTheta
                for r = 1:nTheta
                    d2Vy_th(i,k,s,r) = sum( d2eiv_(:,s,r).*Lmatrix(:,i).*Lmatrix(:,k) + ...
                                            deivs(:,r).*(dLmatrix(:,i,s).*Lmatrix(:,k)+Lmatrix(:,i).*dLmatrix(:,k,s)) + ...
                                            deivs(:,s).*(dLmatrix(:,i,r).*Lmatrix(:,k)+Lmatrix(:,i).*dLmatrix(:,k,r)) + ...
                                            eivs.*(d2Lmatrix(:,i,s,r).*Lmatrix(:,k)+Lmatrix(:,i).*d2Lmatrix(:,k,s,r) + ...
                                                    dLmatrix(:,i,s).*dLmatrix(:,k,r)+dLmatrix(:,i,r).*dLmatrix(:,k,s)) );
                end
            end
        end
    end
    

    d2Vy_th = d2Vy_th + d2Vy_thu(theta,P,n);


    
end