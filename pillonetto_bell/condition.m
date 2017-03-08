function bool = condition(epsilon,T,lumF,lumC,i)


    fiFsq = @(t) feval(lumF(i),t).^2;
    fiFminusfiCsq = @(t) (feval(lumF(i),t)-feval(lumC(i),t)).^2;

    bool = epsilon^2 * integrate(fiFsq,0,T) >= integrate(fiFminusfiCsq,0,T);
    
end