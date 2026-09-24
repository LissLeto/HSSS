function para=Getpara2(LAB,max_l)
[L,A,B]=imsplit(LAB);
varL=var(L(:));varA=var(A(:));varB=var(B(:));
Lrate=2*varL/(varA+varB);
alpha=min(Lrate,max_l);
ABrate=max(varA,varB)/min(varA,varB);
rest=(3-alpha)/(log(ABrate)+1);
if varB>varA
    gemma=rest;
    beta=3-alpha-gemma;
else
    beta=rest;
    gemma=3-alpha-beta;
end
para=[alpha,beta,gemma];
%para=[2,0.5,0.5];