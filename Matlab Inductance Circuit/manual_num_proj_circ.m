R1=0.1; 
R2=10; 
C=0.5; 
L1=3; 
L2=5; 
M=0.8;

ZRC=R1+(-1i*C);

e=1;    % in case of this kind of calculation, e=1 means e=1*sin(t)

I = [-ZRC-(1i*L1)  1i*M
    -1i*M          R2+(1i*L2)] \ [-e;0]

UL1 = e*(1i*L1)/(ZRC+1i*L1)

Uc = e-(R1*(I(1))-((-1i*C)*I(1))-UL1)

abs(Uc)


abs(I(1))
abs(I(2))

ratioI1toI2 = abs(I(1))/abs(I(2))
