ui = [ 20 50 100 150 200 250 280 300 ];
Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18]; 

ud= -600:1:600;
Md = []
M_spline3 = []
M_Poly3 = []
M_Poly5 = []


h= 0.001;    % so high cause required sampling for high frequency signal
t= 0: h: 30; 
f = 50;   
% e = 240*sin(t) ;
e = 120*sin(2*pi*f*t);

y= zeros( 3, length(t) ); 

%%%%%% RK4 method
    for i=1:length(t)-1 
     k1 = F(t(i), y(:,i),f);
     k2 = F(t(i)+h/2, y(:,i) + h/2*k1,f); 
     k3 = F(t(i)+h/2, y(:,i) + h/2*k2,f);
     k4 = F(t(i)+h, y(:,i) + h*k3,f);
      y(:, i+1) = y(:,i) + h/6*(k1+2*k2+2*k3+k4);
    end 

 figure(22)
    %hold on
    plot(t,Md,t,M_sp,t,M_P3,t,M_P5); grid %state variable in 3rd row of matrix   
    title("e and capacitor voltages for three methods"); xlabel('t[s]'); ylabel('U[V]');  
    legend('Lagrange','e')
   
 figure(33)
    %hold on
    plot(t,Md,t,M_sp,t,M_P3,t,M_P5);grid  %second row
    title("i_1 [A] for e(t) = 120sin(2πft), for f = 50Hz "); xlabel('t[s]'); ylabel('i[A]');  
    legend('i_1')


function dy = F(t,y,f)
% dy=[ f1(t, y(1), y(2), y(3))]

R1=0.1; 
R2=10; 
C=0.5; 
L1=3; 
L2=5; 
% M=0.8;

% e=100*sin(t);
% f = 0.68;
% f=50;
% e =  240*sin(t);

% h= 0.01;  
% t= 0: h: 30;
e = 210*sin(2*pi*f*t);
% T = 3;    
% e = 120 * (rem(t, T) < T/2);

% di1dt = 1/(L1/M-M/L2)*(-R1/M*y(1)+R2/L2*y(2)-1/M*y(3)+1/M*e)
% di2dt = 1/(M/L1-L2/M)*(-R1/L1*y(1) + R2/M*y(2) - 1/L1*y(3)+1/L1*e)
% uL1 = L1*di1dt + M*di2dt;

i1= y(1);
uC= y(3);

uL1 = e - R1*i1 - uC;
% e = sin(t);
M = fM(uL1);
% M=0.8;

dy = [ 1/(L1/M-M/L2)*(-R1/M*y(1)+R2/L2*y(2)-1/M*y(3)+1/M*e)
       1/(M/L1-L2/M)*(-R1/L1*y(1) + R2/M*y(2) - 1/L1*y(3)+1/L1*e)
       1/C * y(1) ];
end



for i = 1:length(ud)
    Md(i) = fM_lagr(t,y(1, :),ud(i));
end

for i = 1:length(ud)
    M_sp(i) = fM_spline3(t,y(1, :),ud(i));
end

for i = 1:length(ud)
    M_P3(i) = fM_Poly3(t,y(1, :),ud(i));
end

for i = 1:length(ud)
    M_P5(i) = fM_Poly5(t,y(1, :),ud(i));
end

% for i = 1:length(ud)
%     s(i) = fM_spline(ud(i));
% end

% function M_s = fM_spline(uL1)   %default matlab spline function for verification purposes
% ui = [ 20 50 100 150 200 250 280 300 ];
% Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];
% 
% uL1 = abs(uL1);
%     if uL1>300
%         M_s = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
%                    % for voltage greater than 300V mutual inductance is invalid
%                    % [typical problem of extrapolation (determining values 
%                    % of function outside the given range of measurements)]
%     elseif uL1 < 20
%         M_s = 0.46;
%     else 
% 
%         M_s = spline(ui, Mi, uL1); %default matlab spline to verify
% 
%     end
% end


% figure(66)
% plot(ui,Mi, 'o-', ud, Md, ud, s);grid;
% title("Lagrange interpolation on mutual inductance"); xlabel('uL1[V]'); ylabel('M[H]');  
% legend('Table data','Lagrange interpolation','s')
% %xline([4])
function M_lagr = fM(uL1)   
ui = [ 20 50 100 150 200 250 280 300 ];
Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M_lagr = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    else 
          M_lagr = mylagr(ui,Mi,uL1);
    end
end

function M_lagr = fM_lagr(ui,Mi,uL1)   
% ui = [ 20 50 100 150 200 250 280 300 ];
% Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M_lagr = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    else 
          M_lagr = mylagr(ui,Mi,uL1);
    end
end

function M_spline3 = fM_spline3(uL1)   
% ui = [ 20 50 100 150 200 250 280 300 ];
% Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M_spline3 = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    elseif uL1 < 20
        M_spline3 = 0.46; %limiting the spline from extrapolating unreasonable and invalid values
    else 
    
        M_spline3 = spline3(ui, Mi, uL1);
    
    end
end

function M_Poly3 = fM_Poly3(uL1)   
% ui = [ 20 50 100 150 200 250 280 300 ];
% Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M_Poly3 = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    else 
        M_Poly3 = Poly3(ui,Mi,uL1);
    end
end

function M_Poly5 = fM_Poly5(uL1)   
% ui = [ 20 50 100 150 200 250 280 300 ];
% Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M_Poly5 = 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    else 
        M_Poly5 = Poly5(ui,Mi,uL1);
    end
end

figure(77);
plot(ui, Mi, 'o', ud, Md,'-', ud, M_sp, ud, M_P3, ':', ud, M_P5, ':');
title('Interpolation and Polynomial Approximations');
xlabel('uL1[V]');
ylabel('M[H]');
legend('Data Points','Lagrange interpolation', 'Spline Interpolation', 'Degree 3 Polynomial', 'Degree 5 Polynomial','matlab spline');
grid on;

function y = mylagr(xi,yi,x) %goal: be able to implement this
% Lagrange polynomial interpolation from memory 
y=0;
n=length(xi);
for i = 1:n
    L=1;
    for j=1:n
        if i ~= j
            L= L*(xi(j) - x)/(xi(j) - xi(i));
        end
    end
    y= y+yi(i)*L;
end
end

function yy = spline3(x,y,xx)
yy = zeros(length(xx),1);
n = length(x)-1;
A = zeros(n*4);
b = zeros(n*4,1);
    for r = 1:n
    i = (r-1)*4 + 1;
    A(i,i) = 1;
    A(i,i+1) = x(r);
    A(i,i+2) = x(r)^2;
    A(i,i+3) = x(r)^3;
    b(i) = y(r);
    A(i+1,i) = 1;
    A(i+1,i+1) = x(r+1);
    A(i+1,i+2) = x(r+1)^2;
    A(i+1,i+3) = x(r+1)^3;
    b(i+1) = y(r+1);
        if (r ~= n)
        A(i+2,i+1) = 1;
        A(i+2,i+2) = 2*x(r+1);
        A(i+2,i+3) = 3*x(r+1)^2;
        A(i+2,i+1+4) = -1;
        A(i+2,i+2+4) = -2*x(r+1);
        A(i+2,i+3+4) = -3*x(r+1)^2;
        A(i+3,i+2) = 2;
        A(i+3,i+3) = 6*x(r+1);
        A(i+3,i+2+4) = -2;
        A(i+3,i+3+4) = -6*x(r+1);
        end
    end
A(end,end-2) = 1;
A(end,end-1) = 2*x(end);
A(end,end) = 3*x(end)^2;
b(end) = -0.001;                 %decreasing these values fixed overfitting (Runge's phenomenon)
A(end-1,2) = 1;
A(end-1,3) = 2*x(1);
A(end-1,4) = 3*x(1)^2;
b(end-1) = 0.001;                %decreasing these values fixed overfitting (Runge's phenomenon)
a = A\b
% A
    for k = 1:length(xx)
    i = findidx(x, xx(k));
    r = (i-1)*4+1;
    yy(k) = a(r) + a(r+1) * xx(k) + a(r+2) * xx(k)^2 + a(r+3) * xx(k)^3;
    
 
    end
    
end

function i = findidx(x, xx)
    for i = 1:length(x)
        if x(i) > xx
            i = i - 1;
        break;
        end
    end

    if i == length(x)
    i = i - 1;
    end
end

function y = Poly3(xi, yi, x)

    n = length(xi);     
    A = ones(n, 4);     % vandermonde matrix
        for i = 1:3
            A(:, i+1) = xi'.^i;
        end
   
 coefs = (A' * A) \ (A' * yi');  % solving the normal equations
        
 y = zeros(size(x)); % evaluating the polynomial at the given x
        for i = 0:3
            y = y + coefs(i+1) * x.^i;
        end
end

function y = Poly5(xi, yi, x)    

    n = length(xi);
    A = ones(n, 6);  
    for i = 1:5
        A(:, i+1) = xi'.^i;
    end

    coefs = (A' * A) \ (A' * yi');    
 
    y = zeros(size(x));     
    for i = 0:5
        y = y + coefs(i+1) * x.^i;
    end
end


