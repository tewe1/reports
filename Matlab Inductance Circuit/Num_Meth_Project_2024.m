global n_eval
n_eval = 0

% this function is just for purpose of testing:
ui = [ 20 50 100 150 200 250 280 300 ];
Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18]; 

ud= -600:1:600;
Md = []
for i = 1:length(ud)
    Md(i) = fM(ud(i));
end

figure(1)
plot(ui,Mi, 'o-', ud, Md);grid;
title("Lagrange interpolation of mutual inductance"); xlabel('uL1[V]'); ylabel('M[H]');  
legend('Table data','Lagrange interpolation')
%xline([4])


% % Part 4 (determine freq for 406W)

% Bisection method
% a,b,c: frequencies f of the input signal we are searching for !!!!!!!!
a = 0.6;
b = 1;      % a and b are some reasonable guesses for the boundaries
P_a = power_simulator (a)
FF_a = 406 - P_a            % we're looking for a frequency that will 
                            % give us 406W on the resistors

P_b = power_simulator (b)
FF_b = 406 - P_b

for i=1:8  % for 20 it's slow (take for 10)
    c = (a+b)/2;
    P_c = power_simulator(c);
    FF_c = 406 - P_c;

    if FF_a * FF_c < 0 
        b = c;  %if there is sign change, take c as the value of b
        FF_b = FF_c;
    else 
        a = c;
        FF_a = FF_c;
    end
end
 c

% % Secant method
% % a, b: initial guesses for the frequency f
% a = 0.6;
% b = 1;
% 
% P_a = power_simulator(a);    % Calculate function values at initial guesses
% FF_a = 406 - P_a;
% P_b = power_simulator(b);
% FF_b = 406 - P_b;
% 
% for i = 1:6       
%     c = b - FF_b * (b - a) / (FF_b - FF_a);
% 
%     P_c = power_simulator(c);
%     FF_c = 406 - P_c;
% 
%     a = b;
%     FF_a = FF_b;
%     b = c;
%     FF_b = FF_c;
% 
%     % if abs(FF_c) < 1e-6  % tolerance (if < some small value)
%     %     break;
%     % end
% end

% % Quasi-Newton Method 
% f0 = 0.7;             % initial guess foe freq
% delta_c = 1e-6;      % small increment for numerical derivative
% tolerance = 1e-2;    % tol for convergence
% max_iter = 25;      % max number of iterations
% 
% c = f0;
% for i = 1:max_iter
%     Pi = power_simulator(c);
%     Fi = Pi - 406;
% 
%     if abs(Fi) < tolerance
%         break;  % convergence criterion met
%     end
% 
%     P_prime = power_simulator(c + delta_c);
%     F_prime = P_prime - 406;
% 
%     dF_df = (F_prime - Fi) / delta_c;  % numerical derivative
%                                        % ( the delta(f) from project instructions )  
%     c = c - Fi / dF_df; % updating the frequency
% end
% c
% i


function Power = power_simulator(f)
%%%%%% here our simulator starts(finds ODE solutions using Euler's method)
%
h= 0.01;    % so high cause required sampling for high frequency signal
t= 0: h: 30; 
f = 50;    %  if used, not using the power simulator (we have given freq 5/50)
% e = 240*sin(t) ;
e = 120*sin(2*pi*f*t);
%
y= zeros( 3, length(t) ); % circuit off-state initial values (zeros)
yh= zeros( 3, length(t) );
yRK= zeros( 3, length(t) ); % generate empty matrix with 3 rows and same columns as
% y= [ 0;
%      0;
%      0 ]; %we fill in initial values(redundant)    
    for i = 1: length(t) - 1
           % y(1, i+1) = y(1, i) + h * f1(t(i),y(1,i),y(2,i),y(3,i));
           % y(2, i+1) = y(2, i) + h * f2(t(i),y(1,i),y(2,i),y(3,i));
           % y(3, i+1) = y(3, i) + h * f3(t(i),y(1,i),y(2,i),y(3,i));        
        y(:, i+1) = y(:, i) + h * F(t(i), y(:, i), f); %Euler method
    end   
 % figure(2)
 %    %hold on
 %    plot(t,y(3, :));grid %state variable in 3rd row of matrix
 %    title("u_C [V] for Euler's method")
 % figure(3)
 %    %hold on
 %    plot(t, y(1, :), t, y(2, :));grid  %second row
 %    title("i_1, i_2 [A] for Euler's method")
   
    %%%%%%%%%%%%%%%%%%%%%
    %Heun's method
    for i= 1: length(t) - 1    
        p= yh(:,i) + h* F(t(i), yh(:,i), f);        
        yh(:, i+1) = yh(:,i) + (h/2) *( F(t(i),yh(:,i),f) + F(t(i+1), p, f) );        
    end    

    %%%%%%%%%%%%%%%%%%%%%
    %RK4 method
    for i=1:length(t)-1 
     k1 = F(t(i), yRK(:,i),f);
     k2 = F(t(i)+h/2, yRK(:,i) + h/2*k1,f); 
     k3 = F(t(i)+h/2, yRK(:,i) + h/2*k2,f);
     k4 = F(t(i)+h, yRK(:,i) + h*k3,f);
      yRK(:, i+1) = yRK(:,i) + h/6*(k1+2*k2+2*k3+k4);
    end 

    % % T = 3;    
    % % e = 120 * (rem(t, T) < T/2);

 figure(2)
    %hold on
    plot(t,y(3, :),t,yh(3, :), t,yRK(3, :),t,e); grid %state variable in 3rd row of matrix   
    title("e and capacitor voltages for three methods"); xlabel('t[s]'); ylabel('U[V]');  
    legend('Euler','Heun','RK4','e')
   
 figure(3)
    %hold on
    plot(t, y(1, :), t, y(2, :),t,yh(1, :), t, yh(2, :),t,yRK(1, :), t, yRK(2, :));grid  %second row
    title("i_1, i_2 [A] for e(t) = 120sin(2πft), for f = 50Hz "); xlabel('t[s]'); ylabel('i[A]');  
    legend('i_1 Euler','i_2 Euler','i_1 Heun','i_2 Heun','i_1 RK4','i_2 RK4')

 % Calculating the power (heat dissipated on the resistors)

    i1_t = y(1,:); %samples in time (vector of the values(heights of the rectangles in integration))
    i2_t = y(2,:);
    R1=0.1; 
    R2=10;   
    Power=0;
    
% composite rectangles rule
    for i= 1:length(t)-1
            Power = Power + (R1*i1_t(i)^2 + R2*i2_t(i)^2) * h;    
    end

% Trapezoidal rule
% for i = 1:length(t)-1
%     Power = Power + ((R1*i1_t(i)^2 + R2*i2_t(i)^2) + (R1*i1_t(i+1)^2 + R2*i2_t(i+1)^2)) * h / 2; 
% end
    

    for i= 1:length(t)-1
        pow(i+1) = R1*i1_t(i)^2 + R2*i2_t(i)^2 ;
    end
pow
    figure(99)
    plot(t,pow); grid
    
Power

   % Power_simple = sum(R1*i1_t.^2 + R2*i2_t.^2)*h %same but simpler implementation
end


% Part 2
% Objective: Analyze how different interpolation / approximation method 
% change the solution of the simulation.
%
% coupling M between the primary and secondary circuits is non-linear 
% and depends on the voltage uL1 across the inductance L1 in the primary circuit, as shown in this table. 
% values in the table from the manufacturer sheet or from our laboratory experiments. 
% key concept here is: 
% We need to determine the value of the mutual inductance for each step of the time simulation.

% To avoid ‘extrapolation’ problem for the negative voltages, 
% the program should find the absolute value of the voltage. 
% The formula for mutual inductance should be modified then to MN(|uL1|).

%%% determine the traces of currents and voltages in the analyzed system for the given excitations
%%% and time periods the same as in the first part but only for the RK4 algorithm!

% mutual inductance is non-linear, so the variable M is replaced with a function

function M = fM(uL1)   
ui = [ 20 50 100 150 200 250 280 300 ];
Mi = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];

uL1 = abs(uL1);
    if uL1>300
        M= 0.18;   % fixing the values of M for uL1 > 300 to 0.18 [H]
                   % for voltage greater than 300V mutual inductance is invalid
                   % [typical problem of extrapolation (determining values 
                   % of function outside the given range of measurements)]
    else 
        M = mylagr(ui,Mi,uL1);
    end
end

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


function dy = F(t,y,f)
% dy=[ f1(t, y(1), y(2), y(3))]
    global n_eval;
    n_eval=n_eval+1;

R1=0.1; 
R2=10; 
C=0.5; 
L1=3; 
L2=5; 
% M=0.8;

% e=100*sin(t);
% f = 0.68;
 f=5;
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
  %M = fM(uL1);
 M=0.8;

dy = [ 1/(L1/M-M/L2)*(-R1/M*y(1)+R2/L2*y(2)-1/M*y(3)+1/M*e)
       1/(M/L1-L2/M)*(-R1/L1*y(1) + R2/M*y(2) - 1/L1*y(3)+1/L1*e)
       1/C * y(1) ];
end

n_eval




%function dy = f1(t, y1, y2, y3)
%R1=0.1; R2=10; C=0.5; L1=3; L2=5; M=0.8; e=sin(t);
%dy=1/(L1/M-M/L2)*(-R1/M*y1+R2/L2*y2-1/M*y3+1/M*e);
%end

%function dy = f2(t, y1, y2, y3)
%R1=0.1; R2=10; C=0.5; L1=3; L2=5; M=0.8; e=sin(t);
%dy=1/(M/L1-L2/M)*(-R1/L1*y1 + R2/M*y2 - 1/L1*y3+1/L1*e);
%end

%function dy = f3(t, y1, y2, y3)
%R1=0.1; R2=10; C=0.5; L1=3; L2=5; M=0.8; e=sin(t);
%dy= 1/C * y1;
%end

% y(1,2)=y(1,1)+h*f1(t0,y(1,1),y(2,1),y(3,1))
% y(2,2)=y(2,1)+h*f2(t0,y1,y2,y3)
% y(3,2)=y(3,1)+h*f3(t0,y1,y2,y3)

%!!!!!!! Please provide waveform diagrams for the following
% excitations
% on your own for the report!!!!!!!!!!!!!!!!:
% 1 [V], T = 3 [s]
% e(t) = 240sin(t),
% e(t) = 210sin(2πft), for f = 5 Hz
% e(t) = 120sin(2πft), for f = 50Hz
% e(t) = 120 if t < T/2 else 0

        % krotkie wprowadzenie o czym to jest
        % krotko o metodach, ich cechy charakkterystyczne
        % ZERO kodu w raporcie (raport merytoryczne, zalaczenie plikow zrodlowych jako archiwum zip)
        % na jednym rysunku trzy wykresy
        % znalezc takie h przy ktorym beda zauwazalne roznice (dla malych krokow h roznice sa nieznaczace, 
        % np nie oplaca sie uzywac 4 stopnia Runge Kutty)

