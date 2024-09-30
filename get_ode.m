%depending on the task, use this function to call ODE system. If "param" is not
%provided, we work with baseline paramaters , otherwise 
% "param" represents a single parameter set coming from virtual patients.
function dydt = get_ode(t,y,param)
    dydt = zeros(4,1);
    if nargin==2
        r1 = 0.18;
        r2 = 0.21;
        K1 = 5.1e+3;
        l = 1.56;
        dC = 0.41;  
        muC = 0.6;
        s = 3.05e-1;
        gammaC = 2.93e-2;
        k = 2.019e-7;
        wC = 3e-5;
        K2 = 1.65e3;
        dB = 0.3;
        muB = 0.89;
        b = 1.4e-3;
        gammaB = 7e-3;
        wB = 3.42e-6;
    elseif nargin==3
        %TsTr params
        r1 = param(1);
        r2 = param(2);
        K1 = param(3);
        %C params
        l = param(4);
        dC = param(5);           
        muC = param(6);                                                                                                                                                                                                                                                                  
        s = param(7);
        gammaC = param(8);
        k = param(9);
        wC = param(10);
        K2 = param(11);
        %B params
        dB = param(12);
        muB = param(13);
        b = param(14);
        gammaB = param(15);
        wB = param(16);
    end
    Ts = y(1); Tr = y(2); C = y(3); B = y(4); 
    DC = dC*(C/Ts)^l/(s+(C/Ts)^l);
    DS = dB*(B/Ts)^l/(s+(B/Ts)^l);
    DR = dB*(B/Tr)^l/(s+(B/Tr)^l);

    dydt(1) = r1*Ts*(1-(Ts+Tr)/K1) - DC*Ts - DS*Ts; %Ts
    dydt(2) = r2*Tr*(1-(Ts+Tr)/K1)- DR*Tr; %Tr
    dydt(3) = - gammaC*C -muC*log((B + C)/K2)*(DC^2/(k+DC^2))*C - wC*C*Ts;%C
    dydt(4) = b - gammaB*B - muB*log((C+B)/K2)*DS^2/(k+DS^2)*C-wB*(Tr+Ts);%B

end
