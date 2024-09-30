clearvars;clc;close all;
%Use this code to gradually optimize ODE sets. Based on odetype, decide 
% what to optimize, bounds and init guess. Modify the bounds (lb,ub),
% init guesses init_val. By default, we optimize r1 and K in Ts, 
% r2 in TsTr and l,dC,jC in TsTrC and dB,jB in TsTrCB equations. 
% To add more variables to optimize, modify get_odetype accordingly.

%time-parameters
time_max = 30;dt = 0.1;tspan = 3:dt:time_max;
cart_dose = 0.1;

%datatype: control,treatment,bystander
%odetype: ["Ts","TsTr","TsTrC","TsTrCB"].
odetype = "TsTrCB";
sp = .9; datatype="bystander";
[exp_days,data] = get_data(sp,datatype);

%%r1,K1-->Ts
% lb = [0,1e3];
% ub = [0.18,6e3];
% init_val = [0.1,1e3];
% [sol,sumsq] = fit_param(lb,ub,init_val,tspan,exp_days,cart_dose,data,sp,odetype);
% fprintf("r1: %0.3f K1: %0.1d\n",sol.s(1),sol.s(2))
% fprintf('mae: %0.2f\n',sumsq);


%%r2-->TsTr
% lb = [0.10];
% ub = [0.21];
% init_val = [0.1];
% [sol,sumsq] = fit_param(lb,ub,init_val,tspan,exp_days,cart_dose,data,sp,odetype);
% fprintf("r2: %0.3f\n",sol.s(1))
%fprintf('mae: %0.2f\n',sumsq);


%l,dC,muC-->TsTrC
% lb = [1.55,0.41,0.60];
% ub = [1.57,0.41,0.62];
% init_val = [1.5,0.4,0.1];
% [sol,sumsq] = fit_param(lb,ub,init_val,tspan,exp_days,cart_dose,data,sp,odetype);
% fprintf("l: %0.3f dC: %0.3f muC: %0.3f ",sol.s(1),sol.s(2),sol.s(3))
%fprintf('mae: %0.2f\n',sumsq);


%optimize dB,muB-->TsTrCB
% lb = [0.30,0.89];
% ub = [0.30,0.89];
% init_val = [0.28,0.88];
% [sol,sumsq] = fit_param(lb,ub,init_val,tspan,exp_days,cart_dose,data,sp,odetype);
% fprintf("dB: %0.3f muB: %0.3f ",sol.s(1),sol.s(2))
%fprintf('mae: %0.2f\n',sumsq);




%plot the results
[tumor_total,Ts,Tr,C,B] = solve_ode(sol.s,tspan,exp_days,cart_dose,sp,odetype);
relative_errors = abs((tumor_total(1:end) - data(1:end))./data(1:end));
mre = mean(relative_errors);
fprintf('\nmre: %0.2f\n',mre);
%plot(exp_days,tumor_total,'r-')
plot(tspan,Ts+Tr,'r-')

hold on;
plot(exp_days,data,'bo','MarkerSize',10)
legend({'model','data'},Location='best')
title(sprintf('model: %s --> sp: %0.0f%%',odetype,100*sp))
xticks(1:2:time_max);
xlim([0,time_max]);
%ylim([0,1000])
grid on;

% figure();
% plot(tspan,Ts,'b--')
% hold on
% plot(tspan,Tr,'r-')
% xticks(1:2:time_max);
% xlim([0,time_max]);
% grid on;
% legend('sensitive','resistant',Location='northeast')



function [optim_sol,sumsq] = fit_param(lb,ub,init_val,tspan,exp_days,cart_dose,data,sp,odetype)
    s = optimvar('s',length(ub),"LowerBound",lb,"UpperBound",ub);
    s0.s = init_val;
    [tumor_total,~,~] = fcn2optimexpr(@solve_ode,s,tspan,exp_days,cart_dose,sp,odetype);
    obj = sum(sqrt((tumor_total(2:end) - data(2:end)).^2))/length(data(2:end));
    prob = optimproblem("Objective",obj);
    options = optimoptions(prob);
    options.ConstraintTolerance = 1e-15;
    options.Display = 'iter';
    [optim_sol,sumsq] = solve(prob,s0,'Options',options);
end



function [tumor_total,Ts,Tr,C,B] = solve_ode(s,tspan,exp_days,cart_dose,sp,odetype)
    %we return Ts,Tr,C,B all but some are just dummies in many case.
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [~,indices] = ismember(exp_days,tspan);
    if odetype=="Ts"  
        [~,y] = ode45(@(t,y) get_odetype(t,y,s,odetype),tspan,y0,opts);
        Ts = y(:,1);
        Tr = 0;
        C = 0;
        B = 0;
        tumor_total = Ts(indices);
    
    elseif odetype=="TsTr"
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        [~,y] = ode45(@(t,y) get_odetype(t,y,s,odetype),tspan,[Ts0,Tr0],opts);
        Ts = y(:,1);
        Tr = y(:,2);
        C = 0;
        B = 0;
        tumor_total = Ts(indices)+Tr(indices);

    elseif odetype=="TsTrC"
        index4 = find(tspan==4);
        index6 = find(tspan==6);
        %tumor grows without intervention until day 7
        t1 = tspan(1:index4);
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        C0 = 0.;
       
        [~,y] = ode45(@(t,y)get_odetype(t,y,s,odetype),t1,[Ts0,Tr0,C0]);
        Ts = y(:,1);
        Tr = y(:,2);
        C = y(:,3);
        
        % day = 7 first injection
        t2 = tspan(index4+1:index6);
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = cart_dose;
        
        [~,y] = ode15s(@(t,y)get_odetype(t,y,s,odetype),t2,[Ts0,Tr0,C0],opts);
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        
        % day = 9, second injection
        t3 = tspan(index6+1:end);
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = C(end) + cart_dose;
     
        [~,y] = ode15s(@(t,y) get_odetype(t,y,s,odetype),t3,[Ts0,Tr0,C0],opts);
       
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = 0;
        
        %disp([length(t1),length(t2),length(t3)])
        tumor_total = Ts(indices)+Tr(indices);
   
    elseif odetype=="TsTrCB"
        y0 = 50;
        opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
        [~,indices] = ismember(exp_days,tspan);
        index4 = find(tspan==4);
        index6 = find(tspan==6);
        %tumor grows without intervention until day 3
        t1 = tspan(1:index4);
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        C0 = 0.;
        B0 = 0.4;
    
        [~,y] = ode45(@(t,y) get_odetype(t,y,s,odetype),t1,[Ts0,Tr0,C0,B0]);
        Ts = y(:,1);
        Tr = y(:,2);
        C = y(:,3);
        B = y(:,4);
        

        % day = 7
        t2 = tspan(index4+1:index6);
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = cart_dose;
        B0 = B(end);
        
        [~,y] = ode23s(@(t,y) get_odetype(t,y,s,odetype),t2,[Ts0,Tr0,C0,B0],opts);
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = [B ; y(:,4)];
        
        % day = 9, inject CART-cells
        t3 = tspan(index6+1:end);
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = C(end) + cart_dose;
        B0 = B(end);

        [~,y] = ode23s(@(t,y) get_odetype(t,y,s,odetype),t3,[Ts0,Tr0,C0,B0],opts);
    
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = [B ; y(:,4)];
        
        %disp([length(t1),length(t2),length(t3)])

        tumor_total = Ts(indices)+Tr(indices);
    
    else
        disp("enter a valid ode type")
    end

end



function dydt = get_odetype(t,y,param,odetype)
    if odetype == "Ts"
        Ts = y(1);
        r1 = param(1);
        K1 = param(2);
        dydt = r1*Ts*(1-Ts/K1); %Ts

    elseif odetype == "TsTr"
        dydt = zeros(2,1);
        Ts = y(1);Tr = y(2);
        %already optimized
        r1 = 0.18;K1 = 5.1e3;
        %to be optimized
        r2 = param(1);
        
        %ode sytems
        dydt(1) = r1*Ts*(1 - (Ts + Tr)/K1); %TS
        dydt(2) = r2*Tr*(1 - (Ts + Tr)/K1); %TR
    elseif odetype == "TsTrC"
        dydt = zeros(3,1);
        Ts = y(1);Tr = y(2);C = y(3);
        %already optimized
        r1 = 0.18;K1 = 5.1e03;r2 = 0.21;
        %to be optimized
        l = param(1);dC = param(2);muC = param(3);
        %default
        s = 3.05e-1; gammaC = 2.93e-2; k = 2.019e-7;
        wC = 3e-5; K2 = 1.65e3;
        
        %ode sytems
        DC = dC*(C/Ts)^l/(s+(C/Ts)^l);
        dydt(1) = r1*Ts*(1-(Ts+Tr)/K1) - DC*Ts; %Ts
        dydt(2) = r2*Tr*(1-(Ts+Tr)/K1); %Tr
        dydt(3) = - gammaC*C -muC*log((5e-2 + C)/K2)*(DC^2/(k+DC^2))*C...
            - wC*C*Ts;%C
            
    elseif odetype == "TsTrCB"
        dydt = zeros(4,1);
        Ts = y(1);Tr = y(2);C = y(3);B = y(4); 
        %default or alredy optimized(!) parameters here
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
        %====optimize dB and jB====%
        dB = param(1);
        muB = param(2);
        b = 5e-2;
        gammaB = 2e-2;       
        wB = 3.42e-6;
        
    
        DC = dC*(C/Ts)^l/(s+(C/Ts)^l);
        DS = dB*(B/Ts)^l/(s+(B/Ts)^l);
        DR = dB*(B/Tr)^l/(s+(B/Tr)^l);

        dydt(1) = r1*Ts*(1-(Ts+Tr)/K1) - DC*Ts - DS*Ts; %Ts
        dydt(2) = r2*Tr*(1-(Ts+Tr)/K1)- DR*Tr; %Tr
        dydt(3) = - gammaC*C -muC*log((B + C)/K2)*(DC^2/(k+DC^2))*C - wC*C*Ts;%C
        dydt(4) = b - gammaB*B - muB*log((C+B)/K2)*DS^2/(k+DS^2)*C-wB*(Tr+Ts);%B

    else
        disp("enter a valid ode type")
    end
end
