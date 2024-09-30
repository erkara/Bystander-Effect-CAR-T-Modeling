%solve system based on odetype and display the solution with data.
%beware that "odetype" and "datatype" must be compatible.

clearvars;clc;close all;
%time-parameters
time_max = 30;dt = 0.05;tspan = 3:dt:time_max;
%translated from the actual experiment
cart_dose = 0.1;

%odetype: ["Ts","TsTr","TsTrC","TsTrCB"]
%datatype: control,treatment,bystander


odetype = "TsTrC";
sp = .9; datatype="treatment";
[exp_days,data] = get_data(sp,datatype);

%solve ode and report root-mean-square error
[tumor_total,Ts,Tr,C,B] = solve_ode(tspan,exp_days,cart_dose,sp,odetype);
%rmse = sqrt(sum((tumor_total - data).^2)/length(data));
relative_errors = abs((tumor_total(1:end) - data(1:end)) ./ data(1:end));
mre = mean(relative_errors);
fprintf('\nmre: %0.2f\n',mre);

%plot the results
plot(tspan,Ts+Tr,'r-','LineWidth',3)
hold on;
plot(exp_days,data,'ks','MarkerSize',13,'MarkerFaceColor','k')
ax = gca;
ax.FontSize = 12;
figure_name = sprintf('figures//%s_%i.png',datatype,100*sp);

if strcmp(odetype, 'TsTrC') || strcmp(odetype, 'TsTrCB')
    index4 = find(tspan == 4);
    index6 = find(tspan == 6);
    yValueAt4 = Ts(index4) + Tr(index4);
    yValueAt6 = Ts(index6) + Tr(index6);
    plot([4, 6], [0, 0], 'bd', 'MarkerSize', 12,'MarkerFaceColor','b'); 
    legend({'model','data','{CAR-T injection}'},Location="best",FontSize=15)
else
    legend({'model','data'},Location="best",FontSize=15)
end
xticks(1:2:time_max);
xlim([0,max(exp_days)+2]);
xlabel('time(days)',FontSize=20)
ylabel('tumor volume(mm^3)',FontSize=20)
title(sprintf('antigen-positive ratio = %i%%',100*sp),'FontSize', 20)
grid on;
exportgraphics(ax,figure_name,'Resolution',450)

% figure;
% plot(tspan,Ts,'r','LineWidth',2)
% hold on;
% plot(tspan,Tr,'b','LineWidth',2)
% hold on;
% plot(tspan,C,'g','LineWidth',2)
% hold on;
% plot(tspan,B,'m','LineWidth',2)
% grid on;
% legend('Ts','Tr','C','B')





function [tumor_total,Ts,Tr,C,B] = solve_ode(tspan,exp_days,cart_dose,sp,odetype)
    %we return Ts,Tr,C,B all but some are just dummies in many case.
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [~,indices] = ismember(exp_days,tspan);
    
    %tumor grows until day-4
    index4 = find(tspan==4);
    t1 = tspan(1:index4);
    
    %first cart injection at day-6
    index6 = find(tspan==6);
    t2 = tspan(index4+1:index6);

    %second cart injection at day-
    t3 = tspan(index6+1:end);


    if odetype=="Ts"  
        [~,y] = ode45(@(t,y) get_odetype(t,y,odetype),tspan,y0,opts);
        Ts = y(:,1);
        Tr = 0;
        C = 0;
        B = 0;
        tumor_total = Ts(indices);
    
    elseif odetype=="TsTr"
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        [~,y] = ode45(@(t,y) get_odetype(t,y,odetype),tspan,[Ts0,Tr0],opts);
        Ts = y(:,1);
        Tr = y(:,2);
        C = 0;
        B = 0;
        tumor_total = Ts(indices)+Tr(indices);

    elseif odetype=="TsTrC"
        
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        C0 = 0.;
       
        [~,y] = ode45(@(t,y)get_odetype(t,y,odetype),t1,[Ts0,Tr0,C0]);
        Ts = y(:,1);
        Tr = y(:,2);
        C = y(:,3);
        
        % day = 4 first injection
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = cart_dose;
        
        [~,y] = ode15s(@(t,y)get_odetype(t,y,odetype),t2,[Ts0,Tr0,C0],opts);
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        
        %day = 6, second injection
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = C(end) + cart_dose;
     
        [~,y] = ode15s(@(t,y) get_odetype(t,y,odetype),t3,[Ts0,Tr0,C0],opts);
       
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = 0;
        
        tumor_total = Ts(indices)+Tr(indices);
   
    elseif odetype=="TsTrCB"
        
        Ts0 = sp*y0;
        Tr0 = (1 - sp)*y0;
        C0 = 0.;
        B0 = 0.1;
    
        [~,y] = ode45(@(t,y) get_odetype(t,y,odetype),t1,[Ts0,Tr0,C0,B0]);
        Ts = y(:,1);
        Tr = y(:,2);
        C = y(:,3);
        B = y(:,4);
        
        % first injection day=4
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = cart_dose;
        B0 = B(end);
        
        [~,y] = ode23s(@(t,y) get_odetype(t,y,odetype),t2,[Ts0,Tr0,C0,B0],opts);
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = [B ; y(:,4)];
        
        % second injection day=9
        Ts0 = Ts(end);
        Tr0 = Tr(end);
        C0 = C(end) + cart_dose;
        B0 = B(end);
        [~,y] = ode23s(@(t,y) get_odetype(t,y,odetype),t3,[Ts0,Tr0,C0,B0],opts);
    
        Ts = [Ts ; y(:,1)];
        Tr = [Tr ; y(:,2)];
        C = [C ; y(:,3)];
        B = [B ; y(:,4)];
      
        tumor_total = Ts(indices)+Tr(indices);
    
    else
        disp("enter a valid ode type")
    end

end

function dydt = get_odetype(t,y,odetype)
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
       
    if odetype == "Ts"
        Ts = y(1);
        dydt = r1*Ts*(1-Ts/K1); %Ts

    elseif odetype == "TsTr"
        dydt = zeros(2,1);
        Ts = y(1); Tr = y(2);
        dydt(1) = r1*Ts*(1 - (Ts + Tr)/K1); %Ts
        dydt(2) = r2*Tr*(1 - (Ts + Tr)/K1); %Tr
    elseif odetype == "TsTrC"
        dydt = zeros(3,1);
        Ts = y(1); Tr = y(2); C = y(3);
        %ode sytems
        DC = dC*(C/Ts)^l/(s+(C/Ts)^l);
        dydt(1) = r1*Ts*(1-(Ts+Tr)/K1) - DC*Ts; %Ts
        dydt(2) = r2*Tr*(1-(Ts+Tr)/K1); %Tr
        dydt(3) = - gammaC*C -muC*log((b + C)/K2)*(DC^2/(k+DC^2))*C - wC*C*Ts;%C
            
    elseif odetype == "TsTrCB"
        dydt = zeros(4,1);
        Ts = y(1); Tr = y(2); C = y(3); B = y(4); 
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







