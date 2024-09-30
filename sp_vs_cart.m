% In this file, we investigate how Responder/Non-Responder categories 
% change for varying CAR-T infusions vs positive antigent rates.
clc;clearvars;
N = 50;
cart_doses = linspace(0.01,0.5,N);
sp_values = linspace(0.10,0.95,85);
m = length(cart_doses);
n = length(sp_values);
%response = zeros(m*n);
time_max = 30;dt = 0.1;tspan = 3:dt:time_max;
%initial seed is 2mm3 from mouse paper.
lower_criteria = 2;
upper_criteria = 50;

%save the tumor size on these days
checkpoint_days = [9,12,15,18,21];
day_indices = arrayfun(@(d) find(abs(tspan - d) < dt/2, 1), checkpoint_days);


%write to csv
headers = {'cart_dose', 'sp_value', 'response',...
    't9', 't12', 't15', 't18', 't21','t_final_30'};
file2write = 'meta_data/cart_sp_response.csv';
if exist(file2write, 'file')==2
  delete(file2write);
end
writecell(headers, file2write);

k=0;
progress = waitbar(0,'starting...');
for cart_dose=cart_doses
    for sp=sp_values
        k=k+1;
        waitbar(k/(m*n),progress,sprintf('%0.2f %%',100*(1-k/(m*n))));
        [Ts,Tr,C,B] = solve_ode(tspan,cart_dose,sp);
        T = Ts+Tr;
        T_final = T(end);
        T_days = T(day_indices);
        
        if max(T_final) < lower_criteria
            response = "R";
        elseif lower_criteria<max(T_final) && max(T_final)<upper_criteria
            response ="PR";
        else
            response = "NR";
        end
        
        writematrix([cart_dose,sp,response,T_days',T_final],file2write,...
            'WriteMode','append','Delimiter', 'comma');
    end
end
delete(progress)




function [Ts,Tr,C,B] = solve_ode(tspan,cart_dose,sp)
    %we return Ts,Tr,C,B all but some are just dummies in many case.
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    %[~,indices] = ismember(exp_days,tspan);
    
    %tumor grows until day-3, first car injection day-4
    index4 = find(tspan==4);
    t1 = tspan(1:index4);
    
    %second cart injection at day-6
    index6 = find(tspan==6);
    t2 = tspan(index4+1:index6);
    t3 = tspan(index6+1:end);
    
    Ts0 = sp*y0;
    Tr0 = (1 - sp)*y0;
    C0 = 0.;
    B0 = 0.1;


    %tumor grows from day-3 to day-4
    [~,y] = ode45(@(t,y) get_ode(t,y),t1,[Ts0,Tr0,C0,B0]);
    Ts = y(:,1);
    Tr = y(:,2);
    C = y(:,3);
    B = y(:,4);
    
    % first injection day=4
    Ts0 = Ts(end);
    Tr0 = Tr(end);
    C0 = cart_dose;
    B0 = B(end);
    
    [~,y] = ode23s(@(t,y) get_ode(t,y),t2,[Ts0,Tr0,C0,B0],opts);
    Ts = [Ts ; y(:,1)];
    Tr = [Tr ; y(:,2)];
    C = [C ; y(:,3)];
    B = [B ; y(:,4)];
    
    % second injection day=6
    Ts0 = Ts(end);
    Tr0 = Tr(end);
    C0 = C(end) + cart_dose;
    B0 = B(end);
    [~,y] = ode23s(@(t,y) get_ode(t,y),t3,[Ts0,Tr0,C0,B0],opts);

    Ts = [Ts ; y(:,1)];
    Tr = [Tr ; y(:,2)];
    C = [C ; y(:,3)];
    B = [B ; y(:,4)];
  
    %tumor_total = Ts(indices)+Tr(indices);

end



