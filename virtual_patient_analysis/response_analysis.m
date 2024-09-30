% This code loops through each virtual parameter set cretated in
% meta_data/test1 and meta/test2 and classify each outcome as Responder(R),
% Non-Responder(NR) or Partial-Responder(PR) with the final profiles
% For one case (mean trajectories of each cell group), we need the full
% trajectory per vp. Uncomment line:57 to achive this but this will lead to
% a very slow code and large dataset. Outcome for test-1 case is located at
% meta_data/test1/vp_responses_test1_all_vars.csv.

clear all;clc;close all;
%random seed
rng('default');
time_max = 30; dt = 0.01; tspan = 3:dt:time_max;
cart_dose = 0.1;
vp = table2array(readtable("meta_data/test1/vp_full_test1.csv"));
%vp = table2array(readtable("meta_data\test2\vp_full_test2.csv"));
N = size(vp,1);

lower_criteria = 2;
upper_criteria = 50;
progress = waitbar(0,'starting...');


file2write = 'meta_data/test1/vp_responses_test1.csv';
%file2write = 'meta_data/test2/vp_responses_test2.csv';
if exist(file2write, 'file')==2
  delete(file2write);
end
%only for 75%
sp = 0.75;
R_legend_added = false;
PR_legend_added = false;
NR_legend_added = false;


for i=1:10
    waitbar(i/N,progress,sprintf('%0.2f %%',100*(1-i/N)));
    params =  vp(i,:);
    [Ts,Tr,C,B] = solve_ode(params,tspan,cart_dose,sp);
    T = real(Ts+Tr);
    T_final = T(end);

    if max(T_final) < lower_criteria
        response = "R";
    elseif lower_criteria<max(T_final) & max(T_final)<upper_criteria
        response ="PR";
    else
        response="NR";
    end
    
    % fprintf('patient:%i Tmax:%0.1f response: %s\n',i,max(T_final),response)
    
    %%% write only params and response %%%
    varnames = ["r1","r2", "K1", "l","dC","muC","s","gammaC","k",...
                "wC","K2","dB","muB","b","gammaB","wB","response"];
    data2write = [params,response];
    vp_bystander = array2table(data2write,"VariableNames",varnames);

    %% write all trajectories, very slow %%
    %data2write = [params, response,real(T_final),real(Ts'),real(Tr'),real(C'),real(B')];
    %vp_bystander = array2table(data2write);
    writetable(vp_bystander,file2write,'WriteMode','append');


end
% legend(h, labels);
% ax = gca;
% ax.FontSize = 12;
% xlabel('time(days)',FontSize=12)
% ylabel('tumor volume(mm^3)',FontSize=12)
% exportgraphics(ax,'figures/base_response_test2.jpg','Resolution',300) 
delete(progress)



function [Ts,Tr,C,B] = solve_ode(params,tspan,cart_dose,sp)
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    
    index4 = find(tspan==4);
    index6 = find(tspan==6);
    %tumor grows without intervention until day 4
    t1 = tspan(1:index4);
    Ts0 = sp*y0;
    Tr0 = (1 - sp)*y0;
    C0 = 0.;
    B0 = 0.1;

    [~,y] = ode45(@(t,y) get_ode(t,y,params),t1,[Ts0,Tr0,C0,B0]);
    Ts = y(:,1);
    Tr = y(:,2);
    C = y(:,3);
    B = y(:,4);
    

    % day = 6
    t2 = tspan(index4+1:index6);
    Ts0 = Ts(end);
    Tr0 = Tr(end);
    C0 = cart_dose;
    B0 = B(end);
    
    [~,y] = ode23s(@(t,y) get_ode(t,y,params),t2,[Ts0,Tr0,C0,B0],opts);
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

    [~,y] = ode23s(@(t,y) get_ode(t,y,params),t3,[Ts0,Tr0,C0,B0],opts);

    Ts = [Ts ; y(:,1)];
    Tr = [Tr ; y(:,2)];
    C = [C ; y(:,3)];
    B = [B ; y(:,4)];
   
end



