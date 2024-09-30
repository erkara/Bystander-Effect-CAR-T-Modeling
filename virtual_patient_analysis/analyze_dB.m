%check if the percentage of responders go up with higher dB values
clear all;clc;close all;
%random seed
rng('default');
time_max = 30; dt = 0.1; tspan = 3:dt:time_max;
cart_dose = 0.1;
sp = 0.75;
vp = readtable("meta_data/vp_full_test1.csv");
%vp = readtable("meta_data/vp_full_test2.csv"); % test2
vp = table2array(vp(:,1:16));
N = size(vp,1);
%for testing on a small sample
M = N;

lower_criteria = 2;
upper_criteria = 50;
%boost dB--> 10%,20%,30%,40%
boost_percentage = transpose([1.0,1.1,1.2,1.3,1.4]);

progress = waitbar(0,'starting...');
total_iterations = M * length(boost_percentage);

% Rows for each boost_percentage and columns for counts
all_counts = zeros(length(boost_percentage), 3); 

for j = 1:length(boost_percentage) 
    boost = boost_percentage(j);
    counts = zeros(1,3); % Reset counts for each boost
    
    
    for i=1:M
       current_iteration = (j-1)*M + i;
       waitbar(current_iteration / total_iterations, progress, ...
           sprintf('Progress: %0.2f%%', 100 * current_iteration / total_iterations));
        % Boost dB, which is the 12th parameter
        vp(i,12) = boost * vp(i,12);
        params = vp(i,:);
        [Ts,Tr,~,~] = solve_ode(params,tspan,cart_dose,sp);
        T = real(Ts+Tr);
        Tmax = T(end);
        
        if Tmax < lower_criteria
            counts(1) = counts(1) + 1;
            response = "R";
        elseif lower_criteria < Tmax && Tmax < upper_criteria
            counts(2) = counts(2) + 1;
            response = "PR";
        else
            response = "NR";
            counts(3) = counts(3) + 1;
        end
    end
    total_responses = sum(counts);
    percentage_counts = counts / total_responses;
    all_counts(j,:) = percentage_counts;
end

% Writing to CSV
filename = 'meta_data/dB_boosts_test2.csv';
if exist(filename, 'file')==2
  delete(filename);
end

boost_percentage_str = arrayfun(@(x) sprintf('%d%%', round((x - 1) * 100)), ...
    boost_percentage, 'UniformOutput', false);
boost_percentage_str{1} = 'base'; % Replace the first element with 'base'

% Create a table with string boost_percentage and all_counts
tbl = table(boost_percentage_str, all_counts(:,1), all_counts(:,2), all_counts(:,3), ...
            'VariableNames', {'boost', 'R', 'PR', 'NR'});

writetable(tbl, filename);

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



