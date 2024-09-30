% Generate virtual patient corresponding to parameters r1, r2 and K1
% Note that We accept a parameter if it results in relative
% error less than "CRITERIA" for each sp. Uncomment stuff below to save
% the generated data

clearvars;clc;close all;
rng('default');
time_max = 30;dt = 0.1; tspan = 3:dt:time_max;
%test N number of parameter set and accept them based on "criteria"
N = 1000;
vp = GetUniformSample(N);
%this is very important. 
CRITERIA = 0.5;
sp = [1.,0.9,0.75];


tic;
progress = waitbar(0,'starting...');

% file2write = 'meta_data/vp_control.csv';
% if exist(file2write, 'file')==2
%   delete(file2write);
% end


f1 = figure('Name','100%');
f2 = figure('Name','90%');
f3 = figure('Name','75%');
figures = [f1,f2,f3];
c = 0;


data = zeros(6,3);
for k=1:length(sp)
    [exp_days,tumor_size] = get_data(sp(k),"control");
    data(:,k)=tumor_size; 
end

tic;
for i=1:N
    waitbar(i/N,progress,sprintf('%0.2f %%',100*(1-i/N)));
    T = zeros(6,length(sp));
    mre = zeros(6,length(sp));
    k = 0;
    for s_rate=sp
        k=k+1;
        s = vp(i,:);
        [tumor_total,~,~] = solve_ode(s,tspan,exp_days,s_rate);
        mre(:,k) = abs((data(:,k) - tumor_total)./data(:,k));
        T(:,k) = tumor_total;
    end
    
    if max(max(mre))<CRITERIA
        c=c+1;
        fprintf('%i-th patient mre:%0.2f\n',c,max(max(mre)))
        for j=1:length(figures)
            set(0, 'CurrentFigure', figures(j))
            plot(exp_days, T(:,j),'r')
            hold on;
        end
        % vp_control = array2table(s,"VariableNames",["r1","r2","K1"]);
        % writetable(vp_control,file2write,'WriteMode','append');
    end
end

title_names = {'antigen-positive ratio = 100%', ...
    'antigen-positive ratio = 90%', ...
    'antigen-positive ratio = 75%'};
for j=1:length(figures)
    [~,data] = get_data(sp(j),"control");
    set(0, 'CurrentFigure', figures(j))
    plot(exp_days, data,'ks-','MarkerSize',8,'MarkerFaceColor','k')
    xlabel('time(days)',FontSize=20)
    ylabel('tumor volume(mm^3)',FontSize=20)
    title(title_names{j}, 'FontSize', 18)
    grid on;
    %saveas(gcf,sprintf("figures/vp_control_test1_%0.0f.jpg",100*sp(j)),'jpg');
end
toc


delete(progress)
fprintf("detect: %0.4f",c/N)


function [tumor_total,Ts,Tr] = solve_ode(s,tspan,exp_days,sp)
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [~,indices] = ismember(exp_days,tspan);
    TS0 = sp*y0;
    TR0 = (1 - sp)*y0;
    [~,y] = ode23s(@(t,y) get_ode(t,y,s),tspan,[TS0,TR0],opts);
    Ts = y(:,1);
    Tr = y(:,2);
   
    tumor_total = Ts(indices)+Tr(indices);
end


function dydt = get_ode(t,y,s)
    dydt = zeros(2,1);
    Ts = y(1);Tr = y(2);
    r1 = s(1);
    r2 = s(2);
    K1 = s(3);
    dydt(1) = r1*Ts*(1-(Ts+Tr)/K1);%Ts    
    dydt(2) = r2*Tr*(1-(Ts+Tr)/K1);%Tr
end


function vp = GetUniformSample(N)
    %draw sensitive params from a normal dist
    r1 = Random(0.12,0.24,N);%
    r2 = Random(0.14,0.28,N);%
    K1 = Random(1e3,9e3,N);
    vp = horzcat(r1,r2,K1);

end
function z = Random(a,b,N)
    z = a+(b-a)*rand(1,N)';
end