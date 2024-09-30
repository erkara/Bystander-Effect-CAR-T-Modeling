% Use this function to access the clinical data we extracted 
% from the experimental paper. "data" is the volume of the tumor size(mm^3)
% measured on "days". Meausurement days are slighly different for different 
% experiments.
% "sp":  percantage of antigen positive cell population for; 
% "control": control group
% "treatment": only CAR-T cell treatment
% "bystander": CAR-T + CTX inducing bystander effect



function [days,data]  = get_data(sp,datatype)
    data_100 = table2array(readtable("data/AE17_100.csv"));
    data_mixed = table2array(readtable("data/AE17.csv"));
    data_ctx = table2array(readtable("data/ctx.csv"));
    if datatype=="control"
        days = data_mixed(1:6,1);
        if abs(sp-0.5)<1e-10
            data = data_mixed(1:6,2);
        elseif abs(sp-0.75)<1e-10
            data = data_mixed(1:6,4);
        elseif abs(sp-0.9)<1e-10
            data = data_mixed(1:6,6);
        elseif abs(sp-1)<1e-10
            days = data_100(1:6,1);
            data = data_100(1:6,2);
        else
            disp('enter a valid percentage')
        end
    elseif datatype=="treatment"
        days = data_mixed(1:7,1);
        if abs(sp - 0.5)<1e-10
            data = data_mixed(1:7,3);
        elseif abs(sp - 0.75)<1e-10
            data = data_mixed(1:7,5);
        elseif abs(sp - 0.9)<1e-10
            data = data_mixed(1:7,7);
        elseif abs(sp - 1)< 1e-10
            data = data_100(1:7,3);
        end
    elseif datatype=="bystander"
        days = data_ctx(:,1);
        if abs(sp - 0.75)<1e-10
            data = data_ctx(:,2);
        elseif abs(sp - 0.9)<1e-10
            data = data_ctx(:,3);
        end
    end
end