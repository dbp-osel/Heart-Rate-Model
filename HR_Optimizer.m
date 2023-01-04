function [f] = HR_Optimizer(x)

global Hemo_Time Infusion_inp Hemorrhage_inp Urine_inp  SamplingTime Hemo_HR 

G1=x(1);G2=x(2);pow1=x(3);pow2=x(4);G3=x(5);Kp=x(6);Ki=x(7);HR0=x(8); 

options = simset('SrcWorkspace','current');
sim('HR_Model',[],options)

associated_response = interp1(time,HR_sol,Hemo_Time);
error = [Hemo_HR-associated_response];% HR error

f =normlike([0,std(error)],error)+2*0.5*norm((x(8)-Hemo_HR(1)),2)+2*0.1*norm(x(1:7),2); % cost function
