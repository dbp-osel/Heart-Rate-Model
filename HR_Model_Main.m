%{
Authors 
Varun Kanal
Christopher Scully
Ramin Bighamian

For questions, contact ramin.bighamian@fda.hhs.gov

The HR model is described and published in:
Varun Kanal, Pras Pathmanathan, Jin-Oh Hahn, George Kramer, 
Christopher Scully, and Ramin Bighamian, Development and Validation of a 
Mathematical Model of Heart Rate Response to Fluid Perturbation, 
Scientific Reports 12, 21463 (2022).

If you found this software useful, please consider citing our publication above. 

Public domain license

FDA Software Disclaimer: 
This software and documentation (the "Software") were developed at the Food
and Drug Administration (FDA) by employees of the Federal Government in the
course of their official duties. Pursuant to Title 17, Section 105 of 
the United States Code, this work is not subject to copyright protection 
and is in the public domain. Permission is hereby granted, free of charge,
to any person obtaining a copy of the Software, to deal in the Software 
without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, or sell copies of the 
Software or derivatives, and to permit persons to whom the Software is 
furnished to do so. FDA assumes no responsibility whatsoever for use by 
other parties of the Software, its source code, documentation or compiled 
executables, and makes no guarantees, expressed or implied, about its 
quality, reliability, or any other characteristic. Further, use of this 
code in no way implies endorsement by the FDA or confers any advantage in 
regulatory decisions. Although this software can be redistributed and/or 
modified freely, we ask that any derivative works bear some notice that 
they are derived from it, and any modified versions bear some notice that 
they have been modified.
%}

clear
clc
close all

global Hemo_Time Infusion_inp Hemorrhage_inp Urine_inp  SamplingTime Hemo_HR 


load('Example_Data') % load data from one individual
Hemo_Time = Example_Data{1,1}(:,1); % time instants for data 
Hemo_HR = Example_Data{1,1}(:,2); % HR data
Infusion_inp = Example_Data{1,2}(:,1); %infusion profile
Hemorrhage_inp = Example_Data{1,2}(:,2); %hemorrhage profile
Urine_inp = Example_Data{1,2}(:,3); %urine profile
SamplingTime = Example_Data{1,3}(:,1); %sampling time
simulation_time = [0:SamplingTime:180]';
%%
lowlimit = [0 0 0 0 0 0 0 40];
uplimit = [5 5 2 2 1 1 0.1 150];
f = [];
for ii = 1:10 % model calibration
    ii
x0 = lowlimit+(uplimit-lowlimit).*rand(1,length(lowlimit));
options = optimset('Algorithm','interior-point'); 
[x(ii,:),f(ii,1)] = fmincon(@HR_Optimizer,x0,[],[],[],[],lowlimit ...
    ,uplimit,[],options);
end
f_opt = min(f) % calibration cost function
x_opt = x(find(f==min(f),1),:) % identified parameters
%%
G1=x_opt(1);G2=x_opt(2);pow1=x_opt(3);pow2=x_opt(4);G3=x_opt(5);Kp=x_opt(6);Ki=x_opt(7);HR0=x_opt(8); 
sim('HR_Model'); 

figure(1)
subplot(221)
plot(simulation_time,Infusion_inp*1000,'LineWidth',2)
hold on
plot(simulation_time,-Hemorrhage_inp*1000,'--r','LineWidth',2)
set(gca,'XTick',[0:30:180])
xlim([0 180])
ylim([-100 80])
ylabel('Fluid & Hemorrhage [ml/min]','fontsize',12,'fontweight','b')
set(gca,'XTick',[0:30:180],'FontWeight','bold')
set(gca,'YTick',[-100:20:80],'FontWeight','bold')
legend('Fluid','Hemorrhage','Location','SouthEast')
grid on

subplot(222)
plot(Hemo_Time,Hemo_HR,'ro')
hold on
plot(simulation_time,HR_sol)
xlabel('Time [min]')
ylabel('HR [bpm]')
set(gca,'XTick',[0:30:180],'FontWeight','bold')
set(gca,'YTick',[50:25:250],'FontWeight','bold')
legend('True','Estimated','Location','SouthEast')
ylim([50 250])
grid on