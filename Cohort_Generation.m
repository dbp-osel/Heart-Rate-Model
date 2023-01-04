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

load('Calibrated_Parameters') % load parameters identified for 21 subjects
Par_aux(15,:) = [];% remove the test subjct which is suject # 15 for leave-one-out analysis


load('Example_Data') % load data for the test subject
Hemo_Time = Example_Data{1,1}(:,1); % time instants for data 
Hemo_HR = Example_Data{1,1}(:,2); % HR data
Infusion_inp = Example_Data{1,2}(:,1); %infusion profile
Hemorrhage_inp = Example_Data{1,2}(:,2); %hemorrhage profile
Urine_inp = Example_Data{1,2}(:,3); %urine profile
SamplingTime = Example_Data{1,3}(:,1); %sampling time
simulation_time = (0:SamplingTime:180)';
%% Compartment method for virtual cohort generation
[A B C]=ndgrid(1:size(Par_aux,1),1:size(Par_aux,1),1:size(Par_aux,1)); % different combinations of infusion, hemorrhage, and controller compartments
A1=reshape(A,[],1);
B1=reshape(B,[],1);
C1=reshape(C,[],1);
d = [];
d=[A1,B1,C1];
mat_all2 = [];mat_all3 = [];
for ii = 1:length(d)
    mat_all2(ii,:) = [Par_aux(d(ii,1),1) Par_aux(d(ii,2),2) Par_aux(d(ii,1),3) Par_aux(d(ii,2),4:5) Par_aux(d(ii,3),6:7)]; % mixing virtual subject combinations
    mat_all3(ii,:) = (Par_aux(d(ii,1),:)+Par_aux(d(ii,2),:)+Par_aux(d(ii,3),:))/3;% average virtual subject combinations
end

mat_all = [mat_all2;mat_all3];% all possible virtual subject parameter combinations
n_size = length(mat_all);% size of all possible virtual subjects to run

%%% Run simulations below
%%
Physio_subject_count = 0;
HR0 = Hemo_HR(1);% start from the initial heart rate. Uncertainty can be added. 
for ii = 1:n_size
    ii
    G1=mat_all(ii,1);G2=mat_all(ii,2);pow1=mat_all(ii,3);pow2=mat_all(ii,4);G3=mat_all(ii,5);Kp=mat_all(ii,6);Ki=mat_all(ii,7);% model parameters for each virtual subject
    sim('HR_Model');  
    
    if min(HR_sol)>40 && max(HR_sol)<250 % check the criteria for physilogical simulations
        Physio_subject_count = Physio_subject_count+1;% counting the number of physilogical simulations
        response_HR_physio(Physio_subject_count,:) = interp1(simulation_time,HR_sol,Hemo_Time);%down sampled HR data for physiological subjects
        physio_sample(Physio_subject_count)= ii; 
        NRMSE(Physio_subject_count,1) = rms(Hemo_HR-(response_HR_physio(Physio_subject_count,:)'))/mean(Hemo_HR);% normalized root mean square error
    end 
end
%% Plot the prediction envelope
figure(1)
subplot(221)
plot(simulation_time,Infusion_inp*1000,'LineWidth',2)% plot infusion profile
hold on
plot(simulation_time,-Hemorrhage_inp*1000,'--r','LineWidth',2)% plot hemorrhage profile
set(gca,'XTick',[0:30:180])
xlim([0 180])
ylim([-100 80])
ylabel('Fluid & Hemorrhage [ml/min]','fontsize',12,'fontweight','b')
set(gca,'XTick',[0:30:180],'FontWeight','bold')
set(gca,'YTick',[-100:20:80],'FontWeight','bold')
legend('Fluid','Hemorrhage','Location','SouthEast')
grid on
box on

subplot(222)
hold on
plot(Hemo_Time,Hemo_HR,'vb','LineWidth',1.5,'MarkerSize',7)
grid on
xlim([0 180])
set(gca,'XTick',[0:30:180],'FontWeight','bold')

% Filter out data that meets the NRMSE criteria
Relevant_subjects = find(NRMSE<=0.2);% identify relevant virtual subjects with NRMSE<20%
response_HR_relevant = response_HR_physio(Relevant_subjects,:);% down sampled predction data for relevant subjects

P_envelope = [];
first_percentile = [];
second_percentile = [];
for jj = 1:length(Hemo_HR) % compute 95th percentile envelope
    first_percentile(jj,1) = prctile(response_HR_relevant(:,jj),2.5); 
    second_percentile(jj,1) = prctile(response_HR_relevant(:,jj),97.5);
end
P_envelope = [first_percentile second_percentile]; % 95th percentile envelope

hold on
P_envelope_plot=[P_envelope(:,1)' flipud(P_envelope(:,2))'];
fill([Hemo_Time' fliplr(Hemo_Time')], P_envelope_plot, 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4);% plot 95th percentile prediction envelope
grid on
xlim([0 180])
set(gca,'XTick',[0:30:180],'FontWeight','bold')