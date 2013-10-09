%%%%%%%%%%%%%%%%%%%%%%%
% Written by Wylie Ahmed (01/10/2013)
% This script is for plotting the corrected_response.mat data generated
% from GUI_intracellular_response_function_v2.m written
%


%% Load Data 

addpath /Users/wylieahmed/Documents/ResearchBySemester/Fall_2013/Matlab/021013_Oocytes



%% Plot Complex moduli

figure
for i = 1:17
loglog(corrected_response(i).response_f',abs(real(corrected_response(i).G)),'xr')
hold all
loglog(corrected_response(i).response_f',abs(imag(corrected_response(i).G)),'ob')
loglog([100:1000],[100:1000].^.75*.1,'k');
legend('elastic','dissipative')
%hold off
xlabel('Frequency (Hz)')
ylabel('G (Pa)')
% pause
end

% Averaging of values over the entire experiment

avgRealG = 


%% Plot Response functions

figure
for i = 1:17
loglog(corrected_response(i).response_f,corrected_response(i).corr_active_response,'-xr')
hold on
loglog(corrected_response(i).psd_f,corrected_response(i).passive_response,'ob')
legend('Active Response','Passive Response')
hold off
xlabel('Frequency (Hz)')
ylabel('Response Function (m/N)')
pause
end