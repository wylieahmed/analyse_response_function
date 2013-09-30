function [f,alp]=extract_respose_from_two_bead_fluctuation(data1,data2,k1,k2,srt)

%This will get the full response function, using the Kramers Kronig
%relation, and substract the traps influence
%You need to give it the fluctuation data of a two beads, the trap stiffness and the
%scanrate, the rest will be come automatically

%1 I will convert the data to a powerspectrum. 
%2 Then the fluct diss theorem will give the imaginary part of the response
%function
%3. Using the Kramers kronig relation we can then get the real part of the
%response function.
kT=1.38e-23*300;


%to get the powerspectrum we calculate the corremation functions and then
%do a fourier transform
[f,psd_1]=power_sd(data1,srt);
[f,psd_2]=power_sd(data2,srt);
[f,psd_d]=power_sd(data1-data2,srt);
[fc,psd_cross]=cross_power_sd(data1,data2,srt);


%fluct diss theorm gives the imaginary part of the the response function via:
chi_1=kramers_kronig_hil(i*pi*f./(kT).*psd_1);
chi_2=kramers_kronig_hil(i*pi*f./(kT).*psd_2);
chi_d=kramers_kronig_hil(i*pi*f./(kT).*psd_d);
chi_cross=kramers_kronig_hil(i*pi*f./(kT).*psd_cross);

%now we correct for the trap stiffness
alp=chi_cross./(1-k1*chi_1-k2*chi_2-k1*k2*chi_cross.^2+k1*k2*chi_1.*chi_2);





%now I use the method by Pietro
function chi=kramers_kronig_hil(chi)
chi_ii=imag(chi);
inter=-1*chi_ii(end:-1:1);
inter=[inter chi_ii];
hil=hilbert(inter);
hil=hil(end/2+1:end);
chi=imag(hil)+i*chi_ii;

