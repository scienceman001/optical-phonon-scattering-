%Program caclulates reduced fermi energy  For InAs
%Program written by Rezo Kobaidze
%fermi_integral_calculator3(q,k) is in external file it calculates integral
%for any q and k values
h_pl=6.626068e-27;
pi=3.14159265359;
mass=0.025*9.108e-28;   %InAs effective mass
k_bol=1.38065e-16;
epsilon=14.5;
e=4.803e-10;
Q=336.23;  %debey temperature for InAs
%{
my_temp=fopen("temperature_N2_dasxivebuli.dat", "r");
my_conc=fopen("concentration_N2_dasxivebuli.dat", "r");
my_qsi=fopen("qsi_N2_dasxivebuli.dat","w");
my_z=fopen("z_N2_dasxivebuli.dat","w");
%}

my_temp=fopen("temperature_N1_dasxivebuli.dat", "r");
my_conc=fopen("concentration_N1_dasxivebuli.dat", "r");
my_qsi=fopen("qsi_N1_dasxivebuli.dat","w");
my_z=fopen("z_N1_dasxivebuli.dat","w");

% my_temp=fopen("temperature_InAs_digit3.dat", "r");
% my_conc=fopen("concentration_N1_dasxivebuli.dat", "r");
% my_qsi=fopen("qsi_InAs_digit3.dat","w");
% my_z=fopen("z_InAs_digit1.dat","w");


while ~feof(my_temp)
T=str2double(fgets(my_temp));
%calculating z parameter
z_par=Q/T;
fprintf(my_z,'%e\n', z_par);
%%%%%%%%%%%%%
 concentration=str2double(fgets(my_conc));
 concentration=concentration*1e+16;

f_0_5=(concentration*(h_pl^3))/(4*pi*(2*mass*k_bol*T)^1.5);


%calculating degeneracy parameter using bisection method
approach=0.001;
ab=[-100000,100000];  %initial values of [a,b]
key=1;
for q=ab
    f(key)=fermi_integral_calculator3(q,0.5)-f_0_5;
    key=key+1;
end %end of for loop

prod=f(1)*f(2);
if prod<0
    while abs(ab(1)-ab(2))>=approach
        c=(ab(1)+ab(2))/2;
key=1;
for q=ab
    f(key)=fermi_integral_calculator3(q,0.5)-f_0_5;
      key=key+1;
end %end of for loop
f_c=fermi_integral_calculator3(c,0.5)-f_0_5;
prod=f(1)*f_c;
if prod<0 
    ab(2)=c;
else
    ab(1)=c;
end
    end  %end of while loop
    qsi=c
   fprintf(my_qsi,'%f\n', qsi);
else
   "enter different initial values"
end %end of prod<0
end %end of   my_temp
