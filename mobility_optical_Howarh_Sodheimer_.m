%this program caclulates Optical scattering mobility for InAs 
%Program written by Rezo Kobaidze
%this program uses external file fermi_integral_caclulator3(q,k)
%formula (43) and (69) in Howarh and Sodheimer Semi-Conductors
% The Theory of Electronic Conduction in Polar

% my_z=fopen("z_N2_dasxivebuli.dat","r");
% my_temp=fopen("temperature_N2_dasxivebuli.dat", "r");
% my_conc=fopen("concentration_N2_dasxivebuli.dat", "r");
% my_qsi=fopen("qsi_N2_dasxivebuli.dat","r");
% my_G=fopen("G_N2_dasxivebuli.dat", "w");
% my_mobility=fopen("mobility_N2_dasxivebuli(optical).dat","w");
% my_fermiInt_05=fopen("FermiIntegral05_N2_mobility.dat","w");
% my_fermiInt_1_5=fopen("FermiIntegral1_5_N2_mobility.dat","w");


my_z=fopen("z_N1_dasxivebuli.dat","r");
my_temp=fopen("temperature_N1_dasxivebuli.dat", "r");
my_conc=fopen("concentration_N1_dasxivebuli.dat", "r");
my_qsi=fopen("qsi_N1_dasxivebuli.dat","r");
my_G=fopen("G_N1_dasxivebuli.dat", "w");
my_mobility=fopen("mobility_N1_dasxivebuli(optical)Howarh and Sodheimer.dat","w");
my_fermiInt_05=fopen("FermiIntegral05_N1_mobility.dat","w");
my_fermiInt_1_5=fopen("FermiIntegral1_5_N1_mobility.dat","w");

%calculating  fermi integrals from external function
%fermi_integral_caclulator3(q,k)
while ~feof(my_temp)
   T=str2double(fgets(my_temp));
   con=str2double(fgets(my_conc));
   con=con*1e16;
   q=str2double(fgets(my_qsi)); %qsi parameter is taken from external program
   z=str2double(fgets(my_z));  %z parameter is taken from external program
  %%%calculating fermi integrals from external program fermi_integral_calculator3(q,k)
 %%%calculating fermi integrals from external program fermi_integral_calculator3(q,k)
k=[0.5,1.5];
key=1;
for pow=k
    pow;
fermi_integral(key)=fermi_integral_calculator3(q,pow);
key=key+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%calculating D00 value
int_1=0;
int_2=0;
int_d_1=0;
int_d_2=0;
D_00=0;
my_c=fopen("gauss_legandre_weights_c.dat","r");
my_x=fopen("gauss_legandre_weights_x.dat", "r");
while ~feof(my_c)
        c=str2double(fgets(my_c));
		x=str2double(fgets(my_x));
        A=0.5+0.5*x; %caclulating functions argument
		B=(exp(q-A)+1)*(exp(A-q)+exp(-z));
		d1=2*sqrt(A*(A+z))/B;   
        int_1=c*0.5*d1; 
        int_d_1=int_d_1+int_1; %the first part of integral
        
        A2=1/A;
		B2=(exp(q-A2)+1)*(exp(A2-q)+exp(-z));
		d2=((1/A)^2)*2*sqrt(A2*(A2+z))/B2;
		int_2=c*0.5*d2;
		int_d_2=int_d_2+int_2; %the second part of integral
		D_00=int_d_1+int_d_2;
        end %end of  my_c
        fclose(my_c);
        fclose(my_x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caclulating D01 value
int_1=0;
int_2=0;
int_d_1=0;
int_d_2=0;
D_01=0;
my_c=fopen("gauss_legandre_weights_c.dat","r");
my_x=fopen("gauss_legandre_weights_x.dat", "r");
while ~feof(my_c)
		c=str2double(fgets(my_c));
		x=str2double(fgets(my_x));

		A=0.5+0.5*x; % caclulating functions argument
		B=(exp(q-A)+1)*(exp(A-q)+exp(-z));
      %M=sqrt(A/z)+(2*A+z)*sqrt(A*(A+z)); % sinh argument
        M=sqrt(A/z); 
        
        d1=(((z^2)/sinh(M))+(2*A+z)*sqrt(A*(A+z)))/(B);%D01
        %d1=(z^2)/(B*sinh(M));  %D01
		int_1=c*0.5*d1;
		int_d_1=int_d_1+int_1;%calculates the first part of integral
        
		A2=1/A;
		B2=(exp(q-A2)+1)*(exp(A2-q)+exp(-z));
		%M2=sqrt(A2/z)+(2*A2+z)*sqrt(A2*(A2+z));
         M2=sqrt(A2/z);
        %A^2 rom gavyavi integrals amoxsnis xerxia
        d2=((1/A)^2)*((((z^2)/sinh(M2))+(2*A2+z)*sqrt(A2*(A2+z)))/(B2));%D01
		%d2=(z^2)/((A^2)*B2*sinh(M2));
		int_2=c*0.5*d2; 
		int_d_2=int_d_2+int_2; %calculates the second part of integral
		D_01=int_d_1+int_d_2; %combine both parts
end %end of my_c 
	fclose(my_c);
	fclose(my_x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caclulating D11 value
int_1=0;
int_2=0;
int_d_1=0;
int_d_2=0;
D_11=0;
my_c=fopen("gauss_legandre_weights_c.dat","r");
my_x=fopen("gauss_legandre_weights_x.dat", "r");
while ~feof(my_c)
		c=str2double(fgets(my_c));
		x=str2double(fgets(my_x));
        A=0.5+0.5*x;  %caclulating functions argument
		B=(exp(q-A)+1)*(exp(A-q)+exp(-z));
        %M=sqrt(A/z)+(A^1.5)*((A+z)^1.5);
        M=sqrt(A/z);
		%d1=2*(z^2)*(2*A+z)/(B*sinh(M));   
        d1=2*(((z^2)*(2*A+z)/(sinh(M)))+((A^1.5)*((A+z)^1.5)));   
		int_1=c*0.5*d1;
		int_d_1=int_d_1+int_1; %calculates the first part of integral
        
        A2=1/A;
		B2=(exp(q-A2)+1)*(exp(A2-q)+exp(-z));
		%M2=sqrt(A2/z)+(A2^1.5)*((A2+z)^1.5);
        M=sqrt(A2/z); 
		%d2=2*(z^2)*(2*A2+z)/((A^2)*B2*sinh(M2));
        d2=((1/A)^2)*(2*(((z^2)*(2*A2+z)/(sinh(M2)))+((A2^1.5)*((A2+z)^1.5))));
		int_2=c*0.5*d2;
		int_d_2=int_d_2+int_2;%calculates the second part of integral
		D_11=int_d_1+int_d_2;  %combining integrals and caclulating D_11
end %end of my_c
fclose(my_c);
fclose(my_x);
%%%%%%%%%%%%%%%%%%%%%%
%%%calculating G function
p1=25*(fermi_integral(2)^2)*D_00+9*(fermi_integral(1)^2)*D_11-30*fermi_integral(1)*fermi_integral(2)*D_01;
p2=4*(D_00*D_11-(D_01^2));
G=p1/p2;
fprintf(my_G,'%f\n', G); 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%caclulating Mobility  for InAs sample iridated by neutrons
%all values are in CGS
a=6.058e-8;
M=0.753e-22;
epsilon0=14.5; %low frequency dielectric constant
epsilon=11.6;% hight frequency dielecric constant
h_plank=6.626068e-27;
k_bol=1.38065e-16;
wl=4.4e+13;
e=4.803e-10;
A=(16*(a^3)*(M*wl)*(k_bol^2))/(300*3*(e^3)*con*(h_plank^2)*2*pi);
mobility=A*G*(T^2)*(exp(z)-1) 
fprintf(my_mobility,'%e\n', mobility); 
end  %end of  my_temp
