%this program caclulates fermi integral for any k an q value using Gauss-Legendre method
%program written by Rezo Kobaidze
function fermi_integral=fermi_integral_calculator3(q,k)
my_c=fopen("gauss_legandre_weights_c.dat","r");
my_x=fopen("gauss_legandre_weights_x.dat", "r");
%reset parameters
int_1=0;
int_2=0;
int1=0;
int2=0;

while ~feof(my_c)
    c=str2double(fgets(my_c));
    x=str2double(fgets(my_x));
    A=0.5+0.5*x;
    f1=((A^k))/(exp(A-q)+1);
    int_1=c*0.5*f1;
    int1=int1+int_1;
    
    A2=1/A;
    f2=((1/A)^2)*((A2^k))/(exp(A2-q)+1);
    int_2=c*0.5*f2;
    int2=int2+int_2;
    fermi_integral=int1+int2;
end
fclose(my_c);
fclose(my_x);
end%end of function