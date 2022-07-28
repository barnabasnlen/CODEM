exceldocument=input('enter the name of document excel containing data <= 5000 :');
[DATA,TILTE,DOCUMENT]=xlsread(exceldocument);
N=size(DATA,1);
m=input('enter the degree of polynomial separation:');
longx=input('enter the length in km along x:');
longy=input('enter the length in km along y:');
nx=input('enter the  number of data along x:');
ny=input('enter the  number of data along y:');
c=input('enter the acceptance criteria <= -0.04:');
pmax=input('enter the maximum depth of upward continuation:');
pas=input('enter the upward continuation step in km:');
ax=input('enter the minimum value of longitude in km:');
ay=input('enter the minimum value of latitude in km:');
n1=(pmax/pas)+1;
pasx=longx/(nx-1);
pasy=longy/(ny-1); 
X=zeros(N,1);	
Y=zeros(N,1);
B=zeros(N,1);
g7ix=zeros(nx,ny);
g7iy=zeros(nx,ny);
g7x=zeros(nx,ny);
g7y=zeros(nx,ny);
g71=zeros(N,1);
g72=zeros(N,1);
g=zeros(nx,ny);
g8=zeros(nx,ny);
g1=zeros(nx,ny);
g7=zeros(nx,ny);
g2=zeros(N,1);
g9=zeros(N,1);
g3=zeros(N,1);
g4=zeros(N,1);
head1={'Reg';'Res'}';
head2={'upwcon';'gradh';'xmax';'ymax';'gmax'}';
for i=1:N
    X(i)=DATA(i,1);
    Y(i)=DATA(i,2);
    B(i)=DATA(i,3);
end
m2=(m+1)*(m+2)/2;
G=zeros(m2,m2);
E=zeros(m2,1);
Re=zeros(N,1);
Ri=zeros(N,1);
t1=0;
for j=0:m
for l=j:-1:0
    t2=0;
    t1=t1+1;
for k=0:m
for s=k:-1:0
        t2=t2+1;
for i=1:N
    G(t1,t2)=G(t1,t2)+(X(i)^(l+s))*(Y(i)^(k+j-l-s));
end
end
end
for i=1:N
    E(t1)=E(t1) +B(i)*((X(i))^l)*(Y(i)^(j-l));
end
end
end
    C=G\E;
for i=1:N
    k=0;
for j=0:m
for l=j:-1:0
    k=k+1;
    Re(i)=Re(i)+C(k)*(X(i)^l)*(Y(i)^(j-l));
end
end
end
for i=1:N
    Ri(i)=B(i)-Re(i);
end
A=zeros(N,2);
for i=1:N
A(i,1)=Re(i,1);
A(i,2)=Ri(i,1);
end
xlswrite(exceldocument,A,'D2:E5000');
xlswrite(exceldocument,head1,'D1:E1');
k=1;
for i=1:nx
for j=1:ny
    g(i,j)= Ri(k);
    k=k+1;
end
end
gr=fft2(g,nx,ny);
xl=0;
for e=1:n1
u=xl;
for i=1:nx
for j=1:ny
        h=exp(-2*pi*u*sqrt(((i-1)/longx)^2+((j-1)/longy)^2));
        g8(i,j)=h;
end
end
gi=gr.*g8;
g5=ifft2(gi,nx,ny);
k=1;
for i=1:nx
for j=1:ny
        g2(k)=g5(i,j);
        k=k+1;
end
end
g6=real(g2);
k=1; 
for i=1:nx
for j=1:ny
       g7(i,j)=g6(k);
       g7x(i,j)=DATA(k,1);
       g7y(i,j)=DATA(k,2);
        k=k+1;
end
end
for i=2:nx-1
for j=2:ny-1
g7ix(i,j)= (g7(i-1,j)-g7(i+1,j))/(2*pasx);
g7iy(i,j)=(g7(i,j+1)-g7(i,j-1))/(2*pasy);
end
end
for i=2:nx-1
for j=2:ny-1
        g1(i,j)=sqrt(g7ix(i,j)^2+g7iy(i,j)^2);
end
end  
k=1;
for i=2:nx-1
for j=2:ny-1
        g9(k)=g1(i,j);
        g71(k)=g7x(i,j);
        g72(k)=g7y(i,j);
        k=k+1;
end
end
 z=1;
 X1=zeros(N,1);
 Y1=zeros(N,1);
 G1=zeros(N,1);
for j=2:ny-1
for i=2:nx-1
    n=0;
    v1=1;
        xmax=zeros(4,1);
        ymax=zeros(4,1);
        gmax=zeros(4,1);
if g1(i,j)>((g1(i-1,j)+g1(i+1,j))/(2))   
            d=pasx;
            a=(g1(i-1,j)-2*g1(i,j)+g1(i+1,j))/((2)*(d^2));
            b=(g1(i+1,j)-g1(i-1,j))/((2)*(d));    
            if (abs((-1)*(b)/(2*a)))<=(d)
            x=(-1)*(b)/(2*a);
            if x>((g7x(i,j)-g7x(i-1,j)))/(2)&&x>((g7x(i+1,j)-g7x(i,j)))/(2)
            if x>((g7y(i,j)-g7y(i-1,j)))/(2)&&x>((g7y(i+1,j)-g7y(i,j)))/(2)
            gmax(v1)=(a)*(x^2)+(b)*x+g1(i,j);
            xmax(v1)=pasx*(i-1)+x+ax;
            ymax(v1)=pasy*(j-1)+ay;
            v1=v1+1;
            n=n+1;
            end
            end
            end
end
if  g1(i,j)>((g1(i,j-1)+g1(i,j+1))/(2))
            d=pasy;
            a=(g1(i,j-1)-2*g1(i,j)+g1(i,j+1))/((2)*(d^2));
            b=(g1(i,j+1)-g1(i,j-1))/((2)*(d));
            if (abs((-1)*(b)/(2*a)))<=(d)
            y=(-1)*(b)/(2*a);
            if y>((g7x(i,j)-g7x(i,j-1)))/(2)&&y>((g7x(i,j+1)-g7x(i,j)))/(2)
            if y>((g7y(i,j)-g7y(i,j-1)))/(2)&&y>((g7y(i,j+1)-g7y(i,j)))/(2)
            gmax(v1)=(a)*(y^2)+(b)*y+g1(i,j);
            xmax(v1)=pasx*(i-1)+ax;
            ymax(v1)=pasy*(j-1)+y+ay;
            v1=v1+1;
            n=n+1;
            end
            end
            end
end
if g1(i,j)>((g1(i-1,j-1)+g1(i+1,j+1))/(2))
            d=sqrt(pasy^2+pasx^2);
            a=(g1(i-1,j-1)-2*g1(i,j)+ g1(i+1,j+1))/((2)*(d^2));
            b=(g1(i+1,j+1)-g1(i-1,j-1))/((2)*(d));
            if (abs((-1)*(b)/(2*a)))<=(d)
            x=(-1)*(b)/(2*a);
            if x>((g7x(i,j)-g7x(i-1,j-1)))/(2)&&x>((g7x(i+1,j+1)-g7x(i,j)))/(2)
            if x>((g7y(i,j)-g7y(i-1,j-1)))/(2)&&x>((g7y(i+1,j+1)-g7y(i,j)))/(2)
            gmax(v1)=(a)*(x^2)+(b)*x+g1(i,j);
            xmax(v1)=pasx*(i-1)+x/sqrt(2)+ax;
            ymax(v1)=pasy*(j-1)+x/sqrt(2)+ay;
            v1=v1+1;
            n=n+1;
            end
            end
            end
end
if g1(i,j)>((g1(i-1,j+1)+g1(i+1,j-1))/(2))
            d=sqrt(pasy^2+pasx^2);
            a=(g1(i-1,j+1)-2*g1(i,j)+g1(i+1,j-1))/((2)*(d^2));
            b=(g1(i+1,j-1)-g1(i-1,j+1))/((2)*(d));            
            if (abs((-1)*(b)/(2*a)))<=(d)
            x=(-1)*(b)/(2*a);
            if x>((g7x(i,j)-g7x(i-1,j+1)))/(2)&&x>((g7x(i+1,j-1)-g7x(i,j)))/(2)
            if x>((g7y(i,j)-g7y(i-1,j+1)))/(2)&&x>((g7y(i+1,j-1)-g7y(i,j)))/(2)       
            gmax(v1)=(a)*(x^2)+(b)*x+g1(i,j);
            xmax(v1)=pasx*(i-1)+x/sqrt(2)+ax;
            ymax(v1)=pasy*(j-1)+x/sqrt(2)+ay;
            v1=v1+1;
            n=n+1;
            end
            end
            end
end
if n==1 
    q=(2*a)/max(g1(i,j));
    if q<c 
     for t=1:n
        G1(z)=gmax(t);
        X1(z)=xmax(t);
        Y1(z)=ymax(t);
        z=z+1;
     end
    end
end
if n>1
     for t=1:n
        G1(z)=gmax(t);
        X1(z)=xmax(t);
         Y1(z)=ymax(t);
        z=z+1;

     end

end
end
end
B=zeros(N,5);
for i=1:N
B(i,1)=g6(i,1);
B(i,2)=g9(i,1);
B(i,3)=X1(i,1);
B(i,4)=Y1(i,1);
B(i,5)=G1(i,1);
end
xlswrite(exceldocument,B,e+1,'A2:E5000');
xlswrite(exceldocument,head2,e+1,'A1:E1');
xl=xl+pas;
end 

