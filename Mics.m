function file = Mics(R,x,y,n)
Rx=R*x;
Ry=R*y;
R=10*R;
x0= Rx/2;
y0= Ry/2;
R2=R+(Rx/2);
ang1=linspace(-pi/2,pi/2,n);
%arco davanti propeller1
x1=cos(ang1).*R;
y1= sin(ang1).*R;
%arco grande
ang2=ang1-(pi/2);
x2=cos(ang2).*R2;
x2=x2-x0;
y2=sin(ang2).*R2;
y2=y2+y0;
%arco dietro
ang3=ang1+pi;
x3=cos(ang3).*R;
x3=x3-Rx;
y3=sin(ang3).*R;
y3=y3+Ry;

file = fopen('Observer_Locations.dat','a+');
fprintf(file,'%i\n',n*3);
for i=1:n
            fprintf(file,'%.5f %.5f 0.0\n',x1(i)-0.19,y1(i));

end
for i=1:n
            fprintf(file,'%.5f %.5f 0.0\n',x2(i)-0.19,y2(i));

end
for i=1:n
            fprintf(file,'%.5f %.5f 0.0\n',x3(i)-0.19,y3(i));

end
fprintf(file,'\r\n\r\n');
fclose(file);
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
LW = 2.5;
FS = 30;
figure;
plot(x1,y1,'-*','Linewidth',LW,'MarkerSize',8);
hold on
grid on
grid minor
plot(x2,y2,'-*','Linewidth',LW,'MarkerSize',8);
plot(x3,y3,'-*','Linewidth',LW,'MarkerSize',8);
ylim([-2.2 1.8])
xlim([-3 2])
legend("External Arch", "Front Arch", "Rear Arch");
title("Mics Positions");
xlabel('x')
ylabel("y")

end

