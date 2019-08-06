% Initial Ar gas velocity according to the Maxwell distribution (negative z direction)
close all
clear
% HI I MADE this CHANGE
kB=1.38064852*10^(-23); %[SI unit]
mg=39.948/1000/(6.02*10^23); %Ar molecular mass [kg]
Tg=295; %[K]

addmol=300; %total number of molecules want to be added

vx=zeros(1,addmol); %velocity components
vy=zeros(1,addmol);
vz=zeros(1,addmol);

%------sampling--------%
for i=1:addmol
rnd_num = rand(1,9);
beta = sqrt(mg/(2*Tg*kB));
vx(i) = sqrt(-log(rnd_num(:,1)))./beta.*sin(2*pi.*rnd_num(:,2));
vy(i) = sqrt(-log(rnd_num(:,3)))./beta.*sin(2*pi.*rnd_num(:,4));
vz(i) = -sqrt(-log(rnd_num(:,5)))./beta;
end

vx=vx*10^10; %convert velocity magnitudes [Angstrom/s]
vy=vy*10^10;
vz=vz*10^10;

%-----verification--------%
figure(1)
v_x=-1000:1:1000;
f_vx=(mg/(2*pi*kB*Tg))^0.5.*exp(-mg*v_x.^2./(2*kB*Tg)); %theoretical distribution for tangential component (vx)
histogram(vx./(10^10),25,'Normalization','pdf') %actual distribution (vx)
hold on
plot(v_x,f_vx)
hold off
xlabel('v_x (m/s)');
ylabel('f(v_x)')
legend('Sampled velocities','Maxwell distribution','Location','Best');
saveas(gcf,'Velocity distribution for tangential.fig')

figure(2)
v_z=-1000:1:0;
sigma=sqrt(kB*Tg/mg);
f_vz=-v_z.*exp(-v_z.^2./(2*sigma^2))./(sigma^2); %theoretical distribution for normal component (vz)
histogram(vz./(10^10),25,'Normalization','pdf') %actual distribution (vz)
hold on
plot(v_z,f_vz)
hold off
xlabel('v_z (m/s)');
ylabel('f(v_z)')
legend('Sampled velocities','Rayleigh distribution','Location','Best');
saveas(gcf,'Velocity distribution for normal.fig')

%-----saving velocities------%
v=[vx;vy;vz];

%save the velocities according to ADF format requirements
fid=fopen('addmol.vel','w+');
fprintf(fid,'v[%d]="%24.17E%24.16E%24.16E"\r\n',[1:1:addmol;v]);
fclose(fid);
