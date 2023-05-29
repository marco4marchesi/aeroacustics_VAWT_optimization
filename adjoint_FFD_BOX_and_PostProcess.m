%{

FFD_BOX computation, plot and write config lines

set the parameters for the FFD box and it will print what you need to paste
in the config file

example
FFD_DEFINITION = (airfoil, -0.02, 0.74, 0.0, 0.0575, 0.74, 0.0, 0.0575, 0.76, 0.0, -0.02, 0.76, 0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0 )

%}


clear; close all; clc; 
folder_varie = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\varie\";
simulationsAdjointPath = "\\wsl.localhost\Ubuntu-20.04\home\marco\testAdjoint_momentZ_30pti2_deform6\DESIGNS\";
list = dir(simulationsAdjointPath+"DSN*");
% simulationsDeformPath = simulationsAdjointPath + "DSN_" + compose("%03d",length(list)-1)+ "\DEFORM\";
simulationsDeformPath = simulationsAdjointPath + "DSN_" + compose("%03d",15)+ "\DEFORM\";
simulationsBaseGeoPath = erase(simulationsAdjointPath,"DESIGNS\");

matlab_graphics;

user = 'marco';
user_settings;

%% user defined properties
% extract profile points
profile = readmatrix(folder_varie+"NACA0021.txt");
profile = profile(2:801, :);

% number of intermediate points per dimension
N_x = 30;
N_y = 2; 
N_z = 1;

flag_save = false;

%% box generation script
xleft = min(profile(:,1));
xright = max(profile(:,1));
xlen = xright-xleft;

ybott = min(profile(:,2));
yup = max(profile(:,2));
ylen = yup - ybott;

x_offRatio = 0.02; %we want to stay 2% of the chord distant
x_offset = x_offRatio*xlen;
xleft_box = xleft-x_offset;
xleft_box = round(xleft_box*100)/100;
x_offset = xleft-xleft_box;
xright_box = xright + x_offset;

y_offRatio = 0.02; %we want to stay 2% of the chord distant
y_offset = y_offRatio*ylen;
ybott_box = ybott-y_offset;
ybott_box = round(ybott_box*100)/100;
y_offset = ybott-ybott_box;
yup_box = yup + y_offset;

% define the points of the FFD box
P1 = [xleft_box, ybott_box, 0];
P2 = [xright_box, ybott_box, 0];
P3 = [xright_box, yup_box, 0];
P4 = [xleft_box, yup_box, 0];
boxPoints = [P1;P2;P3;P4];
boxPoints_T = boxPoints';


% define plot lines ( this is good until N_y becomes larger than 2 )
lines_bott = [linspace(xleft_box,xright_box,N_x)',ones(N_x,1)*ybott_box];
lines_up = [linspace(xleft_box,xright_box,N_x)',ones(N_x,1)*yup_box];

% define the markers
marker_computation = 'AIRFOIL1';
marker_FFD = 'airfoil';

% generate the .txt file
fid = fopen('paste_to_config.txt','wt');
if fid ~= -1
    %ffd box definition
    fprintf(fid,'FFD_DEFINITION= (%s, ',marker_FFD);
    if size(boxPoints,1)*size(boxPoints,2) == 12
          
        for i = 1:size(boxPoints,1)*size(boxPoints,2)*2
            if i <=12
                fprintf(fid,'%.4f, ',boxPoints_T(i)');
            elseif i> 12 && i<24
                fprintf(fid,'0.0, ');
            else
                fprintf(fid,'0.0) ');
            end
        end
    else
        for i = 1:size(boxPoints,1)*size(boxPoints,2)
            if i <24
                fprintf(fid,'%.4f, ',boxPoints(i));
            else 
                fprintf(fid,') ');
            end
        end
    
    end
    fprintf(fid,')\n');

    % ffd box degrees of freedom
    fprintf(fid,'FFD_DEGREE= (%d, %d, %d)\n',N_x-1, N_y-1, N_z-1);
    
    % definition dv variables
    fprintf(fid,'DEFINITION_DV= ');
    for i = 1:N_x
        for j = 1:N_y
            fprintf(fid,'(11, 1.0 | %s | %s,  %d, %d, 0, 0.0, 1.0, 0.0 )',marker_computation, marker_FFD, i-1, j-1);
            if i*j < N_x*N_y
                fprintf(fid,';');
            end
        end
    end
    fprintf(fid,'\n');
    
    % close the file
    fclose(fid);
else
    warningMessage = sprintf('Cannot open file %s', filename);
    uiwait(warndlg(warningMessage));
end

%% PLOT THE BOX
fig_boxAdjoint = figure('Position',[100,100,800,500]);
scatter(profile(:,1),profile(:,2))
hold on; 
plot([boxPoints(:,1);boxPoints(1,1)],[boxPoints(:,2);boxPoints(1,2)],'r')
for i = 1:length(lines_bott)
    plot([lines_bott(i,1),lines_up(i,1)],[lines_bott(i,2),lines_up(i,2)],'r--','LineWidth',1)
end
axis equal
xlim([-0.03, 0.07])
if flag_save
    exportgraphics(fig_boxAdjoint,folder_varie+"adjoint Box.png")
end

%% retrieve deformed profile
fid= fopen(simulationsDeformPath+"surface_deformed.dat");
i = 0;
while i<1000
        i = i+1;
    a = fgets(fid);
    if a == -1
        continue
    end
    array{i} = a;
end
fclose(fid)

array = array(1,4:end/2);
% for now works only in 2D
pointsDeformed = zeros(length(array),2);
for i = 1:length(array)
    x = str2double(array{i}(1:end/2));
    y = str2double(array{i}(end/2:end));
    pointsDeformed(i,:) = [x,y];
end

fid= fopen(simulationsBaseGeoPath+"surface_deformed.dat");
i = 0;
while i<1000
        i = i+1;
    a = fgets(fid);
    if a == -1
        continue
    end
    array{i} = a;
end
fclose(fid);

array = array(1,4:end/2);
% for now works only in 2D
pointsBaseProfile = zeros(length(array),2);
for i = 1:length(array)
    x = str2double(array{i}(1:end/2));
    y = str2double(array{i}(end/2:end));
    pointsBaseProfile(i,:) = [x,y];
end


% sort points
mat = [pointsBaseProfile,pointsDeformed];
mat(:,2) = mat(:,2)-0.75;
[~,idx] = sort(mat(:,2));
mat = mat(idx,:);
l = size(mat,1);
mat_l = mat(1:l/2,:);
mat_u = mat(l/2+1:end,:);
[~,idx] = sort(mat_l(:,1));
mat_l = flip(mat_l(idx,:));
[~,idx] = sort(mat_u(:,1));
mat_u = mat_u(idx,:);
mat = [mat_l;mat_u];
mat(:,2) = mat(:,2)+0.75;

% recall points
pointsDeformed = mat(:,3:4)-[0,0.6];
pointsBaseProfile = mat(:,1:2)-[0,0.6];

% retrieve coefficients
history_adj = readmatrix(simulationsDeformPath+"../DIRECT/history_direct.dat");
CMz_adj = history_adj(end,16);
history_start = readmatrix(simulationsAdjointPath+"DSN_001/DIRECT/history_direct.dat");
CMz_start = history_start(end,16);

disp("CMz start = "+num2str(CMz_start))
disp("CMz adjoint = "+num2str(CMz_adj))

figure
plot(pointsBaseProfile(:,1),pointsBaseProfile(:,2),"DisplayName","Initial")
hold on;
plot(pointsDeformed(:,1),pointsDeformed(:,2),"DisplayName","Adjoint result")
axis equal
for i = 1:length(lines_bott)
    plot([lines_bott(i,1),lines_up(i,1)],[lines_bott(i,2),lines_up(i,2)],'r--','LineWidth',1,'HandleVisibility','off')
end
legend

%% generate 3 profile mesh
N_profili = 3;
theta_vec = linspace(0,-2*pi,N_profili+1); % senso orario per consistenza
Rot = @(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
for i = 1:N_profili
    profili{i} = [Rot(theta_vec(i))*[pointsDeformed, zeros(size(pointsDeformed,1),1)]']';
    profili_nom{i} = [Rot(theta_vec(i))*[pointsBaseProfile, zeros(size(pointsBaseProfile,1),1)]']';
end

figure
hold on; grid on;
for idx = 1:N_profili
    plot(profili{idx}(:,1),profili{idx}(:,2),'k',"DisplayName",num2str(idx))
end
axis equal 

fig_posizione3profili = figure();
hold on; 
for idx = 1:N_profili
    plot(profili_nom{idx}(:,1),profili_nom{idx}(:,2),'Color',[0,0,128]./256,"DisplayName",num2str(idx))
end
scatter(0,0)
axis equal
axis off
exportgraphics(fig_posizione3profili, '01_posizione3profili.emf')

fid = fopen([num2str(N_profili),'_profile_generator.txt'],'w');
for j = 1:length(profili)
    for i = 1:size(profili{1},1)
        fprintf(fid,'Point(%d) = {%.16f, %.16f, 0.0, h};\n',(j-1)*1000+i,profili{j}(i,1),profili{j}(i,2));
    end
end
fclose(fid);

