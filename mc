format long;
% Define the file path and file properties

file_path = '/Users/jenniferlundberg/Desktop/Monte-carlo/RFA341/FullPhantom.bin';
x = 256; % Replace with your desired dimensions
y = 256;
z = 1200;
voxelSize = 0.15; % cm
% Open the binary file for reading
fid = fopen(file_path, 'rb');
if fid == -1
    error('Could not open the file.');
end
% Read the binary data
raw_data = fread(fid, x * y * z, 'single', 0, 'ieee-le');
% Close the file
fclose(fid);
% Reshape the raw data into a 3D array
array_3d = reshape(raw_data, [x, y, z]);



A = table2array(readtable("Attenueringsdata.xlsx"));  % Läser in tabelldata
A(:,1) = A(:,1)*1000; 
energivarden_attenuering = A(:,1);
kidney_varden = A(:,10); %Attenuering i kidney
opts = detectImportOptions("SOFT_TISSUE_TVARSNITT.xlsx");
opts.PreserveVariableNames = true;
T =table2array(readtable("SOFT_TISSUE_TVARSNITT.xlsx",opts));
energivarden_tvarsnitt = T(:,1); 
rayleigh_varden = T(:,2); %Värden på de olika tvärsnitten
compton_varden = T(:,3);
fotoel_varden = T(:,4);
fotoner = 10;  % Antal fotoner

% Koordinat-vektorer
x = zeros(fotoner);  % X-koordinat
y = zeros(fotoner);  % Y-koordinat
z = zeros(fotoner);  % Z-koordinat

% Simulerings-loop för varje foton
for foton = 1:fotoner
    theta_old = acos(2*rand-1); % Polarvinkel
    fi_old = rand*2*pi; % Azimuthalvinkel
    
    [fotonenergi] = bestamma_fotonenergi(); %Funktion som bestämmer fotonens energi
 
 steg = 1; %Visar vilket steg den är på
 % Alla fotoner börjar i origo
 x(foton, 1) = 0;
 y(foton, 1) = 0;
 z(foton, 1) = 0;
 
 % Första riktningen
 u_0 = sin(theta_old)*cos(fi_old);
 v_0 = sin(theta_old)*sin(fi_old);
 w_0 = cos(theta_old);
 
 %Fotonerna tar en slumpmässig riktning i början
 [theta, phi] = vinkel();
 [attenuering_kidney] = attenuering(fotonenergi, AA, energivarden_attenuering, kidney_varden);
 stegl = -log(rand())/attenuering_kidney;
 
 
 % Första steget
 x_0 = stegl*u_0;
 y_0 = stegl*v_0;
 z_0 = stegl*w_0;
 
 % Fotonen flyttas
 x(foton, steg+1) = x_0;
 y(foton, steg+1) = y_0;
 z(foton, steg+1) = z_0;
 Rotation = [cos(theta)*cos(phi) -sin(phi) sin(theta)*cos(phi);
            cos(theta)*sin(phi) cos(phi) sin(theta)*sin(phi);
            -sin(theta) 0 cos(theta)];

 
 steg = steg + 1;
 while fotonenergi > 50 %Sätter gränsen på 50 eV. Alltså försvinner fotonen när den har mindre energi än så
     [tvarsnitt_Compton, tvarsnitt_Rayleigh, tvarsnitt_fotoel] = tvarsnitt(fotonenergi, T, energivarden_tvarsnitt, rayleigh_varden, compton_varden, fotoel_varden);
     summa = tvarsnitt_Compton + tvarsnitt_Rayleigh + tvarsnitt_fotoel;
     R = rand();
     if R <= tvarsnitt_fotoel/summa
         %fotoelektrisk effekt
         fotonenergi = 0;
         
     elseif R <= (tvarsnitt_fotoel + tvarsnitt_Compton)/summa
         %Compton
         [fotonenergi, costheta] = Compton(fotonenergi);
         theta_new = costheta;
         fi_new = rand*2*pi;
         [phi, theta] = vinkel();
         
         %Fotonen flyttas
         [attenuering_kidney] = attenuering(fotonenergi, AA, energivarden_attenuering, kidney_varden);
         stegl = -log(rand())/attenuering_kidney;
         
         % Ekvation 7.13 från Bielajew.
         u_1 = u_0*cos(theta_new)+sin(theta_new)*(w_0*cos(fi_new)*cos(fi_old)-sin(fi_new)*sin(fi_old));
         v_1 = v_0*cos(theta_new)+sin(theta_new)*(w_0*cos(fi_new)*sin(fi_old)+sin(fi_new)*cos(fi_old));
         w_1 = w_0*cos(theta_new)-sin(theta_new)*sin(theta_old)*cos(fi_new);
         
         % Nästa position
         x_1 = x_0 + stegl*u_1;
         y_1 = y_0 + stegl*v_1;
         z_1 = z_0 + stegl*w_1;
         
         % Skriv över gammal vinkel och riktningar
         u_0 = u_1;
         v_0 = v_1;
         w_0 = w_1;
         
         x_0 = x_1;
         y_0 = y_1;
         z_0 = z_1;
         
         theta_old = theta_new;
         fi_old = fi_new;

         
         x(foton, steg+1) = x_0;
         y(foton, steg+1) = y_0;
         z(foton, steg+1) = z_0;
         
         steg = steg + 1;
     else
         %Rayleigh
         steg = steg + 1;
     end
     
 end

end
% Plotta fotonernas rörelse i 3D
figure;
for foton = 1:fotoner
    plot3(x(foton, 1:steg), y(foton, 1:steg), z(foton, 1:steg), 'b.-');
    hold on;
end

% Lägger en röd prick i origo
scatter3(0, 0, 0, 25, 'r', 'filled'); % storlek = 25

hold off;


axis equal;
title('Punktkälla med fotoner (3D)');
xlabel('X-koordinat');
ylabel('Y-koordinat');
zlabel('Z-koordinat');


function [fotonenergi] = bestamma_fotonenergi()
energibestamning = rand();
if 0 < energibestamning & energibestamning < 0.00854;
    fotonenergi = 71600;  %Olika sannolikheter för olika fotonenergier enl. sönderfallsdata
elseif 0.00854 < energibestamning & energibestamning < 0.36354;
    fotonenergi = 113000;
elseif 0.36354 < energibestamning & energibestamning < 0.3662;
    fotonenergi = 137000;
elseif 0.3662 < energibestamning & energibestamning < 0.9762;
    fotonenergi = 208000;
elseif 0.9762 < energibestamning & energibestamning < 0.995;
    fotonenergi = 250000;
elseif 0.995 < energibestamning & energibestamning < 1;
    fotonenergi = 321000;
end
end
function [theta, phi] = vinkel()
     theta = acos(1-2*rand());  % Cosinus av polarvinkeln mellan -1 och 1
     phi = 2 * pi * rand();       % asimutalvinkeln mellan 0 och 2pi
end
%Funktion som linjärinterpolerar tvärsnitten
function [tvarsnitt_Rayleigh, tvarsnitt_Compton, tvarsnitt_fotoel] = tvarsnitt(fotonenergi, T, energivarden_tvarsnitt, rayleigh_varden, compton_varden, fotoel_varden)
    % Hittar index på de närmsta energierna
    [~, index1] = min(abs(T(:,1) - fotonenergi));
    index2 = index1 + 1;
    % Linjärinterpolation
    x1 = energivarden_tvarsnitt(index1);
    x2 = energivarden_tvarsnitt(index2);
    r1 = rayleigh_varden(index1);  %Värdena på Rayleigh
    r2 = rayleigh_varden(index2);
    c1 = compton_varden(index1); %Värden på Compton
    c2 = compton_varden(index2);
    f1 = fotoel_varden(index1); %Värden på fotoelektrisk effekt
    f2 = fotoel_varden(index2);
    
    tvarsnitt_Rayleigh = r1 + (fotonenergi - x1) * (r2 - r1) / (x2 - x1);
    tvarsnitt_Compton = c1 + (fotonenergi - x1) * (c2 - c1) / (x2 - x1);
    tvarsnitt_fotoel = f1 + (fotonenergi - x1) * (f2 - f1) / (x2 - x1);
end
%Funktion som interpolerar fram attenuering i njuren
function [attenuering_kidney] = attenuering(fotonenergi, AA, energivarden_attenuering, kidney_varden)
     % Hittar index på de närmsta energierna
    [~, index1] = min(abs(AA(:,1) - fotonenergi));
    index2 = index1 + 1;
        % Linjärinterpolation
    x1 = energivarden_attenuering(index1);
    x2 = energivarden_attenuering(index2);
    k1 = kidney_varden(index1);  %Värdena på Kidney
    k2 = kidney_varden(index2);
    attenuering_kidney = k1 + (fotonenergi - x1) * (k2 - k1) / (x2 - x1);
end
%Funktion som ger vinkeln
function [costheta,fotonenergi] = Compton(fotonenergi)
    alpha = fotonenergi/511; %alpha är fotonenergin i elektronmassor
    %Lista från Kahn
    b = 1;
    while b == 1 %Gör om dessa tills något krav är uppfyllt
        R1 = rand();
        R2 = rand();
        R3 = rand();
        if R1 <= (2*alpha+1)/(2*alpha+9)
            eta = 2*alpha*R2;
            
            if R3 <= 4*(eta^(-1)-eta^(-2))
                costheta = 1-2*R2;
                b = 0;
            end
        elseif R1 > (2*alpha+1)/(2*alpha+9)
            eta = (2*alpha+1)/(2*R2*alpha+1);
            costheta = 1-(eta-1)/alpha;
            
            if R3 <= 0.5*(costheta^2+eta^(-1))
                b = 0;
            end
        end
        
    end
    fotonenergi = fotonenergi/eta;
end
