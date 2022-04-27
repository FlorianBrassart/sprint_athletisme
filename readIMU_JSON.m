%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% filename = chemin/nom du fichier.bin
% frhz = fréquence d'acquisition
% rayon_roue = pour le calcul de la vitesse en km/h indiquer
% le rayon de la roue en m
% OUPUT
% DATA
% Colonne 1 : temps en sec
% Colonne 2 : accélération en X
% Colonne 3 : accélération en Y
% Colonne 4 : accélération en Z
% Colonne 5 : gyromètre en X
% Colonne 6 : gyromètre en Y
% Colonne 7 : gyromètre en Z
% Colonne 8 : magnétomètre en X
% Colonne 9 : magnétomètre en Y
% Colonne 10 : magnétomètre en Z
% Colonne 11 : vitesse en km/h
% Colonne 12 : accélération linéaire

% Exemple
% filename = 'C:\Users\Kirespi\Documents\MATLAB\fonction ATN capteurs\S2_E6D1056849DF.bin'
% frhz = 128
% rayon_roue = 0.3
% DATA = readIMU_JSON(filename,frhz,rayon_roue)
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DATA = readIMU_JSON(filename,frhz,rayon_roue)

clear fid data_bin data_bin_brute

%--- Import et lecture des fichier .bin
DATA_struct = jsondecode(fileread(filename));

%--- Stockage des datas dans Matrice

% accélero
DATA(:,2) = DATA_struct.acc_x(10:end-10)*(32/2^16);
DATA(:,3) = DATA_struct.acc_y(10:end-10)*(32/2^16);
DATA(:,4) = DATA_struct.acc_z(10:end-10)*(32/2^16);

% gyro
DATA(:,5) = DATA_struct.gyro_x(10:end-10)*(4000/2^16);
DATA(:,6) = DATA_struct.gyro_y(10:end-10)*(4000/2^16);
DATA(:,7) = DATA_struct.gyro_z(10:end-10)*(4000/2^16);

% magneto
DATA(:,8) = DATA_struct.mag_x(10:end-10);
DATA(:,9) = DATA_struct.mag_y(10:end-10);
DATA(:,10) = DATA_struct.mag_z(10:end-10);

% tps
time_unit = 1/frhz;
time = time_unit : time_unit: time_unit*(size(DATA,1));time = time';
DATA(:,1) = time;
            

if ~isempty(rayon_roue)
    % calcul de la vitesse de déplacement du fauteuil
    clear vitesse
    vitesse = rayon_roue*(deg2rad(DATA(:,7)))*3.6;
    DATA(:,11) = vitesse;
    
    % calcul accélération linéaire à partir de la vitesse
    clear Acc_frm vit_frm
    vit_frm(1,1) = DATA(1,11)-(DATA(2,11)-DATA(1,11));
    vit_frm(2:size(DATA)+1,1) = DATA(1:end,11);
    vit_frm(end+1,1) = DATA(end,11)+(DATA(end,11)-DATA(end-1,11));
    for i = 2:size(vit_frm,1)-1
        Acc_frm(i-1,1) = ((vit_frm(i+1,1)/3.6)-(vit_frm(i-1,1)/3.6))/(2/frhz);
    end
    DATA(:,12) = Acc_frm;  
end

end

