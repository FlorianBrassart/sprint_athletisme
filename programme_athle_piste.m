%% fonction appelée : 
% readIMU_JSON.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc

rr=dir([cd '\Result_athle_piste.mat']);
if ~isempty(rr)
load('Result_athle_piste.mat')
end

%% Selection du dossier 
path=uigetdir(path,'Dossier contenant les fichiers .json'); % permet d'ouvrir une pop up pour choisir le dossier dans l'ordi
dirData_J = dir([path '/*.json']);      %# Get the data for the current directory
fileList = {dirData_J.name}';  %'# Get a list of the files
split_folder=split(path,'\');
nom=split_folder(end);

num_file=5;
%% ----- ouverture des fichers
% IMU
%% Ouverture des fichiers 
FrHZ=128; % fréquence d'acquisition des IMUs
M=80; % masse du système choisie arbitrairement pour l'instant
cc1=contains(fileList,'Droite');
cc2=contains(fileList,'Gauche');
cc3=contains(fileList,'Cale');
file1=strcat(path,'\',char(fileList(cc1)));
file2=strcat(path,'\',char(fileList(cc2)));
file3=strcat(path,'\',char(fileList(cc3)));
if length(fileList)==num_file
    cc4=contains(fileList,'Torse');
    file4=strcat(path,'\',char(fileList(cc4)));
else
    file4=strcat(path,'\',cell2mat(fileList(3,1))); % mise en lien du chemin d'accés + le nom du fichier 2
end
cc5=contains(fileList,'information');
file_infos=strcat(path,'\',char(fileList(cc5)));
INFOS=jsondecode(fileread(file_infos)); % inforamtions du fauteuil utilisées dans les calculs
rayon_roue=(str2num(INFOS.taille_roues)/39.37)/2;% transformation du diamètre en pouce renseigné par le rayon en metre. 
camber_angle=str2num(INFOS.angle_carossage);
% dist_roue=str2num(INFOS.largeur_fauteuil_sol);

DATAR1_struct = jsondecode(fileread(file1));
DATAR2_struct = jsondecode(fileread(file2));
DATAF_av_struct = jsondecode(fileread(file3));
DATADos_struct = jsondecode(fileread(file4));
% fonction permettant de récupérer les données des fichiers .json sous forme de matrice en argument il faut renseigner (chemain d'accès + nom du fichier, la fréquence à laquelle les données ont étaient enregistrées.
DATAR1 = readIMU_JSON(file1,FrHZ,rayon_roue); % Roue droire
DATAR2 = readIMU_JSON(file2,FrHZ,rayon_roue); % Roue Gauche
DATAF_av= readIMU_JSON(file3,FrHZ,rayon_roue); % IMU placée sur le cadre
DATADos= readIMU_JSON(file4,FrHZ,rayon_roue); % IMU placée sur le dos
%%
% DATAR2 = resample (DATAR2, FrHZ, 256);
% DATAR1 = resample (DATAR1, FrHZ, 256);
%%

if mean(DATAR1(:,11))<0
DATAR1(:,2:end)=-DATAR1(:,2:end); % lors des tests le sujet doit aller majoritairement vers l'avant, on s'attend donc à une vitesse de rotation des roues majoritairement positive. 
end
if mean(DATAR2(:,11))<0
DATAR2(:,2:end)=-DATAR2(:,2:end); % si jamais la droite et la gauche on été inversé lors du test.
disp('iversion roues')
end
clear lgt zz
lgt(1,1)=length(DATAR1); % si jamais les matrices ne fond pas toute la même taille (perte de connexion d'une IMU au cours du test)
lgt(2,1)=length(DATAR2);
lgt(3,1)=length(DATAF_av);
lgt(4,1)=length(DATADos);
[zz, ~]=find(lgt==min(lgt));
zz=zz(1,1);
clear DATA
DATA= DATAR1(1:lgt(zz,1),:);
clear DATAR1
DATAR1= DATA;
clear DATA
DATA= DATAR2(1:lgt(zz,1),:);
clear DATAR2
DATAR2 = DATA;
clear DATA
DATA= DATAF_av(1:lgt(zz,1),:);
clear DATAF_av
DATAF_av=DATA;
clear DATA
DATA= DATADos(1:lgt(zz,1),:);
clear DATADos
DATADos=DATA;
clear DATA
sizeDATA=size(DATADos);
% accélération de G en m/s²
DATAF_av(:,2)=DATAF_av(:,2)-(mean(DATAF_av(50:150,2)));
DATAF_av(:,3)=DATAF_av(:,3)-(mean(DATAF_av(50:150,3)));
DATAF_av(:,4)=DATAF_av(:,4)-(mean(DATAF_av(50:150,4)));
DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81;
%% Calcul de la vitesse en prenant en compte l'angle de carrossage d'après fuss et al 2012
dif=DATAR2(:,7)-DATAR1(:,7);
dif(dif==0)=0.0001;                              
signe=dif(:,1)./abs(dif(:,1));
Mwz= deg2rad(DATAR1(:,7));
Mwxy =deg2rad(sqrt((DATAR1(:,5)).^2+(DATAR1(:,6)).^2));
tan0=deg2rad(tan(camber_angle));
Mwz1= deg2rad(DATAR2(:,7));
Mwxy1=deg2rad(sqrt((DATAR2(:,5)).^2+(DATAR2(:,6)).^2));
DATAR1(:,13)=(Mwz-(Mwxy.*tan0.*signe)).*rayon_roue.*3.6;
DATAR2(:,13)=(Mwz1-(Mwxy1.*tan0.*(-signe))).*rayon_roue.*3.6;                                        
%% Filtre DATA IMU
Wn=(10/FrHZ);
[af,bf] = butter(2,Wn,'low');
DATAR1(:,2:end) = filtfilt(af,bf,DATAR1(:,2:end));
DATAR2(:,2:end) = filtfilt(af,bf,DATAR2(:,2:end));
DATADos(:,2:end) = filtfilt(af,bf,DATADos(:,2:end));
DATAF_av(:,2:end) = filtfilt(af,bf,DATAF_av(:,2:end));
%% Calcul de la vitesse et de l'accélération moyenne
DATAmean1=(DATAR1(:,13)+DATAR2(:,13))/2;
DATAmean2=DATAmean1(1,1)-(DATAmean1(2,1)-DATAmean1(1,1));
DATAmean2=[DATAmean2 ;DATAmean1];
Accmean = diff(DATAmean2/3.6)*FrHZ; clear DATAmean2

DATAR1_2=DATAR1(1,13)-(DATAR1(2,13)-DATAR1(1,13));
DATAR1_2=[DATAR1_2 ;DATAR1(:,13)];
AccR1 = diff(DATAR1_2/3.6)*FrHZ; clear DATAR1_2
DATAR2_2=DATAR2(1,13)-(DATAR2(2,13)-DATAR2(1,13));
DATAR2_2=[DATAR2_2 ;DATAR2(:,13)];
AccR2 = diff(DATAR2_2/3.6)*FrHZ; clear DATAR2_2

meandist_tot=cumtrapz(DATAmean1/3.6)/FrHZ;
meandist_R1=cumtrapz(DATAR1(:,13)/3.6)/FrHZ;
meandist_R2=cumtrapz(DATAR2(:,13)/3.6)/FrHZ;
Force=Accmean*M;
Force_eff=Force;
Force_eff(Force_eff<0)=0;
Puissance=Force_eff.*(DATAmean1/3.6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul d'orientation avec les quaternions 
cadre_eul=rad2deg(quat2eul([DATAF_av_struct.quaternion_w DATAF_av_struct.quaternion_x DATAF_av_struct.quaternion_y DATAF_av_struct.quaternion_z],'XYZ'));
cadre_eul=unwrap(cadre_eul);
%%
% DATADos_struct.quaternion_w   = resample (DATADos_struct.quaternion_w, FrHZ, 256);
% DATADos_struct.quaternion_x   = resample (DATADos_struct.quaternion_x, FrHZ, 256);
% DATADos_struct.quaternion_y   = resample (DATADos_struct.quaternion_y, FrHZ, 256);
% DATADos_struct.quaternion_z   = resample (DATADos_struct.quaternion_z, FrHZ, 256);
%%
Dos_eul=rad2deg(quat2eul([DATADos_struct.quaternion_w DATADos_struct.quaternion_x DATADos_struct.quaternion_y DATADos_struct.quaternion_z],'XYZ'));
Dos_eul=-(unwrap(Dos_eul));
% Calcul d'orientation avec la vitesse de rotation du gyrometre
%% Calcule de l'orientation du fauteil

orientationDos=[0 0 0];
orientationDos1=[0 0 0];
orientationCadre=[0 0 0 0 0];
orientationCadre1=[0 0 0 0 0];
DATA(:,1)=DATADos(:,1);
DATA(:,2:4)=DATADos(:,5:7);
DATA(:,5)= DATAF_av(:,1);
DATA(:,6:8)=DATAF_av(:,5:7);
DATA(end+1,1) = DATADos(end,1)+(DATADos(end,1)-DATADos(end-1,1));
DATA(end,2:4) = DATADos(end-1,5:7)+(DATADos(end-1,5:7)-DATADos(end-2,5:7));
DATA(end,5) = DATAF_av(end-1,1)+(DATAF_av(end-1,1)-DATAF_av(end-2,1));
DATA(end,6:8) = DATAF_av(end-1,5:7)+(DATAF_av(end-1,5:7)-DATAF_av(end-2,5:7));

for i=2:sizeDATA
    % capteur sur le dos du joueur
    orientationDos(i,1:3) = DATA(i,2:4)*(1/FrHZ);
    orientationDos1(i,1:3)= orientationDos1(i-1,1:3)+orientationDos(i,1:3);
    % capteur avant cadre
    orientationCadre(i,1:3) = DATA(i,6:8)*(1/FrHZ);
    orientationCadre1(i,1:3)= orientationCadre1(i-1,1:3)+orientationCadre(i,1:3);
end


leftDistance=0;
rightDistance=0;
delta_distance=0;
Pos_x=0;
Pos_y=0;
Pos_x1=0;
Pos_y1=0;
for v=1:numel(DATAR1(:,1))-1
%%%%%%%%% calcul de la position du fauteil sur un plan 
leftDistance(v+1,1)=(meandist_R2(v+1)-meandist_R2(v));
rightDistance(v+1,1)=(meandist_R1(v+1)-meandist_R1(v));
delta_distance(v,1)=(leftDistance(v,1)+rightDistance(v,1))/2;
Pos_x(v+1,1)=cos((deg2rad(cadre_eul(v,3))))*delta_distance(v,1)+Pos_x(v,1);
Pos_y(v+1,1)=sin((deg2rad(cadre_eul(v,3))))*delta_distance(v,1)+Pos_y(v,1);
Pos_x1(v+1,1)=cos(deg2rad(orientationCadre1(v,3)))*delta_distance(v,1)+Pos_x1(v,1);% position selon x calculé avec le clacul de l'orientation du fauteuil par le gyro du cadre
Pos_y1(v+1,1)=sin(deg2rad(orientationCadre1(v,3)))*delta_distance(v,1)+Pos_y1(v,1);% position selon y calculé avec le clacul de l'orientation du fauteuil par le gyro du cadre
end
% hold on
% plot(Pos_x1,Pos_y1,'DisplayName','Trajectoire fauteuil selon IMU cadre')
% % plot(Pos_x,Pos_y,'DisplayName','Trajectoire fauteuil selon euler cadre')
% legend
% ylabel 'Déplacement en m'
% xlabel 'Déplacement en m'
% grid on
%% 
% Dossier de sauvgarde
Info_path=split(path,'\');
info_name = split(Info_path(end-1,1),'_');
name=cell2mat(info_name(1));
name_dos=cell2mat([info_name(1) info_name(2)]);
origindir=strcat(cell2mat(Info_path(1)),'\',cell2mat(Info_path(2)),'\',cell2mat(Info_path(3)),'\',cell2mat(Info_path(4)),'\',cell2mat(Info_path(5)),'\',cell2mat(Info_path(6)),'\');
figurepath = strcat(cell2mat(Info_path(1)),'\',cell2mat(Info_path(2)),'\',cell2mat(Info_path(3)),'\',cell2mat(Info_path(4)),'\',cell2mat(Info_path(5)),'\',cell2mat(Info_path(6)),'\',cell2mat(Info_path(7)),'\');

rr=dir([cd '\Acc_push_piste_' name '.mat']);
if ~isempty(rr)
load([cd '\Acc_push_piste_' name '.mat'])
end

%% Evenements 

% départ du test
seuil = max(DATAmean1)/2;
% DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81;
start1 = find(DATAmean1>seuil);
start = start1(1,1);
while DATAmean1(start,1)>0.21
    start=start-1;
end

                        StartDep=start;
                        % Vérification Graphique du Start des sprints
                        FIG_DECOLLAGE = figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7],'Color',[0.76 0.64 0.58],'NumberTitle','off','Name','Vérification du Point de départ du Déplacement (StartDep)','MenuBar','none');movegui('center');
                        FIGAXIS = axes('Units','Normalized','Color',[1 0.9686 0.9216],'Position',[0.1 0.1 0.65 0.8],'XColor',[.19 .19 .19],'YColor',[.19 .19 .19]);
                        Curve = plot(DATAmean1(:,1),'o','MarkerFaceColor','k','MarkerSize',2); hold on; ZOOM = 0;grid on
                        Marker = plot(StartDep,DATAmean1(StartDep),'o','color',[1 0.3 0],'LineWidth',1.5,'MarkerSize',7);
                        title({[cell2mat(nom) '-Essai : 1 ']}) ;
                        BUT_VALUE_Title = uicontrol('Units','Normalized','BackgroundColor',[.60 .60 .60],'FontSize',12,'Position',[.8 .80 .15 .05],'String','Coordonées','Style','text',...
                            'callback','');
                        BUT_VALUE_1 = uicontrol('Units','Normalized','BackgroundColor',[.80 .80 .80],'FontSize',12,'Position',[.8 .75 .15 .05],'String',num2str(3*FrHZ),'Style','text',...
                            'callback','');
                        BUT_VALUE_2 = uicontrol('Units','Normalized','BackgroundColor',[.80 .80 .80],'FontSize',12,'Position',[.8 .70 .15 .05],'String',num2str(DATAmean1(StartDep)),'Style','text',...
                            'callback','');
                        BUT_ZOOM = uicontrol('Units','Normalized','BackgroundColor',[.60 .60 .60],'FontSize',12,'Position',[.8 .55 .15 .125],'String','ZOOMER','Style','pushbutton',...
                           'callback','if ZOOM == 0;zoom on;ZOOM = 1;set(BUT_ZOOM,''String'',''DEZOOMER'');elseif ZOOM == 1;zoom off;xlim auto;ylim auto;ZOOM = 0;set(BUT_ZOOM,''String'',''ZOOMER'');end');
                        BUT_MODIF = uicontrol('Units','Normalized','BackgroundColor',[0.94 0.67 0],'FontSize',12,'Position',[.8 .30 .15 .125],'String','MODIFIER','Style','pushbutton',...
                            'callback','[Time,Value]=ginput(1);StartDep = round(Time);delete(Marker);Marker = plot(StartDep,DATAmean1(StartDep),''o'',''color'',[1 0.3 0],''LineWidth'',1.5,''MarkerSize'',7);set(BUT_VALUE_1,''String'',num2str(StartDep));set(BUT_VALUE_2,''String'',num2str(DATAmean1(StartDep)));');
                        
                        %                         if E==1
%                             BUT_SEUIL = uicontrol('Units','Normalized','BackgroundColor',[0.94 0.67 0.4],'FontSize',12,'Position',[.8 .45 .15 .1],'String','MODIF SEUIL','Style','pushbutton',...
%                                 'callback','seuil = inputdlg(''Indiquer le seuil pour la detection auto des pics :'',''Modification du seuil'', [1 50]);seuil = str2num(seuil{:});clear Event Events start start1 meandist DATAmean; DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81; start1 = find(DATAmean1(:,1)>seuil); start1 = start1(1,1); while DATAmean1(start1,1)>0.21; start1=start1-1;end; start = start1; Event=[0 0]; Events=0; count=0; count2=0; count1 = 2; for e=start:size(DATAmean1); DATAmean(e,1)=(DATAR1(e,13)+DATAR2(e,13))/2; if DATAmean(e,1) > seuil; count = count+1; Event(count,:) = [DATAmean(e,1) e];  end; end; clear count;count=1; for b=2:size(Event)-1; clear Event1 Ind_event;if Event(b,2)-Event(b-1,2)>7*FrHZ;count=count+1; Event1 = Event(b,2);while DATAmean1(Event1,1)>0.21;Event1=Event1-1;end; Events(count,1) = Event1; end; end; Events(1,1)=start; Events(end+1,1)=length(DATAmean1); meandist(1,1:size(Events)-1)={0};');
%                         end
                        BUT_VALID = uicontrol('Units','Normalized','BackgroundColor',[0.50 0.87 0.30],'FontSize',12,'Position',[.8 .15 .15 .125],'String','VALIDER','Style','pushbutton',...
                            'callback','close');
                        waitfor(FIG_DECOLLAGE);

                        start=StartDep;

% fin du test
End_event1= str2double(cell2mat(inputdlg(['durée de l exercice en secondes : ' name_dos])));
            meandist={0};
        for e=start:sizeDATA
            meandist((e-start)+2,1)={(DATAmean1(e,1)*((1/3600)/FrHZ)*1000)+cell2mat(meandist((e-start)+1,1))};
        end

    mat_dist=cell2mat(meandist(:,1));


if End_event1>29 && End_event1<80
    condition = 'resistance';
    Ind_end=find(mat_dist>=400);  %%%%%%%% Distance du 20m pour détecter la fin du sprint
    End_event=Ind_end(1,1)+start;
elseif End_event1<20
    condition = 'Sprint';
        Ind_end=find(mat_dist>=100);  %%%%%%%% Distance du 20m pour détecter la fin du sprint
    End_event=Ind_end(1,1)+start;
elseif End_event1>100
    condition = 'Endurance';
    End_event = start + round(End_event1*FrHZ);
end


% vitesses 
Vitdroite = DATAR1(start:End_event,13);
Vitgauche = DATAR2(start:End_event,13);
Vittot = DATAmean1(start:End_event,1);
% acceleration
Adroite = AccR1(start : End_event,1);
Agauche = AccR2(start : End_event,1);
Atot = Accmean(start : End_event,1);
% distances
Disttot = cumtrapz(Vittot/3.6)/FrHZ;
Distdroite = cumtrapz(Vitdroite/3.6)/FrHZ;
Distgauche = cumtrapz(Vitgauche/3.6)/FrHZ;
% orientation du dos
mvttronc = -(Dos_eul(start : End_event,2) + 110); 
% time
time_unit = 1/FrHZ;
time = time_unit : time_unit: time_unit*(size(Vitdroite,1));time = time';

% Cadre
cadre_acc = DATAF_av(start:End_event,2:4);
%% Detection des pics 
clear peak_minR1 peak_maxR1 peak_minR2 peak_maxR2 peak_accR1 peak_accR2 DATAfilt_R1 DATAfilt_R2 z x pminr1 pminr2 pmaxr1 pmaxr2 str color VPPP VP z
str = '#BCD4E6';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
fig_2 = figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7],'Color',color,'NumberTitle','off','Name','Vérification du Point de départ du Déplacement (StartDep)','MenuBar','none');movegui('center');
FIGAXIS = axes('Units','Normalized','Color',[1 0.9686 0.9216],'Position',[0.1 0.1 0.65 0.8],'XColor',[.19 .19 .19],'YColor',[.19 .19 .19]);set(gcf, 'WindowState', 'maximized');
% Détection auto des pics
z=1;x(1,1)=start;x(2,1)=End_event; 
VPPP=[Vitdroite Vitgauche]; % regroupement des vitesses droite et gauche dans l'interval du test (entre start et End_event)
%%%%%%% roue droite
[~,peak_maxR1] = findpeaks(VPPP(:,1), 'MinPeakWidth',12, 'MinPeakDistance', 20, 'MinPeakProminence', 0.5);
for i=1:length(peak_maxR1)-1
    [~,peak_minR1(i+1,1)]=min(VPPP(peak_maxR1(i,1):peak_maxR1(i+1,1),1));
    peak_minR1(i+1,1)=peak_minR1(i+1,1)+peak_maxR1(i,1)-1;
end
[~,peak_minR1(1,1)]=min(VPPP(1:peak_maxR1(1),1));
[~,peak_minR1(end+1,1)]=min(VPPP(peak_maxR1(end):end,1));
peak_minR1(end,1)=peak_minR1(end,1)+peak_maxR1(end)-1;
%%%%%% roue gauche
[~,peak_maxR2] = findpeaks(VPPP(:,2), 'MinPeakWidth',12, 'MinPeakDistance', 20, 'MinPeakProminence', 0.5);
for i=1:length(peak_maxR2)-1
    [~,peak_minR2(i+1,1)]=min(VPPP(peak_maxR2(i,1):peak_maxR2(i+1,1),2));
    peak_minR2(i+1,1)=peak_minR2(i+1,1)+peak_maxR2(i,1)-1;
end
[~,peak_minR2(1,1)]=min(VPPP(1:peak_maxR2(1),2));
[~,peak_minR2(end+1,1)]=min(VPPP(peak_maxR2(end):end,2));
peak_minR2(end,1)=peak_minR2(end,1)+peak_maxR2(end)-1;
%%%%%%%%%
pminr1=peak_minR1;pminr2=peak_minR2;pmaxr1=peak_maxR1;pmaxr2=peak_maxR2;
clear peak_minR1 peak_maxR1 peak_minR2 peak_maxR2
% Vérification et correction

sp1=subplot(2,1,1);
plot(VPPP(:,1),'r','DisplayName','Roue droite')
hold on
ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),'o','MarkerFaceColor','m','MarkerSize',10,'DisplayName','right pic min');
ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),'o','MarkerFaceColor','c','MarkerSize',10,'DisplayName','right pic max');
legend ('Position',[.83 .60 .06 .06])
sp2= subplot(2,1,2);
plot(VPPP(:,2),'DisplayName','Roue gauche')
hold on
ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),'o','MarkerFaceColor','r','MarkerSize',10,'DisplayName','left pic min');
ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),'o','MarkerFaceColor','b','MarkerSize',10,'DisplayName','left pic max');   
legend ('Position',[.83 .18 .06 .06])
% roue droite
BUTADMAXR1 = uicontrol('Units','Normalized','BackgroundColor','c','FontSize',12,'Position',[.0 .835 .12 .08],'String','Ajouter pic max R1','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pmaxr1(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1));pmaxr1(end,1)=pmaxr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pmaxr1(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1));pmaxr1(end,1)=pmaxr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp1;delete(ppmax1);hold on;pmaxr1=sort(pmaxr1) ;sp1=subplot(2,1,1); hold on;ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),''o'',''MarkerFaceColor'',''c'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
BUTADMINR1 = uicontrol('Units','Normalized','BackgroundColor','m','FontSize',12,'Position',[.0 .755 .12 .08],'String','Ajouter pic min R1','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pminr1(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1));pminr1(end,1)=pminr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pminr1(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1));pminr1(end,1)=pminr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp1;delete(ppmin1);hold on;pminr1=sort(pminr1) ;sp1=subplot(2,1,1); hold on;ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),''o'',''MarkerFaceColor'',''m'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');
BUTSUPMAXR1 = uicontrol('Units','Normalized','BackgroundColor','c','FontSize',12,'Position',[.0 .675 .12 .08],'String','Supprimer pic max R1','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)==max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)==max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pmaxr1(find(pmaxr1==indexpos(i,1)),:)=[];end;sp1;delete(ppmax1);hold on;pmaxr1=sort(pmaxr1) ;sp1=subplot(2,1,1); hold on;ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),''o'',''MarkerFaceColor'',''c'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
BUTSUPMINR1 = uicontrol('Units','Normalized','BackgroundColor','m','FontSize',12,'Position',[.0 .595 .12 .08],'String','Supprimer pic min R1','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)==min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)==min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pminr1(find(pminr1==indexpos(i,1)),:)=[];end;sp1;delete(ppmin1);hold on;pminr1=sort(pminr1) ;sp1=subplot(2,1,1); hold on;ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),''o'',''MarkerFaceColor'',''m'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');
% roue gauche
BUTADMAXR2 = uicontrol('Units','Normalized','BackgroundColor','b','FontSize',12,'Position',[.0 .37 .12 .08],'String','Ajouter pic max R2','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pmaxr2(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2));pmaxr2(end,1)=pmaxr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pmaxr2(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2));pmaxr2(end,1)=pmaxr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp2;delete(ppmax2);hold on;pmaxr2=sort(pmaxr2) ;sp2=subplot(2,1,2); hold on;ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),''o'',''MarkerFaceColor'',''b'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
BUTADMINR2 = uicontrol('Units','Normalized','BackgroundColor','r','FontSize',12,'Position',[.0 .29 .12 .08],'String','Ajouter pic min r2','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pminr2(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2));pminr2(end,1)=pminr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pminr2(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2));pminr2(end,1)=pminr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp2;delete(ppmin2);hold on;pminr2=sort(pminr2) ;sp2=subplot(2,1,2); hold on;ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),''o'',''MarkerFaceColor'',''r'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');
BUTSUPMAXR2 = uicontrol('Units','Normalized','BackgroundColor','b','FontSize',12,'Position',[.0 .21 .12 .08],'String','Supprimer pic max R2','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)==max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)==max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pmaxr2(find(pmaxr2==indexpos(i,1)),:)=[];end;sp2;delete(ppmax2);hold on;pmaxr2=sort(pmaxr2) ;sp2=subplot(2,1,2); hold on;ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),''o'',''MarkerFaceColor'',''b'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .18 .06 .06]);');
BUTSUPMINR2 = uicontrol('Units','Normalized','BackgroundColor','r','FontSize',12,'Position',[.0 .13 .12 .08],'String','Supprimer pic min R2','Style','pushbutton',...
    'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)==min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)==min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pminr2(find(pminr2==indexpos(i,1)),:)=[];end;sp2;delete(ppmin2);hold on;pminr2=sort(pminr2) ;sp2=subplot(2,1,2); hold on;ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),''o'',''MarkerFaceColor'',''r'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .18 .06 .06]);');
% validation
BUT_VALID = uicontrol('Units','Normalized','BackgroundColor',[0.50 0.87 0.30],'FontSize',12,'Position',[.0 .01 .12 .10],'String','VALIDER','Style','pushbutton',...
        'callback','close');
waitfor(fig_2);

peak_minR1(:,1) = pminr1(:,1); 
peak_minR1(:,2) = VPPP(peak_minR1(:,1),1);
peak_minR2(:,1) = pminr2(:,1); 
peak_minR2(:,2) = VPPP(peak_minR2(:,1),2);
peak_maxR1(:,1) = pmaxr1(:,1);  
peak_maxR1(:,2) = VPPP(peak_maxR1(:,1),1);
peak_maxR2(:,1) = pmaxr2(:,1); 
peak_maxR2(:,2) = VPPP(peak_maxR2(:,1),2);

clear pminr1 pminr2 pmaxr1 pmaxr2  

%% Calculs 


% Analyse par cycle :
switch condition
    case 'resistance'
        E=2;
    case 'Sprint'
        E=1;
    case 'Endurance'
        E=3;
end
del=length(peak_maxR1)-length(peak_maxR2);
if del <= 0
    sizepeak=length(peak_maxR1);
elseif del>0
    sizepeak=length(peak_maxR2);
end
Acc_push(1,1)={'Acc Poussée Roue 1'};
Acc_push(1,4)={condition};
Acc_push(9,1)={'Acc Poussée Roue 2'};
Acc_push(17,1)={'temps poussée Roue 1'};
Acc_push(17,4)={condition};
Acc_push(25,1)={'temps poussée Roue 2'};
Acc_push(33,1)={'temps recouvrement Roue 1'};
Acc_push(33,4)={condition};
Acc_push(41,1)={'temps recouvrement Roue 2'};
Acc_push(49,1)={'temps cycle Roue 1'};
Acc_push(49,4)={condition};
Acc_push(57,1)={'temps cycle Roue 2'};
Acc_push(65,1) = {'Cadance roue droite '};
Acc_push(65,4)={condition};
Acc_push(73,1) = {'Cadance roue gauche'};
Acc_push(73,4)={condition};
Acc_push(81,1) = {'Vitesse max par poussée roue 1'};
Acc_push(81,4)={condition};
Acc_push(89,1) = {'Vitesse max poussée roue 2'};
Acc_push(97,1) = {'Vitesse min par poussée roue 1 '};
Acc_push(97,4)={condition};
Acc_push(105,1) = {'Vitesse min par poussée roue 2'};
Acc_push(113,4)={condition};
Acc_push(113,1) = {'coefficient de correlation droite gauche '};
Acc_push(121,1) = {'Angle de poussée roue 2'};
for c=1:sizepeak-1
    % temps poussée
    tpsPR1(E,c)={(peak_maxR1(c,1)-peak_minR1(c,1))/FrHZ};
    tpsPR2(E,c)={(peak_maxR2(c,1)-peak_minR2(c,1))/FrHZ};
    tpsP_mean(E,c)=(cell2mat(tpsPR1(E,c))+cell2mat(tpsPR2(E,c)))/2;
    % temps recouvrement
    tpsReR1(E,c)={(peak_minR1(c+1,1)-peak_maxR1(c,1))/FrHZ};
    tpsReR2(E,c)={(peak_minR2(c+1,1)-peak_maxR2(c,1))/FrHZ};
    tpsRe_mean(E,c)=(cell2mat(tpsReR1(E,c))+cell2mat(tpsReR2(E,c)))/2;
    % temps cycles
    tpscycleR1(E,c)={(peak_minR1(c+1,1)-peak_minR1(c,1))/FrHZ};
    tpscycleR2(E,c)={(peak_minR2(c+1,1)-peak_minR2(c,1))/FrHZ};
    tpscycle_mean(E,c)=(cell2mat(tpscycleR1(E,c))+cell2mat(tpscycleR2(E,c)))/2;
    % Asymétrie en %age
        if peak_maxR1(c,2)-peak_maxR2(c,2)>=0
    Asy(E,c)=(peak_maxR1(c,2)-peak_maxR2(c,2))/peak_maxR1(c,2);
        else
    Asy(E,c)=(peak_maxR2(c,2)-peak_maxR1(c,2))/peak_maxR2(c,2);
        end   
   % Acceleration par poussée
   Acc_push(E+1,c)={(((peak_maxR1(c,2)/3.6)-(peak_minR1(c,2)/3.6)))/((peak_maxR1(c,1)-peak_minR1(c,1))/FrHZ)};
   Acc_push(E+9,c)={(((peak_maxR2(c,2)/3.6)-(peak_minR2(c,2)/3.6)))/((peak_maxR2(c,1)-peak_minR2(c,1))/FrHZ)};
   Acc_push(E+17,c)=tpsPR1(E,c);
   Acc_push(E+25,c)=tpsPR2(E,c);
   Acc_push(E+33,c)=tpsReR1(E,c);
   Acc_push(E+41,c)=tpsReR2(E,c);
   Acc_push(E+49,c)=tpscycleR1(E,c);
   Acc_push(E+57,c)=tpscycleR2(E,c);

   Acc_push_mean(E,c) = {(cell2mat(Acc_push(E+1,c))+ cell2mat(Acc_push(E+9,c)))/2};
%    Force_push(E,c) = {cell2mat(Acc_push_mean(E,c))*M};
%    Power_push(E,c) = {cell2mat(Force_push (E,c)) * ((peak_maxR1(c,2)/3.6)-(peak_minR1(c,2)/3.6))};

%                            Angle_pushR1(E,c) = {rad2deg((((peak_maxR1(c,2)/3.6)+(peak_minR1(c,2)/3.6))/2)/rayon_roue)*(cell2mat(tpsPR1(E,c)))};
%                            Angle_pushR2(E,c) = {rad2deg((((peak_maxR2(c,2)/3.6)+(peak_minR2(c,2)/3.6))/2)/rayon_roue)*(cell2mat(tpsPR2(E,c)))};

   Acc_push(E+65,c) = {60/cell2mat(tpscycleR1(E,c))};
   Acc_push(E+73,c) = {60/cell2mat(tpscycleR2(E,c))};

   Acc_push(E+81,c) = {peak_maxR1(c,2)};
   Acc_push(E+89,c) = {peak_maxR2(c,2)}; 
   Acc_push(E+97,c) = {peak_minR1(c,2)};
   Acc_push(E+105,c) = {peak_minR2(c,2)};
%                            Acc_push(E+113,c) = Angle_pushR1(E,c);
%                            Acc_push(E+121,c) = Angle_pushR2(E,c);
    pic_max(c,:)=round((peak_maxR1(c,:)+peak_maxR2(c,:))/2);
    pic_start(c,:)=round((peak_minR1(c,:)+peak_minR2(c,:))/2);
 cadence_evolution(c,:) = ( 60/cell2mat(tpscycleR1(E,c))+60/cell2mat(tpscycleR2(E,c)) )/2 ;
 
 %% rangement de la vitesse par cycle 
 % droite
    clear uuuu bbb b vit_rot_tronc
    uuuu = Vitdroite(peak_minR1(c,1):peak_minR1(c+1,1),1);%-orientation1(pic_start(i,1)+Events(end-1,1),1);
    clear y x xx data_norm ampmaxcpc ampmincpc ampmaxcadreacc
    y = length(uuuu);
                x = 0:1/FrHZ:((y-1)/FrHZ);
                xx = 0:1*x(size(x,2))/100:x(size(x,2));
                data_norm = interp1(x,uuuu,xx);
                %----------- Stockage Normalisation ------------
                Cycle_Norm_Rdroite(1:101,c)=num2cell(data_norm');
   % gauche
    clear uuuu bbb b vit_rot_tronc
    uuuu = Vitgauche(peak_minR2(c,1):peak_minR2(c+1,1),1);%-orientation1(pic_start(i,1)+Events(end-1,1),1);
    clear y x xx data_norm ampmaxcpc ampmincpc ampmaxcadreacc
    y = length(uuuu);
                x = 0:1/FrHZ:((y-1)/FrHZ);
                xx = 0:1*x(size(x,2))/100:x(size(x,2));
                data_norm = interp1(x,uuuu,xx);
                %----------- Stockage Normalisation ------------
                Cycle_Norm_Rgauche(1:101,c)=num2cell(data_norm');
                
 % mesure du taux de ressemblance entre les cyles de propulsion droite et
 % gauche (coeficient de correlation)
 sumx=0;
 sumxC = 0;
 sumy = 0;
 sumyC = 0;
 sumxy=0;
xmean = mean(cell2mat(Cycle_Norm_Rdroite(:,c)));
ymean = mean(cell2mat(Cycle_Norm_Rgauche(:,c)));

 for i=1:101
    xnum=cell2mat(Cycle_Norm_Rdroite(i,c));
    ynum=cell2mat(Cycle_Norm_Rgauche(i,c));
    sumx = (xnum-xmean)+sumx;
    sumy = (ynum - ymean)+sumy;
    sumxy = ((xnum-xmean)*(ynum - ymean))+sumxy;
    sumxC = ((xnum-xmean)^2) + sumxC;
    sumyC = ((ynum-ymean)^2) + sumyC;
 end
 Cc(c) = sumxy/(sqrt(sumxC*sumyC));    
     
%       % Synchro
%       clear acor lag
%     [acor, lag] = xcorr(cell2mat(Cycle_Norm_Rgauche(:,c)),cell2mat(Cycle_Norm_Rdroite(:,c)));
%                 [~, I] = max (abs (acor));
%             lagDiff = lag (I);

Acc_push(E+113,c) = {Cc(c)};
  %%              
end    
E=1;
cadenceR1(E,1)=(sizepeak-1)/(peak_minR1(end,1)/FrHZ)*60;
cadenceR2(E,1)=(sizepeak-1)/(peak_minR2(end,1)/FrHZ)*60;
cadence(E,1)= (cadenceR1(E,1)+cadenceR2(E,1))/2;



% Figures 
% plot(peak_maxR1(1:end-1,1)/FrHZ,cell2mat(Acc_push(E+65,:)))

for i=1:length(pic_max)-1
    pos_peak_max(i,1)=round(((pic_max(i,1)-pic_start(i,1))/(pic_start(i+1,1)-pic_start(i,1)))*100);
%     POS_PM_Stat(i,E)=pos_peak_max(i,1);
    clear uuuu bbb b vit_rot_tronc
    uuuu = mvttronc(pic_start(i,1):pic_start(i+1,1),1);%-orientation1(pic_start(i,1)+Events(end-1,1),1);
    clear y x xx data_norm ampmaxcpc ampmincpc ampmaxcadreacc
    y = length(uuuu);
                x = 0:1/FrHZ:((y-1)/FrHZ);
                xx = 0:1*x(size(x,2))/100:x(size(x,2));
                data_norm = interp1(x,uuuu,xx);
                %----------- Stockage Normalisation ------------
                Cycle_Norm_dos(1:101,i)=num2cell(data_norm');
           ampmaxcpc=max(cell2mat(Cycle_Norm_dos(:,i)));
        ampmincpc=min(cell2mat(Cycle_Norm_dos(:,i)));
        ampcpc(i,1)=abs(abs(ampmaxcpc)-abs(ampmincpc));
        ampmaxcadreacc = max(cadre_acc(pic_start(i,1):pic_start(i+1,1),3));  
         ampmincadreacc = min(cadre_acc(pic_start(i,1):pic_start(i+1,1),3));
        ampcadre_acc(i,1)=abs(abs(ampmaxcadreacc)+abs(ampmincadreacc));
         
       
            
% result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).(strcat('amplitude_cycle_sprint_',num2str(E))).data=Cycle_Norm_dos;
end

% en fonction de l'exercice 

        % general :

        [Vmax, tpsframeVmax] = max(Vittot);
        Vmean = mean(Vittot);
        tpsVmax = tpsframeVmax/FrHZ;
        DistVmax = Disttot(tpsframeVmax);
%         Dist_test = Disttot(end,1);
        tpstest = time(end);
%         Dist_test_droite = Distdroite(end,1);
%         Dist_test_gauche = Distgauche(end,1);        
        Ameantot = mean(Atot);
        AtoVmax = mean(Atot(1:tpsframeVmax));
% Figures
        figure(1)
set(gcf, 'WindowState', 'maximized');
hold on 
plot(time,Vitdroite,'DisplayName','Vitesse Roue droite','LineWidth',3)
plot(time,Vitgauche,'DisplayName','Vitesse Roue gauche','LineWidth',3)
legend('location','southeast')
set(gca,'fontsize',20);
ylabel('vitesse en Km/h','FontSize',25)
xlabel('Temps (s)','FontSize',25)
switch condition
    case 'resistance'
        xticks(0:2:60)
        xlim([0 60])
    case 'Sprint'
        xticks(0:1:15)
        xlim([0 15])
    case 'Endurance'
        xlim([0 300])
        xticks(0:15:300)
end
grid on 
hold off 
saveas(figure(1),[figurepath name 'courbe_vitesse_' condition],'fig')
saveas(figure(1),[figurepath name 'courbe_vitesse_' condition],'jpg')
imageFileName2= [figurepath name 'courbe_vitesse_' condition '.jpg'];
close all
      


switch condition
    case 'resistance'
E=2;
ind_100 = find(Disttot>=100);ind_100=ind_100(1,1);
tps100m = time(ind_100,1);
ind_200 = find(Disttot>=200);ind_200=ind_200(1,1);
tps200m = time(ind_200,1);

        % Indice de fatigue
        if End_event1>59
            Indfatigue60 = ((Vmax-(peak_maxR1(end,2) + peak_maxR2(end,2))/2)/Vmax);
            pos_peak30 = find(peak_maxR1(:,1)<30*FrHZ); pos_peak30=pos_peak30(end);
            Indfatigue30 = ((Vmax-(peak_maxR1(pos_peak30,2) + peak_maxR2(pos_peak30,2))/2)/Vmax);
        else
            Indfatigue30 = ((Vmax-(peak_maxR1(end,2) + peak_maxR2(end,2))/2)/Vmax);
        end
        % Cinematique du tronc
        [~,ind_cyclemax] = max(pic_max(:,2));
        Cycle_dos_start(:,1)=mean(cell2mat(Cycle_Norm_dos(:,1:ind_cyclemax)),2);
        Cycle_dos_start(:,2)=std(cell2mat(Cycle_Norm_dos(:,1:ind_cyclemax))');
        Cycle_dos_stab(:,1)=mean(cell2mat(Cycle_Norm_dos(:,ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3)))),2);
        Cycle_dos_stab(:,2)=std(cell2mat(Cycle_Norm_dos(:,ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3))))');
        Cycle_dos_end(:,1)=mean(cell2mat(Cycle_Norm_dos(:,end-10:end)),2);
        Cycle_dos_end(:,2)=std(cell2mat(Cycle_Norm_dos(:,end-10:end))');
        T=0:1:100;
        mean_pos_start=round(mean(pos_peak_max(1:ind_cyclemax)));
        mean_pos_stab=round(mean(pos_peak_max(ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3)))));
        mean_pos_end=round(mean(pos_peak_max(end-10:end)));
        % figure 
          figure(4)
          set(gcf, 'WindowState', 'maximized');
        hold on
        grid on 
        boundedline (T,Cycle_dos_start(:,1),Cycle_dos_start(:,2),'r','alpha') % phase de demarrage
        boundedline (T,Cycle_dos_stab(:,1),Cycle_dos_stab(:,2),'g','alpha') % phase stabilisé
        boundedline (T,Cycle_dos_end(:,1),Cycle_dos_end(:,2),'b','alpha') % phase stabilisé
        xline(mean_pos_start,'r','LineWidth',4)
        xline(mean_pos_stab,'g','LineWidth',4)        
        xline(mean_pos_end,'b','LineWidth',4)
        legend('Amplitude du tronc du début à Vmax','Pic max du cycle au début','Amplitude du tronc au milieu du test','Pic max du cycle au mileu','Amplitude du tronc à la fin','Pic max du cycle à la fin','location','southeast')
        set(gca,'fontsize',20);
        ylabel('Angle en °','FontSize',25)
        xlabel('pourcentage du cycle de propulsion','FontSize',25)
        hold off
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'fig')
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'jpg')
        imageFileName1= [figurepath name 'courbe_amplitudetronc_' condition '.jpg'];
        close all

        % amplitude tronc 
        ampmaxstart=max(Cycle_dos_start(:,1));
        ampminstart=min(Cycle_dos_start(:,1));
        ampstart=abs(abs(ampmaxstart)-abs(ampminstart));
        ampmaxstab=max(Cycle_dos_stab(:,1));
        ampminstab=min(Cycle_dos_stab(:,1));
        ampstab=abs(abs(ampmaxstab)-abs(ampminstab));
        ampmaxend=max(Cycle_dos_end(:,1));
        ampminsend=min(Cycle_dos_end(:,1));
        ampend=abs(abs(ampmaxend)-abs(ampminsend));
        
                % temporel
        tps_cycle = mean(tpscycle_mean(E,:));
        tps_push = mean(tpsP_mean(E,:));
        tps_reco = mean(tpsRe_mean(E,:));
       
        tps_cyclestart = mean(tpscycle_mean(E,1:ind_cyclemax));
        tps_pushstart = mean(tpsP_mean(E,1:ind_cyclemax));
        tps_recostart = mean(tpsRe_mean(E,1:ind_cyclemax));

        tps_cyclestab = mean(tpscycle_mean(E,ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3))));
        tps_pushstab = mean(tpsP_mean(E,ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3))));
        tps_recostab = mean(tpsRe_mean(E,ind_cyclemax:ind_cyclemax+(round(length(pic_max)/3))));

        tps_cycleend = mean(tpscycle_mean(E,end-10:end));
        tps_pushend = mean(tpsP_mean(E,end-10:end));
        tps_recoend = mean(tpsRe_mean(E,end-10:end));
        
        % RECAP resistance
        RECAP_sprint(:,1) = Vmax;
        RECAP_sprint(:,2) = Vmean;
        RECAP_sprint(:,3) = (peak_maxR1(end,2) + peak_maxR2(end,2))/2;
        RECAP_sprint(:,4) = Indfatigue60;
        RECAP_sprint(:,5) = (peak_maxR1(pos_peak30,2) + peak_maxR2(pos_peak30,2))/2;
        RECAP_sprint(:,6) = Indfatigue30;
        RECAP_sprint(:,7) = tpsVmax;
        RECAP_sprint(:,8) = DistVmax;
        RECAP_sprint(:,9) = tps100m;
        RECAP_sprint(:,10)= tps200m;
        RECAP_sprint(:,11)= tpstest;
        RECAP_sprint(:,12)= tps_pushstart;
        RECAP_sprint(:,13)= tps_pushstab;
        RECAP_sprint(:,14)= tps_pushend;
        RECAP_sprint(:,15)= tps_recostart;
        RECAP_sprint(:,16)= tps_recostab;
        RECAP_sprint(:,17)= tps_recoend;
        RECAP_sprint(:,18)= tps_cyclestart;
        RECAP_sprint(:,19)= tps_cyclestab;
        RECAP_sprint(:,20)= tps_cycleend;
        RECAP_sprint(:,21)= ampstart;
        RECAP_sprint(:,22)= ampstab;        
        RECAP_sprint(:,23)= ampend;        
        RECAP_sprint=RECAP_sprint';
        
        
x = Pos_x;
y = -Pos_y;
y(end) = NaN;
c = DATAmean1;
[~, Vmaxtot]= max(DATAmean1);
pos_Vmax_x = x(Vmaxtot);
pos_Vmax_y = y(Vmaxtot);
pos_start_x = x(start);
pos_start_y = y(start);
figure(5)
patch(x, y, c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat')

colorbar;
hold on 
plot(pos_Vmax_x,pos_Vmax_y,'Marker','o','MarkerFaceColor','y','linewidth',10)
plot(pos_start_x,pos_start_y,'Marker','o','MarkerFaceColor','b','linewidth',10)
xlabel('Distance en m')
ylabel('Distance en m')
text(pos_Vmax_x,pos_Vmax_y,'\leftarrow Vmax','FontSize',14)
text(pos_start_x,pos_start_y,'\leftarrow Départ','FontSize',14)
        saveas(figure(5),[figurepath name 'courbe_trajectoire_' condition],'fig')
        saveas(figure(5),[figurepath name 'courbe_trajectoire_' condition],'jpg')
        imageFileName3= [figurepath name 'courbe_trajectoire_' condition '.jpg'];
close all


                    clear N FileName_Rapport
                    copyfile([cd '\Model_rapport_athle_piste.xlsx'],origindir);                    
                    FileName_Rapport = strcat(origindir,name,'_Rapports','.xlsx');
                    movefile (strcat(origindir,'\','Model_rapport_athle_piste.xlsx'),FileName_Rapport);   
                    xlswrite(FileName_Rapport,RECAP_sprint,'DATAS','F10');
                    xlswrite(FileName_Rapport,Acc_push,'Recap_poussées');
    xlswrite(FileName_Rapport,cadence_evolution','DATAS','A44');

                    imageSize=[482, 326];

                           xlsputimage(FileName_Rapport, imageFileName1, 'Rapport_resistance', 'A101', imageSize)
                                               imageSize=[482, 300];

                    xlsputimage(FileName_Rapport, imageFileName2, 'Rapport_resistance', 'A126', imageSize)
                    xlsputimage(FileName_Rapport, imageFileName3, 'Rapport_resistance', 'A148', imageSize)

        
    case 'Sprint'
          
              
        % Cinematique du tronc
        ind_cyclemax = 3;
        Cycle_dos_start(:,1)=mean(cell2mat(Cycle_Norm_dos(:,1:ind_cyclemax)),2);
        Cycle_dos_start(:,2)=std(cell2mat(Cycle_Norm_dos(:,1:ind_cyclemax))');
        Cycle_dos_end(:,1)=mean(cell2mat(Cycle_Norm_dos(:,end-5:end)),2);
        Cycle_dos_end(:,2)=std(cell2mat(Cycle_Norm_dos(:,end-5:end))');
        T=0:1:100;
        mean_pos_start=round(mean(pos_peak_max(1:ind_cyclemax)));
        mean_pos_end=round(mean(pos_peak_max(end-5:end)));
        figure(4)
        set(gcf, 'WindowState', 'maximized');
        hold on
        grid on
        boundedline (T,Cycle_dos_start(:,1),Cycle_dos_start(:,2),'r','alpha') % phase de demarrage
        boundedline (T,Cycle_dos_end(:,1),Cycle_dos_end(:,2),'b','alpha') % phase stabilisé
        xline(mean_pos_start,'r','LineWidth',4)
        xline(mean_pos_end,'b','LineWidth',4)
        legend('Amplitude du tronc au démarrage','Pic max du cycle au démarrage','Amplitude du tronc à la fin','Pic max du cycle à la fin','location','southeast')
        set(gca,'fontsize',20);
        ylabel('Angle en °','FontSize',25)
        xlabel('pourcentage du cycle de propulsion','FontSize',25)
        hold off
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'fig')
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'jpg')
        imageFileName1= [figurepath name 'courbe_amplitudetronc_' condition '.jpg'];
        close all
        
                % amplitude tronc 
        ampmaxstart=max(Cycle_dos_start(:,1));
        ampminstart=min(Cycle_dos_start(:,1));
        ampstart=abs(abs(ampmaxstart)-abs(ampminstart));
        ampmaxend=max(Cycle_dos_end(:,1));
        ampminsend=min(Cycle_dos_end(:,1));
        ampend=abs(abs(ampmaxend)-abs(ampminsend));

          % temporel        
        E=1;
        RECAP_cycle(E,1:3)   = tpsP_mean(E,1:3);
        RECAP_cycle(E,4)     = mean(tpsP_mean(E,1:3));
        RECAP_cycle(E,5)     = mean(tpsP_mean(E,end-5:end));
        RECAP_cycle(E,6)     = mean(tpsP_mean(E,:));
        RECAP_cycle(E,7:9)   = tpsRe_mean(E,1:3);
        RECAP_cycle(E,10)    = mean(tpsRe_mean(E,1:3));
        RECAP_cycle(E,11)    = mean(tpsRe_mean(E,end-5:end));
        RECAP_cycle(E,12)    = mean(tpsRe_mean(E,:));
        RECAP_cycle(E,13:15) = tpscycle_mean(E,1:3);
        RECAP_cycle(E,16)    = mean(tpscycle_mean(E,1:3));
        RECAP_cycle(E,17)    = mean(tpscycle_mean(E,end-5:end));
        RECAP_cycle(E,18)    = mean(tpscycle_mean(E,:));
        RECAP_cycle(E,19:21) = Asy(E,1:3);
        RECAP_cycle(E,22)    = mean(Asy(E,1:3));
        RECAP_cycle(E,23)    = mean(Asy(E,end-5:end));
        RECAP_cycle(E,24)    = cadence(E,1);
        
                %% RECAP : 
                            % vitesse max 3 1er
        vmaxR1start1(E,1)=peak_maxR1(1,1)/FrHZ;
        vmaxR1start1(E,2)=peak_maxR1(1,2);
        vmaxR2start1(E,1)=peak_maxR2(1,1)/FrHZ;
        vmaxR2start1(E,2)=peak_maxR2(1,2);
        Vmoymax1(E,1)=(peak_maxR1(1,2)+peak_maxR2(1,2))/2;

        vmaxR1start2(E,1)=peak_maxR1(2,1)/FrHZ;
        vmaxR1start2(E,2)=peak_maxR1(2,2);
        vmaxR2start2(E,1)=peak_maxR2(2,1)/FrHZ;
        vmaxR2start2(E,2)=peak_maxR2(2,2);        
        Vmoymax2(E,1)=(peak_maxR1(2,2)+peak_maxR2(2,2))/2;

        vmaxR1start3(E,1)=peak_maxR1(3,1)/FrHZ;
        vmaxR1start3(E,2)=peak_maxR1(3,2);
        vmaxR2start3(E,1)=peak_maxR2(3,1)/FrHZ;
        vmaxR2start3(E,2)=peak_maxR2(3,2);  
        Vmoymax3(E,1)=(peak_maxR1(3,2)+peak_maxR2(3,2))/2;
                                    % acc moy 3 1er
        accR1start(E,1) =mean(Adroite(1:peak_minR1(4,1),1));
        accR2start(E,1) =mean(Agauche(1:peak_minR2(4,1)));
        accstart(E,1) =mean(Atot(1:(peak_maxR2(3,1)+peak_maxR1(3,1))/2,1));

                                % Calcul entre les 5 avant-dreniers pics (phase stabilisée)
            % vitesse moy 5 Der
        vmoyR1stab(E,1) =mean(Vitdroite(peak_minR1(end-5,1):peak_minR1(end,1),1));
        vmoyR2stab(E,1) =mean(Vitgauche(peak_minR2(end-5,1):peak_minR2(end,1),1));
        vmoystart(E,1)=(Vmoymax1(E,1)+Vmoymax2(E,1)+Vmoymax3(E,1))/3;
            % vitesse max 5 Der
        vmaxR1stab(E,1) = mean(peak_maxR1(end-5:end,2));
        vmaxR2stab(E,1) = mean(peak_maxR2(end-5:end,2));
        vmaxmoystab(E,1) = (vmaxR1stab(E,1)+vmaxR2stab(E,1))/2;
            % acc moy 5 Der
        accR1stab(E,1) = mean(Adroite(peak_minR1(end-5,1):peak_minR1(end,1),1));
        accR2stab(E,1) = mean(Agauche(peak_minR2(end-5,1):peak_minR2(end,1),1));
        accstab(E,1)   = mean(Atot(peak_minR2(end-5,1):peak_minR2(end,1),1));

        RECAP_sprint(:,1) = vmaxR1start1(:,2);
        RECAP_sprint(:,2) = vmaxR2start1(:,2);
        RECAP_sprint(:,3) = Vmoymax1;
        RECAP_sprint(:,4) = vmaxR1start2(:,2);
        RECAP_sprint(:,5) = vmaxR2start2(:,2);
        RECAP_sprint(:,6) = Vmoymax2;
        RECAP_sprint(:,7) = vmaxR1start3(:,2);
        RECAP_sprint(:,8) = vmaxR2start3(:,2);
        RECAP_sprint(:,9) = Vmoymax3;
        RECAP_sprint(:,10)= vmoystart;
        RECAP_sprint(:,11)= vmaxR1stab;
        RECAP_sprint(:,12)= vmaxR2stab;
        RECAP_sprint(:,13)= vmaxmoystab;
        RECAP_sprint(:,14)= Vmax;
        RECAP_sprint(:,15)= Vmean;
        RECAP_sprint(:,16)= accstart;
        RECAP_sprint(:,17)= max(Atot);
        RECAP_sprint(:,18)= mean(Atot);
        RECAP_sprint(:,19)= accstab;
        RECAP_sprint(:,20)= DistVmax;        
        RECAP_sprint(:,21)= tpstest;
        RECAP_sprint(:,22)= tpsVmax;        
        RECAP_sprint(:,23:46)=RECAP_cycle;
        RECAP_sprint(:,47)= ampstart;
        RECAP_sprint(:,48)= ampend; 
x = Pos_x;
y = -Pos_y;
y(end) = NaN;
c = DATAmean1;
[~, Vmaxtot]= max(DATAmean1);
pos_Vmax_x = x(Vmaxtot);
pos_Vmax_y = y(Vmaxtot);
pos_start_x = x(start);
pos_start_y = y(start);
figure(5)
patch(x, y, c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat')

colorbar;
hold on 
plot(pos_Vmax_x,pos_Vmax_y,'Marker','o','MarkerFaceColor','y','linewidth',10)
plot(pos_start_x,pos_start_y,'Marker','o','MarkerFaceColor','b','linewidth',10)
xlabel('Distance en m')
ylabel('Distance en m')
text(pos_Vmax_x,pos_Vmax_y,'\leftarrow Vmax','FontSize',14)
text(pos_start_x,pos_start_y,'\leftarrow Départ','FontSize',14)

% coordinates of point 1 
x1 = pos_start_x;
y1 = pos_start_y;
  
% coordinates of point 2
x2 = x(End_event);
y2 = y(End_event);

        saveas(figure(5),[figurepath name 'courbe_trajectoire_' condition],'fig')
        saveas(figure(5),[figurepath name 'courbe_trajectoire_' condition],'jpg')
        imageFileName3= [figurepath name 'courbe_trajectoire_' condition '.jpg'];
close all


                    FileName_Rapport = strcat(origindir,name,'_Rapports','.xlsx');
                    xlswrite(FileName_Rapport,RECAP_sprint,'DATAS','F3');
                    xlswrite(FileName_Rapport,Acc_push,'Recap_poussées');
    xlswrite(FileName_Rapport,cadence_evolution','DATAS','A42');

                    imageSize=[482, 326];
                    xlsputimage(FileName_Rapport, imageFileName1, 'Rapport_sprint', 'A69', imageSize)
                    xlsputimage(FileName_Rapport, imageFileName2, 'Rapport_sprint', 'A100', imageSize)

    case 'Endurance'
        E=3;
                % Cinematique du tronc
        Cycle_dos_start(:,1)=mean(cell2mat(Cycle_Norm_dos(:,10:end)),2);
        Cycle_dos_start(:,2)=std(cell2mat(Cycle_Norm_dos(:,10:end))');
        T=0:1:100;
        mean_pos_start=round(mean(pos_peak_max(10:end)));
          figure(4)
        set(gcf, 'WindowState', 'maximized');
        hold on
        grid on
        boundedline (T,Cycle_dos_start(:,1),Cycle_dos_start(:,2),'r','alpha') % phase de demarrage
        xline(mean_pos_start,'r','LineWidth',4)
        legend('Amplitude du tronc','Pic max du cycle','location','southeast')
        set(gca,'fontsize',20);
        ylabel('Angle en °','FontSize',25)
        xlabel('pourcentage du cycle de propulsion','FontSize',25)
        hold off
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'fig')
        saveas(figure(4),[figurepath name 'courbe_amplitudetronc_' condition],'jpg')
        imageFileName1= [figurepath name 'courbe_amplitudetronc_' condition '.jpg'];
        close all
                % amplitude tronc 
        ampmaxstart=max(Cycle_dos_start(:,1));
        ampminstart=min(Cycle_dos_start(:,1));
        ampstart=abs(abs(ampmaxstart)-abs(ampminstart));
       
        % temporel
        tps_cycle = mean(tpscycle_mean(E,:));
        tps_push = mean(tpsP_mean(E,:));
        tps_reco = mean(tpsRe_mean(E,:));
        
        RECAP_sprint(:,1) = tps_push;
        RECAP_sprint(:,2) = tps_reco;
        RECAP_sprint(:,3) = tps_cycle;
        RECAP_sprint(:,4) = cadence;
        RECAP_sprint(:,5) = ampstart;
RECAP_sprint=RECAP_sprint';

        
        FileName_Rapport = strcat(origindir,name,'_Rapports','.xlsx');
%                     movefile (strcat(origindir,'\','Model_rapport_athle.xlsx'),FileName_Rapport);   
    xlswrite(FileName_Rapport,RECAP_sprint,'DATAS','H21');
    xlswrite(FileName_Rapport,cadence_evolution','DATAS','A40');
    xlswrite(FileName_Rapport,Acc_push,'Recap_poussées');

                    imageSize=[482, 326];

                           xlsputimage(FileName_Rapport, imageFileName1, 'Rapport_endurance', 'A53', imageSize)
                                               imageSize=[482, 300];

                    xlsputimage(FileName_Rapport, imageFileName2, 'Rapport_endurance', 'A78', imageSize)


end

Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Recap = RECAP_sprint;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Vitesses = [time VPPP Vittot];
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Distances = [Distdroite Distgauche Disttot];
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).trunk_amp = mvttronc;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).recap_push = Acc_push;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Mvtcadre = cadre_acc;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).trunk_amp_cycles = Cycle_Norm_dos;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Cycle_roue_droite = Cycle_Norm_Rdroite;
Result_athle_piste.(strcat('Sujet_',name)).(strcat('Test_',condition)).Cycle_roue_gauche = Cycle_Norm_Rgauche;
save(strcat(cd,'\','Result_athle_piste.mat'),'Result_athle_piste');





save([cd '/Acc_push_piste_' name '.mat'],'Acc_push');
