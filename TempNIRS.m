

%% Load data
clear all  
folder='C:Data'
d=dir(folder)
subf=d([d.isdir])
im=[]
for k=3:numel(subf)
  subfolder=subf(k).name
  subf1=fullfile(folder,subfolder)
  f=dir([subf1 '\*NIRSCYRIL*.mat'])
  for ii=1:numel(f)
      file=fullfile(subf1,f(ii).name);
      importfile(file);
  end
   f=dir([subf1 '\*Temperature*.mat'])
  for ii=1:numel(f)
      file=fullfile(subf1,f(ii).name);
      importfile(file);
  end

  
%% Load spectra
load spectra_hhb.mat
load spectra_o2hb.mat
load spectra_water.mat
load spectra_temp.mat

%% Wavelength selection
load wavelengths.mat
wavelengths_selected=wavelengths;minwave=830;maxwave=850;

Spectra_water=vq_water(wavelengths_selected(:,1)>minwave&wavelengths_selected(:,1)<maxwave,:);
Spectra_o2hb=vq_o2hb(wavelengths_selected(:,1)>minwave&wavelengths_selected(:,1)<maxwave,:);
Spectra_hhb=vq_hhb(wavelengths_selected(:,1)>minwave&wavelengths_selected(:,1)<maxwave,:);
Spectra_temp=vq_temp(wavelengths_selected(:,1)>minwave&wavelengths_selected(:,1)<maxwave,:);
wavelengths=wavelengths(wavelengths_selected (:,1)>minwave&wavelengths_selected (:,1)<maxwave,:);

%% Savitzsky-Golay filter
B = sortrows([Temperature,Abs_Attenuation'],1) ;B=B';
for ii=1:size(B,2)
    y = sgolayfilt(B(2:end,ii),1,81); x(:,ii)=y;
end
Measured_Temperature=B(1,:);
x=x(wavelengths_selected (:,1)>minwave&wavelengths_selected (:,1)<maxwave,:);

%% Concentration calculation
x=x'+1;ppf=60; dpf=6; L=3;
for ll=1:size(x,2)
    dOD_w(:,ll)= -log(x(:,ll)/mean(x(:,ll)));
end
for ll=1:size(x,2)
dOD_L (:,ll)= dOD_w(:,ll)* ppf/(dpf*L);        %Pre-divide by pathlength (dpf,pvcetc)
end
E=[Spectra_o2hb Spectra_hhb Spectra_water Spectra_temp];
Einv = inv(E'*E)*E';                            %Linear inversion operator (for 2 or more wavelengths)
REG = Einv * dOD_L';                            %Solve for HbO and HbR (This is the least-squares solution for unlimited #of wavelengths)
HbO = REG(1,:);HbR = REG(2,:);Water = REG(3,:);Temp = REG(4,:);Temp=fliplr(Temp);

%% Calibration
i=1;
for delta=10:0.1:60
    X(i,:) = Spectra_water +(delta*Spectra_temp);
    i=i+1;
end
Hollis_Temperature = [10:0.1:60];


%% Select number of PCs
[Yp,Predicted_Temperature,b,T,P,U,R,R2X,R2Y]=plspredict(X,Hollis_Temperature',1,5,x)
Predicted_Temperature=((Predicted_Temperature(:,1)-nanmean(Predicted_Temperature(:,1))).*0.85)+(Measured_Temperature(1));    
%% Plot data

% mdl=fitlm(Measured_Temperature',Predicted_Temperature)
% plot(mdl)
% hold on
% plot(Measured_Temperature,Measured_Temperature,'Color','r','LineWidth',2)
% plot(Predicted_Temperature(:,1))
% hold on
% plot(Measured_Temperature)

IDx= str2num(subf1(end-2:end))
ID=1:length(Measured_Temperature);ID(:)=IDx+0;ID=ID';

%% Assess prediction
mdl1=fitlm(Measured_Temperature',[HbO+HbR]');  
dlmwrite('Table_PCA1_thb_r2.csv',mdl1.Rsquared.Ordinary,'-append')
mdl2=fitlm(Measured_Temperature',Water');  
dlmwrite('Table_PCA1_water_r2.csv',mdl2.Rsquared.Ordinary,'-append')
mdl3=fitlm(Measured_Temperature',Predicted_Temperature);   
dlmwrite('Table_PCA1_temp.csv',[mdl3.Rsquared.Ordinary,r(1,2)],'-append')

clearvars -except subf subf1 file d folder
end
%%



