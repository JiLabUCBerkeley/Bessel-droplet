%% Function Discritption
% 1. flag_saveTemResults=1? save the results during the simulation process.
% 2. flag_GaussianBeam=1; Genrally use the GaussianBeam profile.
% 3. flag_3D=1, save a 3D PSF.
% 4. flag_NormalizePSF=1;
% 5. flag_mask='DR', set the simulation for doule rings.
% 6. dPhase=0 to generate equal phase values within the two rings;
% dPhase=pi to generate pi phase difference between the two rings.
% 7. Use the optimized R2 for coressponding R1 from the
% Bessel_droplet_SIM.m program.

%% Add dependecy path
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

clc;
clear all;
ClearAllFigures;

%% Plot parameters:
lineW=2;
titleSize=10;
labelSize=14;
%% Step 1 . Start the main program “Bessel_droplet_SIM_dPhase.m”
%% Fixed Simulation inputs Pars
  wavelength=0.94;%um, the two photon excitation wavelength;
  K=2*pi/(wavelength*1e-3);% in mm
  BeamDiameter=7.2; %mm, Beam size at SLM1, SLM is of size 7.68mm by 7.68mm
  obj.tubeLength=180;% mm, Olympus 180, Nikon 200, Zeiss 165
  obj.NA=1.05;% NA of objective
  obj.magnification=25;% magnification factor of objective
  obj.refind=1.333;% refractive index of immersion medium  
  
%% ############################################ Simulation Control Bools:
flag_saveTemResults=1;
flag_GaussianBeam=1;
flag_3D=0;
flag_NormalizePSF=1;
flag_mask='DR';
%% ############################################ Phase difference 
dPhase=pi;
%% ############################################ Ring Mask Pars
mag_tem=(350/750)*10;%Obj pupil to Bessel mask magnification

%% Set R1, R2
R1=2.88;%NA=0.4
R2=1.44;%NA=0.4

% R1=3.6;%NA=0.5
% R2=2.088;%NA=0.5

% R1=4.32;%NA=0.6
% R2=2.506;%NA=0.6

% R1=5.04;%NA=0.7
% R2=2.9232;%NA=0.7

% R1=5.76;%NA=0.8
% R2=3.456;%NA=0.8
RingThickness1=0.300;
RingThickness2=0.300;

%% ########### Set result folder based on pars:
dateStr=datestr(now, 'mmddyy-HHMM');
currentFolder = pwd;
resultPath=[currentFolder,'\results\'];
    if exist(resultPath,'dir')
    else
        mkdir(resultPath)
    end


%% ###########    Cal Beam Size:
beamD=BeamDiameter; %mm, 
beadD_Pupil=beamD*(350/750)*10;%mm, 
tubeLength=obj.tubeLength;
NA=obj.NA;
magnification=obj.magnification;
refind=obj.refind;
fOB=tubeLength/magnification; % 
rOB=fOB*NA; %mm


%% Start interation tests:
TrialNum=1;
pk_Value=zeros(1,TrialNum)';
FWHM2=zeros(1,TrialNum)';
FWHMz=zeros(1,TrialNum)';

for tt=1:1:TrialNum
outerD_B=(R2+RingThickness2/2)*2;
innerD_B=(R2-RingThickness2/2)*2;
outerD_A=(R1+RingThickness1/2)*2;
innerD_A=(R1-RingThickness1/2)*2;

%% Step 2. Initialize an optical field (field_Pupil_ring) for the Bessel SLM located at the pupil plane
step_Pupil=3;%um
r_Pupil=0:step_Pupil:rOB*1000;%um
% 
if(flag_GaussianBeam)
Amplitude_Pupil=1*exp(-2*r_Pupil.^2/(beadD_Pupil/2*1000)^2); 
else
Amplitude_Pupil=1*ones(size(r_Pupil)); 
end

%% Apply mask to the optical filed at the back pupil plane
step_mask=3;% micron;
r_mask=0:step_Pupil:rOB*1000; % radial coordinate at mask.
flag_maskA=(r_mask>=innerD_A*1000/2)&(r_mask<=outerD_A*1000/2);
flag_maskB=(r_mask>=innerD_B*1000/2)&(r_mask<=outerD_B*1000/2);
flag_DR=flag_maskA|flag_maskB;
% 
switch flag_mask
    case 'A'
maskProfile=double([fliplr(double(flag_maskA)) flag_maskA]);
Amplitude_Pupil(~flag_maskA)=0; 
    case 'B'
maskProfile=double([fliplr(double(flag_maskB)) flag_maskB]);
Amplitude_Pupil(~flag_maskB)=0; 
    case 'DR'
maskProfile=double([fliplr(double(flag_DR)) flag_DR]);
Amplitude_Pupil(~flag_DR)=0;
end
PhasePhase_Pupil=zeros(size(Amplitude_Pupil));%rad
PhasePhase_Pupil(flag_maskB)=mean(PhasePhase_Pupil(flag_maskA))+dPhase(tt);
field_Pupil=Amplitude_Pupil.*exp(1i*PhasePhase_Pupil);
field_Pupil_ring=field_Pupil;

%% Display masked field at the pupil plane
field_Pupil2D_Amp=Creat3DPSF(r_Pupil,abs(field_Pupil_ring));
field_Pupil2D_Phase=Creat3DPSF(r_Pupil,unwrap(angle(field_Pupil_ring)));%phase in unit pi
field_PupilSection_Amp=[rot90(field_Pupil2D_Amp,2)  rot90(field_Pupil2D_Amp,1); rot90(field_Pupil2D_Amp,3),field_Pupil2D_Amp];
field_PupilSection_Amp(isnan(field_PupilSection_Amp))=0;
field_PupilSection_Phase=[rot90(field_Pupil2D_Phase,2)  rot90(field_Pupil2D_Phase,1); rot90(field_Pupil2D_Phase,3),field_Pupil2D_Phase];
field_PupilSection_Phase(isnan(field_PupilSection_Phase))=0;
Xaxis=[-flip(r_Pupil,2),r_Pupil];

%% Plot asked field at the pupil plane

fig51=figure(51);
ax1=subplot(2,2,1);
imagesc(Xaxis,Xaxis,field_PupilSection_Amp);
h_title=title('Optical Field Amplitude Distribution at Pupil Plane');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
colormap(ax1,jet)
colorbar;

ax2=subplot(2,2,3);
imagesc(Xaxis,Xaxis,field_PupilSection_Phase);
h_title=title('Optical Field Phase Distribution at Pupil Plane (unwrapped)');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
colormap(ax1,jet)
colorbar;

ax3=subplot(2,2,2);
plot(Xaxis,field_PupilSection_Amp(end/2,:),'k','linewidth',1.5);
h_title=title('Optical Field  Amplitude Profile Distribution at Pupil Plane');h_title.FontSize=titleSize;
xlabel('x (um)');
xlim([-10000 10000]);
ylabel('Optical field amplitude (a.u.)');
hold on;
plot([fliplr(-r_mask) r_mask],maskProfile.*max(max(field_PupilSection_Amp)),'-r');
hold off;

ax4=subplot(2,2,4);
plot(Xaxis,field_PupilSection_Phase(end/2,:),'k','linewidth',1.5);
h_title=title('Optical Field Phase Profile at Pupil Plane (unwrapped)');h_title.FontSize=titleSize;
xlabel('x (mm)');
xlim([-10000 10000]);
ylabel('Optical field phasee (a.u.)');
hold on;
plot([fliplr(-r_mask) r_mask],maskProfile.*min(min(field_PupilSection_Phase)),'-r');
hold off;
if(flag_saveTemResults)
filename=[resultPath, '\PupilFiled_masked_',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f')];
saveas(fig51,[filename,'.tif']);
end

%% Step 3. Fourier transform field_Pupil_ring to the SLM conjugated to the objective focal plane
% propogate the filed at obj pupil to mask and do Circular Fourier
% transform
r_PupilMask=r_Pupil./mag_tem;
field_PupilMask=field_Pupil_ring;
step_SLM=2.5;% um
pixel_SLM=15;% um
r_SLM=0:step_SLM:pixel_SLM*512*ceil(sqrt(2));%um, test field after 50mm diameter lens, coordinates extends to 512*sqrt(2)pixels.
f1=200;%mm

field_PupilMask_Cal=field_PupilMask;
r_PupilMask_Cal=r_PupilMask;
Field_SLM=Fourier_CircularLens(field_PupilMask_Cal,wavelength,r_PupilMask_Cal,f1,r_SLM);%r, lambda-um, f-mm

field_Mask2D_Amp=Creat3DPSF(r_SLM,abs(Field_SLM));
field_Mask2D_Phase=Creat3DPSF(r_SLM,angle(Field_SLM));%phase in unit pi
field_MaskSection_Amp=[rot90(field_Mask2D_Amp,2)  rot90(field_Mask2D_Amp,1); rot90(field_Mask2D_Amp,3),field_Mask2D_Amp];
field_MaskSection_Phase=[rot90(field_Mask2D_Phase,2)  rot90(field_Mask2D_Phase,1); rot90(field_Mask2D_Phase,3),field_Mask2D_Phase];
Xaxis=[-flip(r_SLM,2),r_SLM];

SLM1_Xaxis=Xaxis(end/2-3840/step_SLM:end/2+3840/step_SLM);
SLM1_Phase=field_MaskSection_Phase(end/2-3840/step_SLM:end/2+3840/step_SLM,end/2-3840/step_SLM:end/2+3840/step_SLM);
SLM1_Amp=field_MaskSection_Amp(end/2-3840/step_SLM:end/2+3840/step_SLM,end/2-3840/step_SLM:end/2+3840/step_SLM);
%% Display phase pattern on SLM1:

fig52=figure(52);
ax1=subplot(2,2,1);
imagesc(SLM1_Xaxis,SLM1_Xaxis,SLM1_Amp);
h_title=title('Optical Field Amplitude Distribution at SLM');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
xlim([-3840 3840]);
ylim([-3840 3840]);
colormap(ax1,jet)
colorbar;

ax2=subplot(2,2,2);
imagesc(SLM1_Xaxis,SLM1_Xaxis,SLM1_Phase);
h_title=title('Optical Field Phase Distribution at SLM');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
xlim([-3840 3840]);
ylim([-3840 3840]);
colormap(ax2,jet)
colorbar;

ax3=subplot(2,2,3);
plot(SLM1_Xaxis,SLM1_Amp(round(end/2),:));
h_title=title('Optical Field Amplitude Profile at SLM');h_title.FontSize=titleSize;
xlabel('x (um)');
xlim([-3840 3840]);

ax4=subplot(2,2,4);
plot(SLM1_Xaxis,SLM1_Phase(round(end/2),:));
hold on;
plot(SLM1_Xaxis(1:6:end),SLM1_Phase(round(end/2),1:6:end),'.r');
hold off;
h_title=title('Optical Field Amplitude Profile at SLM');h_title.FontSize=titleSize;
xlabel('x (um)');
xlim([-3840 3840]);

% Save Simulation figure
if(flag_saveTemResults)
filename=[resultPath, '\SLM_PhasePattern_',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f')];
saveas(fig52,[filename,'.tif']);
end
% Save mask file for experiemnt
PasePattenSLM1=field_MaskSection_Phase(end/2-3840/step_SLM:end/2+3840/step_SLM,end/2-3840/step_SLM:end/2+3840/step_SLM);
tem_A=PasePattenSLM1;
tem_A(tem_A(:)<-3)=-pi;
tem_A(tem_A(:)>3)=pi;
b_select=(tem_A(:)>-1)&(tem_A(:)<1);
tem_A(b_select)=0;
tem_B=imresize(tem_A,[512 512],'nearest');
tem_C=tem_B-min(min(tem_B));
PasePattenSLM1_512=uint8(tem_C./pi*256/2);

filename=[resultPath, '\SLM_PhaseMask_',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f')];
imwrite(PasePattenSLM1_512,[filename '.bmp']);


%% Step 4. Validation: Propagate the phase patterns on the SLM to the mask (the mask is conjugated to the objective pupil plane) 
Field_SLM_Amp=1*exp(-2*r_SLM.^2/(BeamDiameter/2*1000)^2);
Field_SLM_Phase=angle(Field_SLM);
Field_SLM_Phase=circshift(Field_SLM_Phase,10);
Field_SLM_New=Field_SLM_Amp.*exp(1i*Field_SLM_Phase);

field_PupilMaskNew=Fourier_CircularLens(Field_SLM_New,wavelength,r_SLM,f1,r_PupilMask);
field_PupilMaskNew2D_Amp=Creat3DPSF(r_PupilMask,abs(field_PupilMaskNew));
field_PupilMaskNew2D_Phase=Creat3DPSF(r_PupilMask,angle(field_PupilMaskNew));%phase in unit pi
field_PupilMaskNewSection_Amp=[rot90(field_PupilMaskNew2D_Amp,2)  rot90(field_PupilMaskNew2D_Amp,1); rot90(field_PupilMaskNew2D_Amp,3),field_PupilMaskNew2D_Amp];
field_PupilMaskNewSection_Phase=[rot90(field_PupilMaskNew2D_Phase,2)  rot90(field_PupilMaskNew2D_Phase,1); rot90(field_PupilMaskNew2D_Phase,3),field_PupilMaskNew2D_Phase];
X_PupilMask=[-flip(r_PupilMask,2),r_PupilMask];



fig53=figure(53);
ax1=subplot(2,2,1);
imagesc(X_PupilMask,X_PupilMask,field_PupilMaskNewSection_Amp);
h_title=title('Optical Field Amplitude Distribution at mask');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');

colormap(ax1,jet)
colorbar;

ax2=subplot(2,2,2);
imagesc(X_PupilMask,X_PupilMask,field_PupilMaskNewSection_Phase);
h_title=title('Optical Field Phase Distribution at mask');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');

colormap(ax2,jet)
colorbar;

ax3=subplot(2,2,3);
plot(X_PupilMask,field_PupilMaskNewSection_Amp(round(end/2),:));
hold on;
plot([fliplr(-r_mask) r_mask]/mag_tem,maskProfile.*max(field_PupilMaskNewSection_Amp(round(end/2),:)),'-r');
hold off;
h_title=title('Optical Field Amplitude Profile at mask');h_title.FontSize=titleSize;
xlabel('x (um)');


ax4=subplot(2,2,4);
plot(X_PupilMask,field_PupilMaskNewSection_Phase(round(end/2),:));
hold on;
plot([fliplr(-r_mask) r_mask]/mag_tem,maskProfile.*3.5,'-r');
hold off;
h_title=title('Optical Field Amplitude Profile at mask');h_title.FontSize=titleSize;
xlabel('x (um)');

%% Save confirmation results
if(flag_saveTemResults)
filename=[resultPath, '\PupilConfirm_',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f')];
saveas(fig53,[filename,'.tif']);
end

%% Generate the PSF using Richard Wolf PSF model for high NA objs
field_Pupil_PSFCal=field_Pupil_ring;
x=0;
y=0;
z=-300:1:300; %um, range of z of the PSF relative to the focal plane

field_Pupil_Amp_nor=abs(field_Pupil_PSFCal)./max(abs(field_Pupil_PSFCal));
field_Pupil_Phase_nor=angle(field_Pupil_PSFCal);


r_PupilObj=r_Pupil;
field_PupilObj=field_Pupil_Amp_nor.*exp(1i*field_Pupil_Phase_nor);
flag1=r_PupilObj<=rOB*1000;%coordinates that are inside the objective

r_angle=asin(r_PupilObj/fOB/refind/1000);
r_angle(~flag1)=[];
field_PupilObj(~flag1)=[];
[PSF] = Calc_Annular_Field_Integrals(x, y, z,field_PupilObj,wavelength,r_angle,refind);

PSF=squeeze(PSF);
xlen=length(x); % coodinates(0,0,0)
ylen=length(y);
zlen=length(z);

xyzlen=sort([xlen,ylen,zlen]);
output.maxPSF=max(PSF);
PSFz=PSF;
PSF2=PSF/max(PSF);% Normalized axial PSF
drawnow
%% Plot Figure 56 (a)
fig56=figure(56);
subplot(2,1,1);
plot(z,PSF2,'b','LineWidth',lineW);
text(-250,1,['R1: ', num2str(R1)]);
text(-250,0.8,['d-R1: ', num2str(RingThickness1)]);
text(-250,0.6,['dPhase: ', num2str(dPhase(tt))]);
text(-250,0.4,['d-R2: ', num2str(RingThickness2)]);


box off;
set(gca,'FontName','Arial')
set(gca,'color','none')
set(gca,'TickDir','out');
set(gcf,'PaperPositionMode','auto');
set(gca,'yTick',[0:0.5:1]);     
% ylim([0, 1.2]);
xlabel('z(\mum)','Fontname','Arial','Fontsize',labelSize);
ylabel('Axial PSF','Fontname','Arial','Fontsize',labelSize)
xlim([z(1),z(end)]);
flag=PSF2>=0.5;
II=find(flag);     
FWHMz(tt)=z(II(end))-z(II(1));   
h_title=title(['Axial PSF mm;FWHM=',num2str(FWHMz(tt)),'\mum']);h_title.FontSize=titleSize;
drawnow; 

% calculating lateral FWHM along y
    [~,I]=max(PSF2);
    output.PSFz=[z(:),PSF2(:)];
    z=0;
    x=0;
    y=-5:.02:5;
    [PSFy] = Calc_Annular_Field_Integrals(x, y, z,field_PupilObj,wavelength,r_angle,refind);
    PSFy=squeeze(PSFy);
    PSFy2=PSFy/max(PSFy);
    y2=y;y2(y2==0)=eps;
    w=2*pi*obj.NA/wavelength;
% do fitting to get more accurate lateral FWHM
    fun=@(par,z,r)(z-(par(1)*besselj(1,par(2)*r)./r/2/pi).^4);
    fun2=@(par,r)((par(1)*besselj(1,par(2)*r)./r/2/pi).^4);
    par0=[1,1]*w;
    par=lsqnonlin(fun,par0,[],[],[],PSFy2(:),y2(:));   
    output.PSFy=[y(:),PSFy2(:)];
    step=0.001;
    y_finner=0:step:2;y_finner(1)=eps;
    y_finner=[-y_finner(end:-1:1),y_finner(2:end)];
    PSFy_fit=fun2(par,y_finner);
    flag=PSFy_fit>=0.5;
    FWHM2(tt)=sum(flag)*(y_finner(2)-y_finner(1))*1;  
    output.PSFy_fit=[y_finner(:),y_finner(:)];
    output.FWHMy_fit=FWHM2(tt);
    drawnow
    
%% Plot Figure 52

    [pks,locs] = findpeaks(PSFy2);
    [sortedpks, pksInds] = sort(pks,'descend');
    pk_Value(tt)=sortedpks(2);
    pk_Position=(locs(pksInds(2))-length(y)/2)*mean(diff(y));
    
    subplot(2,1,2);
    plot(y,PSFy2,'b-','LineWidth',lineW);
    hold on;
    plot(pk_Position,pk_Value(tt),'ks','markerfacecolor',[0 0 0]);
    text(pk_Position,pk_Value(tt),['Max: ', num2str(pk_Value(tt))]);
    hold off;
    
    h_title=title(['PSF along y, FWHM=',num2str(FWHM2(tt),'%.3f'),'\mum']);h_title.FontSize=titleSize;
    box off;
    set(gca,'FontName','Arial')
    set(gca,'color','none')
    set(gca,'TickDir','out');
    set(gca,'yTick',[0:0.5:1]);
    set(gca, 'YScale', 'log');
    ylim([0, 1.2])
    xlabel('y(\mum)','Fontname','Arial','Fontsize',labelSize);
    ylabel('Vertical PSF','Fontname','Arial','Fontsize',labelSize)
    xlim([y(1),y(end)]);
     
filename=[resultPath, '\PSF_Profiles_R2_',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f'),'_Z=',num2str(z)];
saveas(fig56,[filename,'.tif']);
 
%% Step 5. (optional) Calculate the 3D PSF using Richard Wolf PSF model for high NA objective.
if(flag_3D)
ySimGrids=-5:0.1:5;
zSimGrids=-100:1:100; %um, range of z of the PSF relative to the focal plane
    
PSF_Stack=zeros(length(zSimGrids),length(ySimGrids));
PSF_Cross_Stacks=zeros(length(zSimGrids),length(ySimGrids),length(ySimGrids));

parfor ii=1:1:length(zSimGrids)
    zPosition=zSimGrids(ii);
    yPosition=ySimGrids;
    xPosition=0;
    [PSFy] = Calc_Annular_Field_Integrals(xPosition, yPosition, zPosition,field_PupilObj,wavelength,r_angle,refind);
    PSFy=squeeze(PSFy);
    %Not simulate on the focal plane, not fitting
    PSF_Stack(ii,:)=PSFy;
    PSF_Cross_Stacks(ii,:,:)=Creat3DPSF(yPosition,PSFy);
end
%% Transfer to tiff image, Normalized to max value
[Max_3DPSF, I]=max(PSF_Cross_Stacks(:));
filename=[resultPath, '\PSF_3D',flag_mask,'_dPhase_',num2str(dPhase(tt),'%.3f'),'_R1_',num2str(R1,'%.3f'),'_dR2_',num2str(RingThickness2,'%.3f'),'_dR1_',num2str(RingThickness1,'%.3f')];
for ii=1:1:length(zSimGrids)
    
     psf_img=uint16(squeeze(PSF_Cross_Stacks(ii,:,:))./Max_3DPSF*65536);
     psf_img_Square=uint16((squeeze(PSF_Cross_Stacks(ii,:,:))./Max_3DPSF).^2*65536);
     
     if (ii==1)
        imwrite(uint16(psf_img), [filename,'.tif'], 'tif', 'WriteMode', 'overwrite','Compression', 'none');
     else
         imwrite(uint16(psf_img), [filename,'.tif'], 'tif', 'WriteMode', 'append','Compression', 'none');
     end 
end
end

end

