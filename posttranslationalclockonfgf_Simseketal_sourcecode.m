clear all;
clc;
set(0,'Defaultlinelinewidth', 2);
set(0, 'DefaultAxesFontWeight', 'normal','DefaultAxesFontSize', 15);
T=304; %% simulation time in mins
delt=0.02;% temporal step size in mins
delx=1;% spatial step size in cell number
N=T/delt; % total no. of iteration steps
xsize=1:40; % cell numbers
tb=6;
L=length(xsize); % length of the tissue in terms of cell no.
tcell=6; % time in mins at which a cell is added in the tissue

rDusp=[1.0000    0.7772    0.4958    0.3118    0.2410    0.1895    0.1712    0.1634]';% experimental data of the Dusp4 mRNA
xr=1:length(rDusp);
fr=fit(xr',rDusp,'exp1');
a=fr.a;
b=fr.b*length(rDusp)/L;
m=a*exp(b*xsize);m=m/max(m);
m(end-L/8:end)=m(end-L/8);

D= 1; % diffusion const
psM=10; % protein synthesis rate of Mek
paC = 0.05; % activation rate of Mek
pdC = 0.5; % degradation rate of Mek
mdL =0.08; % Fgf RNA Decay Rate
psL = 1.5; % protein synthesis rate of Fgf
pdL = 0.4; % protein degradation rate of Fgf
psD6=0.05; % protein synthesis rate of Dusp6
pdD6=0.6; % protein degradation rate of Dusp6
pdD4=0.5; % protein degradation rate of Dusp4
psE=80; % protein synthesis rate of Erk
pdE=0.2; %protein degradation rate of Erk
pdPE=pdE; %protein degradation rate of pErk
pE=0.1; % phosphorylation rate of Erk
dpE=5;% dephosphorylation rate of pErk 
da=1.0; % Association rate of clock and Dusp
dd=0.1; % dissociation rate of clock and Dusp
beta=0.2; % ratio of the degradation rate of clock-bound Dusp to free Dusp
piC=0.1; % feedback coefficient for pMek


TdelML=10; % tranlational time-delay of Fgf
idelML=TdelML/delt;
Tdelact=10; % activation time-delay of pMek
idelact=Tdelact/delt;
TdelC=0; % phophorylation time delay of pErk by pMek
idelC=TdelC/delt;
TdelD6=3; % transcriptional time-delay of Dusp6
idelD6=TdelD6/delt;

f=zeros(1,L); % clock frequency
msL=zeros(1,L); % Fgf mRNA synthesis rate
clock=zeros(L,N);% clock
mLIG=zeros(L,N); % Fgf mRNA
pLIG=zeros(L,N);% Fgfprotein
COMP=zeros(L,N); % Mek
PCOMP=zeros(L,N);% pMek
Erk=zeros(L,N);% Erk
PErk=zeros(L,N);% pErk
Dusp4f=zeros(L,N);% Free Dusp4
Dusp4b=zeros(L,N);% clock-bound Dusp4
Dusp6f=zeros(L,N);% Free Dusp6
Dusp6b=zeros(L,N);% clock-bound Dusp6
phi=zeros(L,N); % clock phase
l=round(T/tcell+5);
PErk_lf=zeros(L+l,N);% pErk in the lab-frame


msL=exp(-xsize/1);msL=msL/max(msL);
msL=3*msL;
msL=0*xsize;
msL(1:tb)=3;
f=1./(1+exp(5*(xsize-35)/L));f=2*pi/30*f; % frequnecy profile of the clock in the tissue
psD4=1*m;% protein synthesis rate of Dusp4
i=1;
count=0;
for t=0:delt:T
    
    %%%%% boundary condition
    diffusion(2:L-1)=D*(pLIG(3:end,i)+pLIG(1:end-2,i)-2*pLIG(2:L-1,i))/(delx^2);
    diffusion(1)=D*(pLIG(2,i)-pLIG(1,i))/(delx^2);% noflux bc
    diffusion(L)=D*(pLIG(L-1,i)-2*pLIG(L,i))/(delx^2);% absorbing bc
    %%%%%%%%%%
    
    if mod(i,tcell/delt)==0
        PErk_lf(end-L+1-count:end-count,1+count*tcell/delt:i)=PErk(:,1+count*tcell/delt:i);
        PErk_mf(:,1+count*tcell/delt:i)=PErk(:,1+count*tcell/delt:i);
        Dusp4_mf(:,1+count*tcell/delt:i)=Dusp4f(:,1+count*tcell/delt:i)+Dusp4b(:,1+count*tcell/delt:i);
        Dusp6_mf(:,1+count*tcell/delt:i)=Dusp6f(:,1+count*tcell/delt:i)+Dusp6b(:,1+count*tcell/delt:i);
        PCOMP_mf(:,1+count*tcell/delt:i)=PCOMP(:,1+count*tcell/delt:i);
        
        phi=circshift(phi,1,1);
        mLIG=circshift(mLIG,1,1);
        pLIG=circshift(pLIG,1,1);
        COMP=circshift(COMP,1,1);
        PCOMP=circshift(PCOMP,1,1);
        Erk=circshift(Erk,1,1);
        PErk=circshift(PErk,1,1);
        
        Dusp4f=circshift(Dusp4f,1,1);
        Dusp4b=circshift(Dusp4b,1,1);
        Dusp6f=circshift(Dusp6f,1,1);
        Dusp6b=circshift(Dusp6b,1,1);
        
        phi(1,:)=phi(2,:);mLIG(1,:)=mLIG(2,:);
        pLIG(1,:)=pLIG(2,:);COMP(1,:)=COMP(2,:);PCOMP(1,:)=PCOMP(2,:);
        Erk(1,:)=Erk(2,:);PErk(1,:)=PErk(2,:);
        Dusp4f(1,:)=Dusp4f(2,:);Dusp4b(1,:)=Dusp4b(2,:);
        Dusp6f(1,:)=Dusp6f(2,:);Dusp6b(1,:)=Dusp6b(2,:);
        
        count=count+1;
    end
    
    Tclock(:,i)=10*(sin(phi(:,i))+1)/2; % total clock 
    phi(:,i+1)=phi(:,i)+delt*f(:);
    Tclock(:,i+1)=10*(sin(phi(:,i+1))+1)/2;
    clock(:,i+1)=(Tclock(:,i+1)-Tclock(:,i))+Tclock(:,i)-Dusp4b(:,i)-Dusp6b(:,i);
    clock(clock<0)=0;
    
    mLIG(:,i+1)=mLIG(:,i)+delt*(msL(:)-mLIG(:,i)*mdL);
    
    if i> idelML && i>idelact
        pLIG(:,i+1)=pLIG(:,i)+delt*(mLIG(:,i-idelML)*psL-pLIG(:,i)*pdL+diffusion(:));
        COMP(:,i+1)=COMP(:,i)+delt*(psM-pLIG(:,i-idelact).*COMP(:,i)*paC-COMP(:,i)*pdC+piC*PCOMP(:,i));
        PCOMP(:,i+1)=PCOMP(:,i)+delt*(pLIG(:,i-idelact).*COMP(:,i)*paC-piC*PCOMP(:,i)-pdC*PCOMP(:,i));
        PCOMP(PCOMP<0)=0;
        COMP(COMP<0)=0;
    else
        pLIG(:,i+1)=pLIG(:,i);
        COMP(:,i+1)=COMP(:,i);
        PCOMP(:,i+1)=PCOMP(:,i);
    end
    
    if i>idelC
        Erk(:,i+1)=Erk(:,i)+delt*(psE+dpE*PErk(:,i).*(Dusp4f(:,i)+Dusp4b(:,i)+Dusp6f(:,i)+Dusp6b(:,i))-pE*Erk(:,i).*PCOMP(:,i-idelC) ...
            -pdE*Erk(:,i));
        PErk(:,i+1)=PErk(:,i)+delt*(-dpE*PErk(:,i).*(Dusp4f(:,i)+Dusp4b(:,i)+Dusp6f(:,i)+Dusp6b(:,i))+pE*Erk(:,i).*PCOMP(:,i-idelC)...
            -pdPE*PErk(:,i));
        PErk(PErk<0)=0;
    end
    
    if i>idelD6
    Dusp6f(:,i+1)=Dusp6f(:,i)+delt*(psD6*PErk(:,i-idelD6)-pdD6*Dusp6f(:,i)-da*clock(:,i).*Dusp6f(:,i)+dd*Dusp6b(:,i));
    Dusp6b(:,i+1)=Dusp6b(:,i)+delt*(-beta*pdD6*Dusp6b(:,i)+da*clock(:,i).*Dusp6f(:,i)-dd*Dusp6b(:,i));
        
    else
        Duspf6f(:,i+1)=Dusp6f(:,i);
        Duspf6b(:,i+1)=Dusp6b(:,i);
    end
    Dusp4f(:,i+1)=Dusp4f(:,i)+delt*(psD4(:)-pdD4*Dusp4f(:,i)-da*clock(:,i).*Dusp4f(:,i)+dd*Dusp4b(:,i));
    Dusp4b(:,i+1)=Dusp4b(:,i)+delt*(-beta*pdD4*Dusp4b(:,i)+da*clock(:,i).*Dusp4f(:,i)-dd*Dusp4b(:,i));

    i=i+1;
end

yP=max(max(PErk_mf));
yD4=max(max(Dusp4f+Dusp4b));
yD6=max(max(Dusp6f+Dusp6b));
yC=max(max(Tclock));
yM=max(max(PCOMP));

% create snapshots from dynamics
%{

mf=6500;
stp=1000;
figure

plot(Tclock(2:end,mf),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1.0,0.4,0.4])
ax = gca;
ax.FontSize = 30;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ylabel('Clock (a.u.)','FontSize',40, 'Color', 'r', 'FontWeight', 'bold');

figure
plot(Tclock(2:end,mf+stp),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1.0,0.1,0.1])
ax = gca;
ax.FontSize = 30;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ylabel('Clock (a.u.)','FontSize',40, 'Color', 'r', 'FontWeight', 'bold');



figure

plot(PErk_mf(2:end,mf),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.1,0.8,0.1],...
    'MarkerFaceColor',[0.4,0.8,0.4])
ax = gca;
ax.FontSize = 30;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ylabel('ppERK (a.u.)','FontSize',40, 'Color', [0.1,0.8,0.1], 'FontWeight', 'bold');
ylim([0 30]);

figure
plot(PErk_mf(2:end,mf+stp),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.1,0.8,0.1],...
    'MarkerFaceColor',[0.1,0.5,0.1])
ax = gca;
ax.FontSize = 30;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ylabel('ppERK (a.u.)','FontSize',40, 'Color', [0.1,0.8,0.1], 'FontWeight', 'bold');
ylim([0 30]);

figure
plot(PErk_mf(2:end,mf+stp),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.1,0.8,0.1],...
    'MarkerFaceColor',[0.1,0.5,0.1])
ax = gca;
ax.FontSize = 30;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ylabel('ppERK (a.u.)','FontSize',40, 'Color', [0.1,0.8,0.1], 'FontWeight', 'bold');
ylim([0 30]);

figure
yyaxis left
plot(Dusp4_mf(2:end,mf),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[1.0,0.1,1.0],...
    'MarkerFaceColor',[1.0,0.7,1.0]);

ax = gca;
ax.YColor= 'k';
ax.FontSize = 30;
ax.FontWeight='bold';
ylabel('Dusp4 (a.u.)','FontSize',40, 'Color', [1.0,0.1,1.0], 'FontWeight', 'bold');
ylim([0 6]);

yyaxis right
plot(Dusp6_mf(2:end,mf),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[1.0,0.45,0.1],...
    'MarkerFaceColor',[1.0,0.75,0.1]);

ax = gca;
ax.FontWeight='bold';
xlabel('PSM (cells)','FontSize',40);
ax.YColor= 'k';
ylabel('Dusp6 (a.u.)','FontSize',40, 'Color', [1.0,0.45,0.1], 'FontWeight', 'bold');
ylim([0 6]);

figure
yyaxis left
plot(Dusp4_mf(2:end,mf+stp),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[1.0,0.1,1.0],...
    'MarkerFaceColor',[0.7,0.3,0.7]);

ax = gca;
ax.YColor= 'k';
ax.FontSize = 30;
ax.FontWeight='bold';
ylabel('Dusp4 (a.u.)','FontSize',40, 'Color', [1.0,0.1,1.0], 'FontWeight', 'bold');
ylim([0 6]);

yyaxis right
plot(Dusp6_mf(2:end,mf+stp),'o','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[1.0,0.45,0.1],...
    'MarkerFaceColor',[1.0,0.25,0.1]);


ax = gca;
ax.FontWeight='bold';
ax.YColor= 'k';
xlabel('PSM (cells)','FontSize',40);

ylabel('Dusp6 (a.u.)','FontSize',40, 'Color', [1.0,0.45,0.1], 'FontWeight', 'bold');
ylim([0 6]);






%{
hold on
plot(Tclock(2:end,mf+stp))
hold on
plot(Tclock(2:end,mf+2*stp))
hold on
plot(Tclock(2:end,mf+3*stp))


figure

plot(Dusp4_mf(2:end,mf))
hold on
plot(Dusp4_mf(2:end,mf+stp))
hold on
plot(Dusp4_mf(2:end,mf+2*stp))
hold on
plot(Dusp4_mf(2:end,mf+3*stp))



figure

plot(PErk_mf(2:end,mf))
hold on
plot(PErk_mf(2:end,mf+stp))
hold on
plot(PErk_mf(2:end,mf+2*stp))
hold on
plot(PErk_mf(2:end,mf+3*stp))



figure

plot(Tclock(2:end,mf))
hold on
plot(Tclock(2:end,mf+stp))
hold on
plot(Tclock(2:end,mf+2*stp))
hold on
plot(Tclock(2:end,mf+3*stp))

%}
%}

% create the video writer with 1 fps
  writerObj = VideoWriter('plotmovies2.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);

for i=100/delt:1/delt:300/delt
       
    time=sprintf('time=%.2f mins',i*delt);
    subplot(3,1,1);
    plot(Tclock(2:end,i),'m');
    ylim([0 yC]);xlabel('cell #','fontsize',20); xlim([-1 length(xsize)]);ylabel('Clock','fontsize',15);
    title(time,'fontsize',10);

    subplot(3,1,2);
    plot(Dusp4_mf(2:end,i));hold on; plot(Dusp6_mf(2:end,i)); ylim([0 yD4]);legend('Dusp4','Dusp6','fontsize',10);legend box off;
    ylabel('Dusp','fontsize',15);xlabel('cell #','fontsize',20);xlim([-1 length(xsize)]);
    
    hold off;
    
    subplot(3,1,3);
    plot(PErk_mf(2:end,i),'g'); 
    ylim([0 yP]);
    legend('ppERK','fontsize',10);legend box off;
    ylabel('ppERK','fontsize',15);xlabel('cell #','fontsize',20);xlim([-1 length(xsize)]);
    
    %%pause(.1);
    frame = getframe(gcf);
    drawnow
    writeVideo(writerObj, frame);
end
close(writerObj);
