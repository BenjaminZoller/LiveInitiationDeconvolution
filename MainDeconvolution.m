%   Copyright (c) 2024, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

%% To generate the input data structure for the deconvolution
clear
clc
close all

%NC = 12;
%NC = 13;
NC = 14;

mypath = '.';
datafolder = [mypath,filesep,'InData'];
resfolder = [mypath,filesep,'OutDeconvRes'];

GList = {'hb_MS2_INTRON','hb_MS2_3UTR','kr_MS2_INTRON','kr_MS2_3UTR','kni_MS2_INTRON','gt_MS2_INTRON_Female_Ant',...
    'gt_MS2_INTRON_Female_Post','gt_MS2_INTRON_Male_Ant','gt_MS2_INTRON_Male_Post',...
    'kni_MS2_INTRON_hb_null','kni_MS2_INTRON_dist_del','hb_MS2_INTRON_dist_del'};

ll = 1;
line = GList{ll};

% Time window for deconvolution
% Deconvolution is only performed when there is transcription,
% namely ~3 min after the onset of each nuclear cycle.
% We stop deconvolution ~4 min before mitosis, 
% as we cannot model the shutting down of transcription, 
% where Pol2 drop off the DNA.
switch NC
    case 12
        NDtresh = 20;
        Vt = 16:54;
    case 13
        Vt = 16:91;
        NDtresh = 100;
    case 14
        Vt = 16:301;
        NDtresh = 100;
end

OutDS = importdata([datafolder,filesep,line,'_NC',num2str(NC),'_OutDS.mat']);

% select AP where there are sufficient data
M1 = vertcat(OutDS.M1);
Nd = arrayfun(@(x) length(x.Traces),OutDS);
Iap = mean(M1(:,Vt),2) > 0.5 & Nd(:) > NDtresh;
Vap = 1:length(OutDS);
Vap = Vap(Iap);

% setup variables
T = OutDS(1).Ti;
AP = [OutDS(Vap).AP];
Nap = length(Vap);
T = T(Vt);
t1 = T(1);
t2 = T(end);
T = OutDS(1).Ti;

% generate data
for i=1:Nap
    Data(i).AP = AP(i);
    Data(i).T = T;
    It = false(size(T));
    It(Vt) = true;
    Data(i).It = It;
    Data(i).M1 = OutDS(Vap(i)).M1; %the mean signal is important for the prior
    
    Traces = OutDS(Vap(i)).Traces;
              
    kb = arrayfun(@(x) find(x.t>=t1,1,'first'),Traces);
    ke = arrayfun(@(x) find(x.t<=t2,1,'last'),Traces);
    emb = [Traces.emb];
    nid = [Traces.nid];
    
    for j=1:length(Traces)
        t = Traces(j).t;
        int = Traces(j).int;
        ap = Traces(j).ap;
        dv = Traces(j).dv;
        rs = Traces(j).rs;

        t = t(kb(j):ke(j));
        int = int(kb(j):ke(j));
        ap = ap(kb(j):ke(j));
        dv = dv(kb(j):ke(j));
        rs = rs(:,kb(j):ke(j));
                
        Data(i).Traces(j).t = t;
        Data(i).Traces(j).int = int;
        Data(i).Traces(j).ap = ap;
        Data(i).Traces(j).dv = dv;
        Data(i).Traces(j).rs = rs;
        Data(i).Traces(j).emb = emb(j);
        Data(i).Traces(j).nid = nid(j);
    end
end

save([resfolder,filesep,'Data_',line,'_Nc',num2str(NC),'.mat'],'Data','-v7.3');

%% Sampling Deconvolution
clear
clc
close all

%NC = 12;
%NC = 13;
NC = 14;

GList = {'hb_MS2_INTRON','hb_MS2_3UTR','kr_MS2_INTRON','kr_MS2_3UTR','kni_MS2_INTRON','gt_MS2_INTRON_Female_Ant',...
    'gt_MS2_INTRON_Female_Post','gt_MS2_INTRON_Male_Ant','gt_MS2_INTRON_Male_Post',...
    'kni_MS2_INTRON_hb_null','kni_MS2_INTRON_dist_del','hb_MS2_INTRON_dist_del'};

ll = 1;
line = GList{ll};

mypath = '.';
datafolder = [mypath,filesep,'InData'];
resfolder = [mypath,filesep,'OutDeconvRes'];

% Kernel built based on construct
% Kernel.X documents the position of each of the 24 ms2 loop and
% the cleavage positon (end of 3'UTR or ideally polyA site)
% Kernel.K provides the normalized signal contribution for 24 loops
Kernel = importdata([datafolder,filesep,'Kernel_',line,'.mat']);
% Our assessement of the 2-photon measurement noise
Noise = importdata([datafolder,filesep,'MeasNoise2PhotonMS2.mat']);
% Elongation rate in bp/min that we measured using dual color lines
Ke = 1800; 

load([resfolder,filesep,'Data_',line,'_Nc',num2str(NC),'.mat'],'Data');

Nap = length(Data);
Vap = 1:Nap;
%Vap = 16; % Only ap=16 is provided as example data

for ap=Vap
    disp(['ap bin #',num2str(ap)])
    Traces = Data(ap).Traces;
    M1 = Data(ap).M1;
    M1(isnan(M1)) = 0;
    T = Data(ap).T;
    It = Data(ap).It;
    
    [OutTr,OutPr] = DeconvSampling(Traces,M1,T,It,Kernel,Noise,Ke);
        
    save([resfolder,filesep,'DECONV_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutTr','-v7.3');
    save([resfolder,filesep,'BPHENO_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutPr','-v7.3');
end

%% Compute Burst Pheno (clustering)
clear
clc
close all

%NC = 12;
%NC = 13;
NC = 14;

GList = {'hb_MS2_INTRON','hb_MS2_3UTR','kr_MS2_INTRON','kr_MS2_3UTR','kni_MS2_INTRON','gt_MS2_INTRON_Female_Ant',...
    'gt_MS2_INTRON_Female_Post','gt_MS2_INTRON_Male_Ant','gt_MS2_INTRON_Male_Post',...
    'kni_MS2_INTRON_hb_null','kni_MS2_INTRON_dist_del','hb_MS2_INTRON_dist_del'};

ll = 1;
line = GList{ll};

mypath = '.';
resfolder = [mypath,filesep,'OutDeconvRes'];

Ns = 1000;
% Parameters for burst calling
Ws = 5/6; %length of the averaging window in min
Ks = 2; %sensitivity threshold in #of ini events

load([resfolder,filesep,'Data_',line,'_Nc',num2str(NC),'.mat'],'Data');
Nap = length(Data);
Vap = 1:Nap;
%Vap = 16; % Only ap=16 is provided as example data

for ap=Vap
    disp(['ap bin #',num2str(ap)])
    load([resfolder,filesep,'DECONV_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutTr');
    if exist([resfolder,filesep,'BPHENO_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'file')
        load([resfolder,filesep,'BPHENO_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutPr');
    end
    
    T = Data(ap).T;
    It = Data(ap).It;
  
    OutPheno = ExtractIniPheno(OutTr,T(It),Ns,Ws,Ks);
    
    OutPr.G1 = OutPheno.G1;
    OutPr.G2 = OutPheno.G2;
    OutPr.NG1 = OutPheno.NG1;
    OutPr.NG2 = OutPheno.NG2;
    OutPr.A1 = OutPheno.A1;
    OutPr.A2 = OutPheno.A2;
    OutPr.NA1 = OutPheno.NA1;
    OutPr.NA2 = OutPheno.NA2;
    OutPr.Gcorr = OutPheno.Gcorr;
    
    % burst pheno
    OutPr.pn = OutPheno.pn;
    OutPr.tn = OutPheno.tn;
    OutPr.bn = OutPheno.bn;
    OutPr.fn = OutPheno.fn;
    OutPr.kn = OutPheno.kn;
    OutPr.rn = OutPheno.rn;
    OutPr.knb = OutPheno.knb;
    OutPr.ttn = OutPheno.ttn;
    OutPr.tcn = OutPheno.tcn;
    
    OutPr.epn = OutPheno.epn;
    OutPr.etn = OutPheno.etn;
    OutPr.ebn = OutPheno.ebn;
    OutPr.efn = OutPheno.efn;
    OutPr.ekn = OutPheno.ekn;
    OutPr.ern = OutPheno.ern;
    OutPr.eknb = OutPheno.eknb;
    OutPr.ettn = OutPheno.ettn;
    OutPr.etcn = OutPheno.etcn;
    
    OutPr.zn = OutPheno.zn;
    
    % distrib
    OutPr.wf = OutPheno.wf;
    OutPr.wt = OutPheno.wt;
    OutPr.wb = OutPheno.wb;
    OutPr.cf = OutPheno.cf;
    OutPr.ct = OutPheno.ct;
    OutPr.cb = OutPheno.cb;
    OutPr.ecf = OutPheno.ecf;
    OutPr.ect = OutPheno.ect;
    OutPr.ecb = OutPheno.ecb;
    
    OutPr.vt = OutPheno.vt;
    OutPr.vb = OutPheno.vb;
    
    % burst pheno smoothed out
    OutPr.spn = OutPheno.spn;
    OutPr.stn = OutPheno.stn;
    OutPr.sbn = OutPheno.sbn;
    OutPr.sfn = OutPheno.sfn;
    OutPr.skn = OutPheno.skn;
    OutPr.srn = OutPheno.srn;
    OutPr.sknb = OutPheno.sknb;
    
    OutPr.espn = OutPheno.espn;
    OutPr.estn = OutPheno.estn;
    OutPr.esbn = OutPheno.esbn;
    OutPr.esfn = OutPheno.esfn;
    OutPr.eskn = OutPheno.eskn;
    OutPr.esrn = OutPheno.esrn;
    OutPr.esknb = OutPheno.esknb;
    
    save([resfolder,filesep,'BPHENO_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutPr','-v7.3');
end

%% Plot Individual Traces
clear
clc
close all

Wi = 470;
Le = 330;

%NC = 12;
%NC = 13;
NC = 14;

GList = {'hb_MS2_INTRON','hb_MS2_3UTR','kr_MS2_INTRON','kr_MS2_3UTR','kni_MS2_INTRON','gt_MS2_INTRON_Female_Ant',...
    'gt_MS2_INTRON_Female_Post','gt_MS2_INTRON_Male_Ant','gt_MS2_INTRON_Male_Post',...
    'kni_MS2_INTRON_hb_null','kni_MS2_INTRON_dist_del','hb_MS2_INTRON_dist_del'};

ll = 1;
line = GList{ll};

% Only ap=16 is provided as example data
ap = 16;

mypath = '.';
resfolder = [mypath,filesep,'OutDeconvRes'];

load([resfolder,filesep,'Data_',line,'_Nc',num2str(NC),'.mat'],'Data');
load([resfolder,filesep,'DECONV_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutTr');
load([resfolder,filesep,'BPHENO_',line,'_Nc',num2str(NC),'_AP',num2str(ap),'.mat'],'OutPr');

Traces = Data(ap).Traces;
T = Data(ap).T;
It = Data(ap).It;
Ti = T(It);
dt = Ti(2)-Ti(1);
cmap = lines(4);

%%% plot time-depedent bursting parameters
H12=figure(12);
set(H12,'position',[50 700 3*Wi 2*Le],'paperpositionmode','auto','color','w');
h12_01=subplot(2,3,2,'parent',H12);
h12_02=subplot(2,3,3,'parent',H12);
h12_03=subplot(2,3,1,'parent',H12);
h12_04=subplot(2,3,5,'parent',H12);
h12_05=subplot(2,3,4,'parent',H12);
h12_06=subplot(2,3,6,'parent',H12);
hold(h12_01,'on')
hold(h12_02,'on')
hold(h12_03,'on')
hold(h12_04,'on')
hold(h12_05,'on')
hold(h12_06,'on')
 
cmap(1,:) = [0,0.5098,0.7843]; %blue

plotEnveloppe(h12_01,Ti,OutPr.pn(2,:),OutPr.epn(2,:),cmap(1,:));
set(h12_01,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_01,'Time [min]')
ylabel(h12_01,'ON-fraction P_{ON}')
xlim(h12_01,[Ti(1),Ti(end)])
Y=ylim(h12_01);
ylim(h12_01,[0,Y(2)])

plotEnveloppe(h12_02,Ti,OutPr.skn(2,:),OutPr.eskn(2,:),cmap(1,:));
set(h12_02,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_02,'Time [min]')
ylabel(h12_02,'Initiation rate K [mRNA/min]')
xlim(h12_02,[Ti(1),Ti(end)])
Y=ylim(h12_02);
ylim(h12_02,[0,Y(2)])

plotEnveloppe(h12_03,Ti,OutPr.srn(2,:),OutPr.esrn(2,:),cmap(1,:));
set(h12_03,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_03,'Time [min]')
ylabel(h12_03,'Transcri. rate R [mRNA/min]')
xlim(h12_03,[Ti(1),Ti(end)])
Y=ylim(h12_03);
ylim(h12_03,[0,Y(2)])

plotEnveloppe(h12_04,Ti,OutPr.stn(2,:),OutPr.estn(2,:),cmap(1,:));
set(h12_04,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_04,'Time [min]')
ylabel(h12_04,'ON-time T_{ON} [min]')
xlim(h12_04,[Ti(1),Ti(end)])
Y=ylim(h12_04);
ylim(h12_04,[0,Y(2)])

plotEnveloppe(h12_05,Ti,OutPr.stn(1,:),OutPr.estn(1,:),cmap(1,:));
set(h12_05,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_05,'Time [min]')
ylabel(h12_05,'OFF-time T_{OFF} [min]')
xlim(h12_05,[Ti(1),Ti(end)])
Y=ylim(h12_05);
ylim(h12_05,[0,Y(2)])

TN = OutPr.stn(1,:).*OutPr.stn(2,:)./(OutPr.stn(1,:)+OutPr.stn(2,:));
ETN = abs(OutPr.stn(2,:)./(OutPr.stn(1,:)+OutPr.stn(2,:) - OutPr.stn(1,:).*OutPr.stn(2,:)./(OutPr.stn(1,:)+OutPr.stn(2,:)).^2)) .*OutPr.estn(1,:) +...
    abs(OutPr.stn(1,:)./(OutPr.stn(1,:)+OutPr.stn(2,:) - OutPr.stn(1,:).*OutPr.stn(2,:)./(OutPr.stn(1,:)+OutPr.stn(2,:)).^2)) .*OutPr.estn(2,:);

plotEnveloppe(h12_06,Ti,TN,ETN,cmap(1,:));
set(h12_06,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h12_06,'Time [min]')
ylabel(h12_06,'Correlation time T_{C} [min]')
xlim(h12_06,[Ti(1),Ti(end)])
Y=ylim(h12_06);
ylim(h12_06,[0,Y(2)])

%%% make raster plots over all nuclei
H17=figure(17);
set(H17,'position',[50 700 Wi 1.2*Le],'paperpositionmode','auto','color','w');
h17=subplot(1,1,1,'parent',H17);
hold(h17,'on')
axis(h17,'square')

H18=figure(18);
set(H18,'position',[50 700 Wi 1.2*Le],'paperpositionmode','auto','color','w');
h18=subplot(1,1,1,'parent',H18);
hold(h18,'on')
axis(h18,'square')

Rrast = nan(length(Traces),length(Ti));
Arast = nan(length(Traces),length(Ti));
%burst calling
wt = 5/6; %averaging window [min]
ki = 2; %at least 2 ini events over wt min
s = wt/2;
dt0 = 0.0150;
vg = -0.5:(0.5+ceil(1/6/dt0));
Pg = zeros(length(vg)-1,1);
T0 = zeros(1000,length(Traces));
for j=1:length(Traces)
    tt = OutTr(j).tt;
    I = false(size(tt));
    for i=1:length(Ti)
        [~,k] = min(abs(tt-Ti(i)));
        I(k) = true;
    end
    ti = tt(I);
        
    Rs = double(OutTr(j).Rs);
    Gs = cumsum(Rs,2);
    Gs = [Gs(:,1),diff(Gs(:,I),1,2)];
    
    pg = histc(Gs(:),vg);
    Pg = Pg + pg(1:(end-1));
    
    dti = [Inf,diff(ti)];
    Gs = Gs./repmat(dti,size(Gs,1),1); %gives R
    
    T = repmat(ti',1,length(ti));
    T = T-T';
    K = exp(-(T/s).^4);
    K = K./repmat(sum(K,1),length(ti),1);
    Ns = nan(size(Gs));
    for i=1:size(Gs,1)
        Ns(i,:) = Gs(i,:)*K;
        kn = find(Ns(i,:),1,'first');
        if isempty(kn)
            kn = length(Ns(i,:));
        end
        T0(i,j) = ti(kn);
    end
    Ns = double(Ns>=ki/wt); %gives On/Off
    
    %average over samples
    Gs = nanmean(Gs,1);
    Ns = nanmean(Ns,1);
    
    if length(ti)<length(Ti)
        Gs = [Gs,nan(1,length(Ti)-length(ti))];
        Ns = [Ns,nan(1,length(Ti)-length(ti))];
    end
    
    Rrast(j,:) = Gs;
    Arast(j,:) = Ns;
end

rmap1 = [linspace(1,0.3,64)',linspace(1,0.3,64)',linspace(1,0.3,64)'];
rmap2 = [linspace(0.65,0,64)',linspace(0.65,0.5098,64)',linspace(0.65,0.7843,64)'];

imagesc(h17,Ti,1:length(Traces),Rrast)
set(h17,'fontsize',22,'linewidth',2,'tickdir','out','YDir','normal','Layer','Top')
xlabel(h17,'Time [min]')
ylabel(h17,'Nuclei#')
xlim(h17,[Ti(1),Ti(end)])
ylim(h17,[1,length(Traces)])
colormap(h17,rmap1)
clim(h17,[0,35])
box(h17,'on')

cb=colorbar(h17,'LineWidth',2);
set(get(cb,'ylabel'),'String','Transcription rate [mRNA/min]','fontsize',22);

imagesc(h18,Ti,1:length(Traces),Arast)
set(h18,'fontsize',22,'linewidth',2,'tickdir','out','YDir','normal','Layer','Top')
xlabel(h18,'Time [min]')
ylabel(h18,'Nuclei#')
xlim(h18,[Ti(1),Ti(end)])
ylim(h18,[1,length(Traces)])
colormap(h18,rmap2)
colormap(h18,[rmap2(1,:);rmap2(end,:)])
clim(h18,[0,1])
box(h18,'on')

cb=colorbar(h18,'LineWidth',2);
set(get(cb,'ylabel'),'String','Locus state','fontsize',22);
cb.Ticks = [0,1];
cb.TickLabels = {'OFF','ON'};

%%% plot deconvolution results for individual nuclei
for j=1:length(Traces)
    t = Traces(j).t;
    int = Traces(j).int;
    if isfield(Traces,'nid')
        rs = Traces(j).rs;
        emb = Traces(j).emb;
        nid = Traces(j).nid;
    else
        rs = nan(4,length(t));
        emb = nan;
        nid = nan;
    end
    
    Rs = OutTr(j).Rs(1,:);
    Ss = OutTr(j).Ss(1,:);
    S = OutTr(j).S;
    tt = OutTr(j).tt;
    II = OutTr(j).I;
    dt0 = tt(2)-tt(1);
    ti = tt(II);
    dti = [Inf,diff(ti)];
    Gs = double(OutTr(j).Gs); %number of ini during dt
    Gs = Gs./repmat(dti,1000,1);

    %burst calling
    wt = 5/6; %averaging window [min]
    ki = 2; %at least 2 ini events over wt min
    s = wt/2;
    T = repmat(ti',1,length(ti));
    T = T-T';
    K = exp(-(T/s).^4);
    K = K./repmat(sum(K,1),length(ti),1);
    Ns = zeros(size(Gs));
    for i=1:size(Gs,1)
        Ns(i,:) = Gs(i,:)*K;
    end
    Ns = double(Ns>=ki/wt);
    Pn = mean(Ns,1);
       
    H1=figure(1);
    set(H1,'position',[50 700 2.3*Wi Le],'paperpositionmode','auto','color','w');
    h1=axes('parent',H1);
    
    yyaxis(h1,'left')
    hold(h1,'on')
    plot(h1,t,int,'-k','linewidth',2);
    tsn = [ti-[0,diff(ti)/2],ti(end)+median(diff(ti))/2];
    yy = [Ss,Ss(end)];
    stairs(h1,tsn,yy,'-','color','r','linewidth',2)
    ylim(h1,[0,60])
    ylabel(h1,'Activity [C.U.]')

    yyaxis(h1,'right')
    bar(h1,tt,Rs,1,'FaceColor',[0.5,0.5,0.5])
    yy = 0.5*[Ns(1,:),Ns(1,end)];
    stairs(h1,tsn,yy,'-','color',cmap(1,:),'linewidth',2);    
    set(h1,'fontsize',22,'linewidth',2,'tickdir','out')
    ylabel(h1,'mRNA initiation')
    ylim(h1,[0,1])
    set(h1,'Ytick',0:1)
    xlim(h1,[Ti(1),Ti(end)])
    xlabel(h1,'Time [min]')
    
    h1.YAxis(1).Color = 'k';
    h1.YAxis(2).Color = [0.5,0.5,0.5];
        
    H2=figure(2);
    set(H2,'position',[50 277 2.3*Wi Le],'paperpositionmode','auto','color','w');
    h2=axes('parent',H2);
    
    yyaxis(h2,'left')
    hold(h2,'on')
    plot(h2,t,int,'-k','linewidth',2);
    plotEnveloppe(h2,ti,S(1,:),S(2,:),[1,0,0]);
    ylim(h2,[0,60])
    ylabel(h2,'Activity [C.U.]')
    
    yyaxis(h2,'right')
    hold(h2,'on')
    plotEnveloppe(h2,ti,mean(Gs,1),std(Gs,1,1),[0.5,0.5,0.5]);  
    plot(h2,ti,(2/5)*Pn/dt0,'-','color',cmap(1,:),'linewidth',2);
    set(h2,'fontsize',22,'linewidth',2,'tickdir','out')
    ylabel(h2,'Transcription rate [mRNA/min]')
    ylim(h2,[0,1/dt0])
    xlim(h2,[Ti(1),Ti(end)])
    xlabel(h2,'Time [min]')
    
    h2.YAxis(1).Color = 'k';
    h2.YAxis(2).Color = [0.5,0.5,0.5];  
    
    H3=figure(3);
    set(H3,'position',[700 0 2.3*Wi Le],'paperpositionmode','auto','color','w');
    h3=axes('parent',H3);
    
    yyaxis(h3,'left')
    hold(h3,'on')
    plot(h3,t,rs(1,:),'-','color',cmap(1,:),'linewidth',2)
    plot(h3,t,rs(2,:),'-','color',cmap(2,:),'linewidth',2)
    plot(h3,t,rs(3,:),'-','color',cmap(3,:),'linewidth',2)
    ylabel('x,y,z')
    
    yyaxis(h3,'right')
    yy = [rs(4,:),rs(4,end)];
    stairs(h3,tsn,yy,'-','color',[0.5,0.5,0.5],'linewidth',2)
    ylabel('Found')
    ylim(h3,[0,1.1])
    set(h3,'Ytick',0:1)
    
    set(h3,'fontsize',22,'linewidth',2,'box','off','tickdir','out')
    xlim(h3,[Ti(1),Ti(end)])
    xlabel(h3,'Time [min]')
    title(h3,['Emb ',num2str(emb),', nuc ',num2str(nid)])
    
    h3.YAxis(1).Color = 'k';
    h3.YAxis(2).Color = [0.5,0.5,0.5];
    
    figure(H1)
    pause
    close(H1)
    close(H2)
    close(H3)
end