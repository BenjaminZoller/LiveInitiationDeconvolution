function OutPheno = ExtractIniPheno(OutTr,Ti,Ns,Ws,Ks)
%ExtractIniPheno performs burst calling based on initiation events and
%computes bursting parameters.
%   OutPheno = ExtractIniPheno(OutTr,Ti,Ns,Ws,Ks) returns the bursting
%   parameters OutPheno. The function takes as input the deconvolved
%   initiation events OutTr, the time Ti, the number of sampled Ns, the
%   averaging window Ws and the threshold Ks for burst calling.
%
%   Copyright (c) 2024, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

COMPCORR = 1;

if nargin < 3 || Ns > 1e3
    Ns = 1e3;
end

db = 1;
dt = Ti(2)-Ti(1);
Nt = length(Ti);
Ntr = length(OutTr);

Lmax = max(arrayfun(@(x) size(x.Gs,2),OutTr));

Tmax = max(arrayfun(@(x) x.tt(end)-x.tt(1),OutTr));
Ndtn = ceil(Tmax/dt);
vdtn = dt*(0:(Ndtn+1))-0.5*dt;
Bmax = 0.5*15*Ndtn;
Ndbn = ceil(Bmax/db);
vdbn = db*(0:(Ndbn+1))-0.5*db;

%for moments
TT = nan(Ntr,Lmax);
for k=1:Ntr
    tt = OutTr(k).tt(OutTr(k).I);
    l = length(tt);
    TT(k,1:l) = tt;
end

I = cell(Nt,1);
W = cell(Nt,1);
for j=1:Nt
    Vj = abs((TT-Ti(j))/dt);
    Ij = Vj <= 1;
    
    Wj = 1-Vj(Ij);
    I{j} = find(Ij); %will be faster for GG(Ij);
    W{j} = Wj;
end

%for correlations
if COMPCORR
    H = logical(toeplitz([1,zeros(1,Nt-1)],ones(1,Nt)));
    hi = repmat((1:Nt)',1,Nt);
    hj = repmat((1:Nt),Nt,1);
    Ihi = hi(H);
    Ihj = hj(H);
    Ntt = length(Ihi);
    GCC = nan(Ns,Ntt);
end

%precompute ini clustering
wt = Ws; %averaging window [min]
kt = Ks; %at least kt ini events per min over wt
s = wt/2;
for k=1:Ntr
    Gs = double(OutTr(k).Gs);
    
    l = size(Gs,2);
    t = TT(k,1:l);
    dtk = [Inf,diff(t)];
    Gs = Gs./repmat(dtk,Ns,1);
    
    T = repmat(t',1,l);
    T = T-T';
    K = exp(-(T/s).^4);
    K = K./repmat(sum(K,1),l,1);
    
    NNs = Gs*K;
    NNs = double(NNs>=kt/wt);
    IniClusters(k).NNs = NNs;
end

%moments
G1 = nan(Ns,Nt);
G2 = nan(Ns,Nt);
A1 = nan(Ns,Nt);
A2 = nan(Ns,Nt);
N1 = nan(Ns,Nt);

%off & on
TTN = nan(2,Ns,Nt);
TCN = nan(2,Ns,Nt);

KN = nan(2,Ns,Nt);
ZN = nan(2,Ns,Nt);
TN = nan(2,Ns,Nt); 
BN = nan(2,Ns,Nt);

%distributions
WT = zeros(2,Ns,Ndtn+1);
WB = zeros(2,Ns,Ndbn+1);
WF = zeros(Ns,Ndtn+1);

% main loop
for s=1:Ns
    disp(s)
    %profile on
    GG = nan(Ntr,Lmax);
    AG = nan(Ntr,Lmax);
    NN = nan(Ntr,Lmax);
    
    if COMPCORR
        CC = nan(Ntr,Ntt);
        Xi = nan(Ntr,Ntt);
        Xj = nan(Ntr,Ntt);
    end
    
    NTT = nan(2,Ntr,Lmax);
    NTC = nan(2,Ntr,Lmax);
    NK = nan(2,Ntr,Lmax);
    NT = nan(2,Ntr,Lmax);
    NB = nan(2,Ntr,Lmax);
    NZ = nan(2,Ntr,Lmax);
    
    Wt = zeros(2,Ndtn+1);
    Wb = zeros(2,Ndbn+1);
    Wf = zeros(1,Ndtn+1);
    
    for k=1:Ntr
        %get all ini samples
        gg = double(OutTr(k).Gs(s,:));
        aa = cumsum(gg);
        l = length(gg);
        t = TT(k,1:l);
        dtk = [Inf,diff(t)];
        gg = gg./dtk;
                
        %get ini clusters
        nn = IniClusters(k).NNs(s,:);
             
        GG(k,1:l) = gg;
        AG(k,1:l) = aa;
        NN(k,1:l) = nn;
        
        %for correlations
        if COMPCORR
            gi = interp1(t,gg,Ti,'nearest');
            
            xi = gi(Ihi);
            xj = gi(Ihj);
            cc = xi.*xj;
            In = isnan(cc);
            xi(In) = nan;
            xj(In) = nan;
            
            Xi(k,:) = xi;
            Xj(k,:) = xj;
            CC(k,:) = cc;
        end
        
        %for off and on clusters
        IN = logical([~nn;nn]);
        DI = diff([zeros(2,1),double(IN),zeros(2,1)],1,2);
        for n=1:2
            %ini rate
            kk = nan(size(gg));
            kk(IN(n,:)) = gg(IN(n,:));
            NK(n,k,1:l) = kk;
            %total duration
            tt = [0,diff(t)];
            tt(~IN(n,:)) = 0;
            NTT(n,k,1:l) = cumsum(tt);
            %burst durations and sizes
            %kn and ke have the same length
            kb = find(DI(n,:) > 0);
            ke = find(DI(n,:) < 0);
            %total burst count
            cc = zeros(1,l);
            cc(kb) = 1;
            NTC(n,k,1:l) = cumsum(cc);

            if isempty(kb)
                continue
            end
            
            kt = kb;
            kb = kb-1;
            ke = ke-1;
            if kb(1)==0
                kb(1)=1;
                if ke(1)==1
                    kt(1) = [];
                    kb(1) = [];
                    ke(1) = [];
                    if isempty(kb)
                        continue
                    end
                end
            end
            
            tn = t(ke)-t(kb);
            bn = aa(ke)-aa(kb);
            
            %intervals
            for i=1:length(kb)
                k1 = kt(i);
                k2 = ke(i);
                NT(n,k,k1:k2) = tn(i);
                NB(n,k,k1:k2) = bn(i);
                NZ(n,k,k1:k2) = 1/length(k1:k2);
            end
 
            %for distributions
            tn = tn(:)';
            ctn = histc(tn,vdtn);
            Wt(n,:) = Wt(n,:) + ctn(1:(end-1));
            bn = bn(:)';
            cbn = histc(bn,vdbn);
            Wb(n,:) = Wb(n,:) + cbn(1:(end-1));
            
            %if off for first passage time
            if n==1
                ct0 = histc(tn(1),vdtn);
                Wf = Wf + ct0(1:(end-1));
            end
        end
    end
    
    %averaging
    for n=1:2
        NTTn = squeeze(NTT(n,:,:));
        NTCn = squeeze(NTC(n,:,:));
        
        NKn = squeeze(NK(n,:,:));
        NTn = squeeze(NT(n,:,:));
        NBn = squeeze(NB(n,:,:));
        NZn = squeeze(NZ(n,:,:));
        
        for j=1:Nt
            Ij = I{j};
            Wj = W{j};
            Zj = sum(Wj);
            
            if n==1
                Gj = GG(Ij);
                Aj = AG(Ij);
                Nj = NN(Ij);
                
                G1(s,j) = sum(Wj.*Gj)/Zj;
                G2(s,j) = sum(Wj.*Gj.^2)/Zj;
                A1(s,j) = sum(Wj.*Aj)/Zj;
                A2(s,j) = sum(Wj.*Aj.^2)/Zj;
                N1(s,j) = sum(Wj.*Nj)/Zj;
            end
            
            TTN(n,s,j) = sum(Wj.*NTTn(Ij))/Zj;
            TCN(n,s,j) = sum(Wj.*NTCn(Ij))/Zj;
            
            NKj = NKn(Ij);
            NZj = double(~isnan(NKj));
            Zj = nansum(Wj.*NZj);
            KN(n,s,j) = nansum(Wj.*NKj)/Zj;
            
            %full
            NTj = NTn(Ij);
            NBj = NBn(Ij);
            NZj = NZn(Ij);
            
            Zj = nansum(Wj.*NZj);
            TN(n,s,j) = nansum(Wj.*NZj.*NTj)/Zj;
            BN(n,s,j) = nansum(Wj.*NZj.*NBj)/Zj;
            ZN(n,s,j) = Zj;
        end
    end
    
    %correlations
    if COMPCORR
        mi = nanmean(Xi,1);
        mj = nanmean(Xj,1);
        si = nanstd(Xi,1,1);
        sj = nanstd(Xj,1,1);
        si(si<1e-14) = 0;
        sj(sj<1e-14) = 0;
        ccij = (nanmean(CC,1)-mi.*mj)./(si.*sj);
        
        ccij(isinf(ccij)) = nan;
        GCC(s,:) = ccij;
    end
    
    %distributions
    WF(s,:) = Wf/sum(Wf);
    for n=1:2
        WT(n,s,:) = Wt(n,:)/sum(Wt(n,:));
        WB(n,s,:) = Wb(n,:)/sum(Wb(n,:));
    end
    
    %profile viewer
end

G2 = G2-G1.^2;
NG2 = G2./G1.^2;
NG1 = sqrt(NG2);

OutPheno.G1 = [nanmean(G1,1);nanstd(G1,0,1)];
OutPheno.G2 = [nanmean(G2,1);nanstd(G2,0,1)];
OutPheno.NG1 = [nanmean(NG1,1);nanstd(NG1,0,1)];
OutPheno.NG2 = [nanmean(NG2,1);nanstd(NG2,0,1)];

A2 = A2-A1.^2;
NA2 = A2./A1.^2;
NA1 = sqrt(NA2);

OutPheno.A1 = [nanmean(A1,1);nanstd(A1,0,1)];
OutPheno.A2 = [nanmean(A2,1);nanstd(A2,0,1)];
OutPheno.NA1 = [nanmean(NA1,1);nanstd(NA1,0,1)];
OutPheno.NA2 = [nanmean(NA2,1);nanstd(NA2,0,1)];

%kernel smoothed pheno
T = repmat(1/6*(0:(Nt-1))',1,Nt);
T = T-T';

%correlation
if COMPCORR
    MCorr = nanmean(GCC,1);
    SCorr = nanstd(GCC,1,1);
    temp1 = nan(Nt,Nt);
    temp2 = nan(Nt,Nt);
    temp1(H) = MCorr;
    temp2(H) = SCorr;
    Gcorr = nan(2,Nt,Nt);
    for j=1:Nt
        for l=0:(Nt-j)
            Gcorr(1,j,l+1) = temp1(j,j+l);
            Gcorr(2,j,l+1) = temp2(j,j+l);
        end
    end
    OutPheno.Gcorr = Gcorr;
end

%burst pheno
PN = nan(2,Ns,Nt);
PN(1,:,:) = 1-N1; %off
PN(2,:,:) = N1; %on
RN = KN.*PN;
FN = PN./TN; %FNb
KNb = BN./TN;

OutPheno.zn = squeeze(nanmean(ZN,2));
OutPheno.pn = squeeze(nanmean(PN,2));
OutPheno.tn = squeeze(nanmean(TN,2));
OutPheno.bn = squeeze(nanmean(BN,2));
OutPheno.fn = squeeze(nanmean(FN,2));
OutPheno.kn = squeeze(nanmean(KN,2));
OutPheno.rn = squeeze(nanmean(RN,2));
OutPheno.knb = squeeze(nanmean(KNb,2));
OutPheno.ttn = squeeze(nanmean(TTN,2));
OutPheno.tcn = squeeze(nanmean(TCN,2));

OutPheno.epn = squeeze(nanstd(PN,1,2));
OutPheno.etn = squeeze(nanstd(TN,1,2));
OutPheno.ebn = squeeze(nanstd(BN,1,2));
OutPheno.efn = squeeze(nanstd(FN,1,2));
OutPheno.ekn = squeeze(nanstd(KN,1,2));
OutPheno.ern = squeeze(nanstd(RN,1,2));
OutPheno.eknb = squeeze(nanstd(KNb,1,2));
OutPheno.ettn = squeeze(nanstd(TTN,1,2));
OutPheno.etcn = squeeze(nanstd(TCN,1,2));

%distributions
OutPheno.wf = nanmean(WF,1);
OutPheno.wt = squeeze(nanmean(WT,2));
OutPheno.wb = squeeze(nanmean(WB,2));

OutPheno.cf = nanmean(cumsum(WF,2),1);
OutPheno.ct = squeeze(nanmean(cumsum(WT,3),2));
OutPheno.cb = squeeze(nanmean(cumsum(WB,3),2));

OutPheno.ecf = nanstd(cumsum(WF,2),1,1);
OutPheno.ect = squeeze(nanstd(cumsum(WT,3),1,2));
OutPheno.ecb = squeeze(nanstd(cumsum(WB,3),1,2));

OutPheno.vt = dt*(0:Ndtn);
OutPheno.vb = db*(0:Ndbn);

%kernel smoothed pheno
ws = 1/2.355; %1min half maximum
K = normpdf(T,0,ws);

sPN = nan(size(PN));
sTN = nan(size(TN));
sBN = nan(size(BN));
sKN = nan(size(KN));
for n=1:2
    sPN(n,:,:) = myfKernelDensity(squeeze(PN(n,:,:)),K);
    sTN(n,:,:) = myfKernelDensity(squeeze(TN(n,:,:)),K);
    sBN(n,:,:) = myfKernelDensity(squeeze(BN(n,:,:)),K);
    sKN(n,:,:) = myfKernelDensity(squeeze(KN(n,:,:)),K);
end
sFN = sPN./sTN;
sRN = sKN.*PN;
sKNb = sBN./sTN;

OutPheno.spn = squeeze(nanmean(sPN,2));
OutPheno.stn = squeeze(nanmean(sTN,2));
OutPheno.sbn = squeeze(nanmean(sBN,2));
OutPheno.sfn = squeeze(nanmean(sFN,2));
OutPheno.skn = squeeze(nanmean(sKN,2));
OutPheno.srn = squeeze(nanmean(sRN,2));
OutPheno.sknb = squeeze(nanmean(sKNb,2));

OutPheno.espn = squeeze(nanstd(sPN,1,2));
OutPheno.estn = squeeze(nanstd(sTN,1,2));
OutPheno.esbn = squeeze(nanstd(sBN,1,2));
OutPheno.esfn = squeeze(nanstd(sFN,1,2));
OutPheno.eskn = squeeze(nanstd(sKN,1,2));
OutPheno.esrn = squeeze(nanstd(sRN,1,2));
OutPheno.esknb = squeeze(nanstd(sKNb,1,2));

end

function Y = myfKernelDensity(X,K)

Y = zeros(size(X));
for i=1:size(X,1)
    xi = X(i,:);
    In = isnan(xi);
    Ki = K;
    xi(In) = 0;
    Ki(repmat(In',1,size(X,2))) = 0;
    Y(i,:) = xi*Ki./sum(Ki,1);
end

end
