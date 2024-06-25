function [OutTr,OutPr] = DeconvSampling(Traces,M1,T,It,Kernel,Noise,Ke)
%DeconvSampling performs the Bayesian deconvolution of initiation events 
%from transcriptional time traces.
%   [OutTr,OutPr] = DeconvSampling(Traces,M1,T,It,Kernel,Noise,Ke) returns 
%   the sampled deconvolution OutTr and other probabilities OutPr using a
%   MCMC approach. The function takes as input the Traces, the mean 
%   activity M1, the time T and its valid subset It, the deconvolution
%   Kernel, the measurement Noise and the elongation rate Ke.
%
%   Copyright (c) 2024, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

FP = 60; %Pol2 footprint
dt0 = 0.5*FP/Ke; %assuming 2 sister chromatids
dt = T(2)-T(1);
gmax = 1+ceil(dt/dt0);

Lg = sum(diff(Kernel.X).*Kernel.K(1:(end-1))); %effective gene length
Te = Lg/Ke; %effective elongation time

ws = 4*round(1/6/dt0);
Ns = 3.5e3; %number of samples
burnin = 5e2;

noise = @(x) sqrt(Noise(x));

Nd = length(Traces);

Ps = cell(Nd,1);
Pg = cell(Nd,1);
Pg3 = cell(Nd,1);
Pg6 = cell(Nd,1);
Pg12 = cell(Nd,1);
Tp = cell(Nd,1);
Wr = cell(Nd,1);
Wf = cell(Nd,1);

Ti = T(It);
T0 = Ti(1);

Ndtr = ceil((Ti(end)-Ti(1))/dt0);
vdtr = dt0*(0:(Ndtr+1))-0.5*dt0;

%for xx=1:Nd
parfor xx=1:Nd
    disp(xx)
    
    t = Traces(xx).t;
    int = Traces(xx).int;
    
    % setup support
    Nt = round((t(end)-t(1))/dt0);
    tt = t(1)+dt0*(0:Nt);
    I = false(size(tt));
    
    for i=1:length(t)
        [~,k] = min(abs(tt-t(i)));
        I(k) = true;
    end
    
    Nm = round((Nt+1)/ws); %for mixing
    
    OutTr(xx).tt = single(tt);
    OutTr(xx).I = I;
    OutTr(xx).int = single(int);
    
    %%% setup variables
    S = zeros(Ns,Nt+1);
    R = zeros(Ns,Nt+1);
    
    %%% ini sampling variables
    % kernel (Ke dependent)
    kk = makeKernel(Kernel,[],Ke,dt0);
    % proposal(Ke dependent)
    p = dt0*smooth(int)/Te;
    p = abs(interp1(tt(I),p,tt+0.5*Te,'makima',mean(p((end-10):end))));
    p(p < 0.05) = 0.05;
    % prior (Ke dependent)
    rr = dt0*M1/Te;
    pr = abs(interp1(T,rr,tt+0.5*Te,'makima',rr(end)));
    % signal
    r = binornd(1,p,1,Nt+1);
    s = conv(r,kk);
    s = s(1:(Nt+1));
    
    L = logLikelihood(int,s(I),noise);

    for ii = 1:Ns
        for jj=1:Nm
            % setting window
            ir = randi(Nt+1);
            wr = ir+(0:(ws-1));
            wr = 1+mod(wr-1,Nt+1);
            
            % propose update
            rprop = r;
            rprop(wr) = binornd(1,p(wr),1,ws);
            sprop = conv(rprop,kk);
            sprop = sprop(1:(Nt+1));
            
            Lprop = logLikelihood(int,sprop(I),noise);
            Q = sum(log( p(wr).^r(wr) .* (1-p(wr)).^(1-r(wr)) ));
            Qprop = sum(log( p(wr).^rprop(wr) .* (1-p(wr)).^(1-rprop(wr)) ));
            P = sum(log( pr(wr).^r(wr) .* (1-pr(wr)).^(1-r(wr)) ));
            Pprop = sum(log( pr(wr).^rprop(wr) .* (1-pr(wr)).^(1-rprop(wr)) ));
            
            ar = min(1,exp(Lprop+Pprop-L-P+Q-Qprop));
            
            % accepting & updating
            if rand(1) <= ar
                L = Lprop;
                r = rprop;
                s = sprop;
            end
            
            % learning proposal
            if ii > 1
                it = jj+Nm*(ii-1);
                lr = 1/(3*sqrt(it));
            else
                lr = 0;
            end
            
            pnew = r(wr);
            pnew(pnew<1) = 0.05;
            p(wr) = p(wr)*(1-lr) + lr*pnew;
        end
        
        R(ii,:) = r;
        S(ii,:) = s;
    end
    
    S = S((1+burnin):end,:);
    R = R((1+burnin):end,:);
    
    G = cumsum(R,2);
    G = [G(:,1),diff(G(:,I),1,2)];

    OutTr(xx).S = single([mean(S(:,I),1);std(S(:,I),0,1)]);
    OutTr(xx).R = single(mean(R,1)); %is an occupancy, std = sqrt(r(1-r))
    OutTr(xx).G = single([mean(G,1);std(G,0,1)]);
    
    % example samples
    Ir = randperm(Ns-burnin,1e3);
    OutTr(xx).Ss = single(S(Ir(1:10),I));
    OutTr(xx).Rs = logical(R(Ir,:));
    OutTr(xx).Gs = uint8(G(Ir,:));
    
    % distributions
    [Ps{xx},Pg{xx},Pg3{xx},Pg6{xx},Pg12{xx},Tp{xx},Wr{xx},Wf{xx}] = getFeatures(R,S,tt,Ti,gmax,vdtr,T0);
end

% Sum up states distribution
Nt = length(Ti);

smax = max(cellfun(@(x) size(x,1),Ps));
gmax3 = round(3*gmax/2 + sqrt(3)*gmax/2);
gmax6 = round(6*gmax/2 + sqrt(6)*gmax/2);
gmax12 = round(12*gmax/2 + sqrt(12)*gmax/2);

%Mixture distribution
PS = zeros(smax,Nt);
PG = zeros(gmax,Nt);
PG3 = zeros(gmax3,Nt);
PG6 = zeros(gmax6,Nt);
PG12 = zeros(gmax12,Nt);
ZZ = zeros(1,Nt);

for xx=1:Nd
    t = Tp{xx};
    ps = Ps{xx};
    pg = Pg{xx};
    pg3 = Pg3{xx};
    pg6 = Pg6{xx};
    pg12 = Pg12{xx};
    si = size(ps,1);
    for j=1:Nt
        vj = abs((t-Ti(j))/dt);
        Ij = vj <= 1;
        wj = 1-vj(Ij);
        PS(1:si,j) = PS(1:si,j) + sum(ps(:,Ij).*repmat(wj,si,1),2);
        PG(:,j) = PG(:,j) + sum(pg(:,Ij).*repmat(wj,gmax,1),2);
        PG3(:,j) = PG3(:,j) + sum(pg3(:,Ij).*repmat(wj,gmax3,1),2);
        PG6(:,j) = PG6(:,j) + sum(pg6(:,Ij).*repmat(wj,gmax6,1),2);
        PG12(:,j) = PG12(:,j) + sum(pg12(:,Ij).*repmat(wj,gmax12,1),2);
        ZZ(j) = ZZ(j) + sum(wj);
    end
end

% Normalization
PS = PS./repmat(ZZ,smax,1);
PG = PG./repmat(ZZ,gmax,1);
PG3 = PG3./repmat(ZZ,gmax3,1);
PG6 = PG6./repmat(ZZ,gmax6,1);
PG12 = PG12./repmat(ZZ,gmax12,1);

% Renormalize for safety
PS = PS./repmat(nansum(PS,1),smax,1);
PG = PG./repmat(nansum(PG,1),gmax,1);
PG3 = PG3./repmat(nansum(PG3,1),gmax3,1);
PG6 = PG6./repmat(nansum(PG6,1),gmax6,1);
PG12 = PG12./repmat(nansum(PG12,1),gmax12,1);

% Sum up time distributions
WR = zeros(1,Ndtr+1);
ZR = zeros(1,Ndtr+1);
WF = zeros(1,Ndtr+1);
ZF = zeros(1,Ndtr+1);

for xx=1:Nd
    l = size(Wr{xx},2);
    WR(1:l) = WR(1:l) + Wr{xx}(1:l);
    ZR(1:l) = ZR(1:l) + 1;
    l = size(Wf{xx},2);
    WF(1:l) = WF(1:l) + Wf{xx}(1:l);
    ZF(1:l) = ZF(1:l) + 1;
end

WR = WR./ZR;
I = isnan(WR);
WR(I) = 0;
WR = WR/sum(WR(2:end));

WF = WF./ZF;
I = isnan(WF);
WF(I) = 0;
WF = WF/sum(WF);

WT = dt0*(0:Ndtr);

% generate output
OutPr.PS = PS;
OutPr.PG = PG;
OutPr.PG3 = PG3;
OutPr.PG6 = PG6;
OutPr.PG12 = PG12;
OutPr.ZZ = ZZ;
OutPr.WR = WR;
OutPr.WF = WF;
OutPr.WT = WT;

end

function L = logLikelihood(int,s,noise)
%loglikelihood of data based on signal and noise
sig = noise(s);
L = -log(sqrt(2*pi)*sig) -0.5*((int-s)./sig).^2;
L = sum(L);
end

function [Ps,Pg,Pg3,Pg6,Pg12,Tp,Wr,Wf] = getFeatures(R,S,tt,Ti,gmax,vdtr,T0)
% setup time grid based on Ti
I = false(size(tt));
i1=find(Ti>tt(1),1,'first');
i2=find(Ti<tt(end),1,'last');
for i=i1:i2
    [~,k] = min(abs(tt-Ti(i)));
    I(k) = true;
end
Tp = tt(I);

Nt = length(Tp);
Ns = size(R,1);

smax = 1+ceil(max(S(:)));
gmax3 = round(3*gmax/2 + sqrt(3)*gmax/2);
gmax6 = round(6*gmax/2 + sqrt(6)*gmax/2);
gmax12 = round(12*gmax/2 + sqrt(12)*gmax/2);

vs = (0:smax)-0.5;
vg = (0:gmax)-0.5;
vg3 = (0:gmax3)-0.5;
vg6 = (0:gmax6)-0.5;
vg12 = (0:gmax12)-0.5;

S = S(:,I);
G = cumsum(double(R),2);
G = [G(:,1),diff(G(:,I),1,2)];

Ps = zeros(smax,Nt);
Pg = zeros(gmax,Nt);
Pg3 = zeros(gmax3,Nt);
Pg6 = zeros(gmax6,Nt);
Pg12 = zeros(gmax12,Nt);

%for states distribution
for i=1:Nt
    sc = histc(S(:,i)',vs);
    Ps(:,i) = sc(1:(end-1))';
    Ps(:,i) = Ps(:,i)/sum(Ps(:,i));
    
    gc = histc(G(:,i)',vg);
    Pg(:,i) = gc(1:(end-1))';
    Pg(:,i) = Pg(:,i)/sum(Pg(:,i));
    
    vt = i:-1:(i-2);
    vt(vt<1) = [];
    gc3 = histc(sum(G(:,vt),2)',vg3);
    Pg3(:,i) = gc3(1:(end-1))';
    Pg3(:,i) = Pg3(:,i)/sum(Pg3(:,i));
    
    vt = i:-1:(i-5);
    vt(vt<1) = [];
    gc6 = histc(sum(G(:,vt),2)',vg6);
    Pg6(:,i) = gc6(1:(end-1))';
    Pg6(:,i) = Pg6(:,i)/sum(Pg6(:,i));
    
    vt = i:-1:(i-11);
    vt(vt<1) = [];
    gc12 = histc(sum(G(:,vt),2)',vg12);
    Pg12(:,i) = gc12(1:(end-1))';
    Pg12(:,i) = Pg12(:,i)/sum(Pg12(:,i));
end

%for r waiting time distribution
Tmax = tt(end)-tt(1);
k = find(vdtr>Tmax,1);
vdtr = vdtr(1:k);

Wr = zeros(1,length(vdtr)-1); %waiting time
Wf = zeros(1,length(vdtr)-1); %first passge
for i=1:Ns
    k = find(R(i,:));
    dtr = tt(k(2:end))-tt(k(1:(end-1)));
    
    wr = histc(dtr,vdtr);
    Wr = Wr + wr(1:(end-1));
    
    if isempty(k)
        dtr = tt(end)-T0;
    else
        dtr = tt(k(1))-T0;
    end
    wf = histc(dtr,vdtr);
    Wf = Wf + wf(1:(end-1));
end
Wr = Wr/Ns;
Wf = Wf/Ns;
end

