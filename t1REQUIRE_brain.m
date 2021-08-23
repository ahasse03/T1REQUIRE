function [t1map_require, brainFinal] = t1REQUIRE_brain(T1w, WM, GM, CSF, infoT1)

%--------------------------------------------------------------------------
%Function to retrospectively estimate T1 using internal references in the
%brain from t1 weighted spin echo image.
%Includes proton density and t2 corrections, along with skull stripping
%Inputs: T1w - t1 weighted image (spin ech or MPRAGE)
%        WM - white matter segmentation
%        GM - gray matter segmentation
%        CSF - cerebrospinal fluid segmentation
%        infoT1 - structure including header info from DICOM/PARREC image
%Outputs: t1map_require - estimated T1 map
%
%Healthy control study has shown this to be effective from ~500 to 3000 ms
%(within 10%)
%
%Developed by Adam Hasse
%University of Chicago Graduate Program in Medical Physics
%University of Chicago Department of Radiology
%Update 28 May 2021
%-------------------------------------------------------------------------

warning off
sz = size(T1w);

%Initialize T1s for GM, WM, and CSF
if round(infoT1.MagneticFieldStrength,1) == 1.5
    wmT1 = 623;
    gmT1 = 1039;
    csfT1 = 2000;
elseif infoT1.MagneticFieldStrength == 3 %NMR relaxation times in the human brain at 3.0 tesla.
    wmT1 = 832;
    gmT1 = 1331;
    csfT1 = 2500;
else
    error('Magnetic field strength is not 1.5 or 3T');
end

%Initialize T2s for GM and WM
gmT2 = 110;
wmT2 = 80;


%Perform skull stripping
brainOnly = (WM>0.9) + (GM>0.9) + (CSF>0.9);
brainDil = [];
brainFill = [];
brainFinal = [];
for j = 1:size(brainOnly,3)
    brainDil(:,:,j) = imdilate(brainOnly(:,:,j), strel('disk',2,0));
    brainFill(:,:,j) = imfill(brainDil(:,:,j), 'holes');
    brainFinal(:,:,j) = imerode(brainFill(:,:,j), strel('disk',2,0));
end

%% Spin Echo Sequence
if strcmp(infoT1.SeriesDescription, 'SE')
    
    %Get variables from infoT1 structure
    TR = infoT1.RepetitionTime;
    TE = infoT1.EchoTime;
    
    
    %Generate T2/PD correction -- assume PD of CSF is 1 and e^(-TE/T2) ~ 1 for
    %CSF
    t2pdCorr = (((GM>0.5).*exp(TE/110)./0.807 + (WM>0.5).*exp(TE/80)./0.679)-1+brainFinal);
    t2pdCorr(t2pdCorr == 0) = 1;
    t2pdCorr = t2pdCorr.*brainFinal;
    t2pdCorr(t2pdCorr == 0) = NaN;
    h = fspecial('gaussian', [13 13], 3);
    for j = 1:size(t2pdCorr,3)
        t2pdCorr(:,:,j) = nanconv(t2pdCorr(:,:,j), h);
    end
    t2pdCorr = t2pdCorr.*brainFinal;
    t1w_t2pdcorr = T1w.*t2pdCorr;
    
    %Run REQUIRE Algorithm on corrected image
    wmValsT1 = (WM>0.9).*t1w_t2pdcorr;
    gmValsT1 = (GM>0.9).*t1w_t2pdcorr;
    csfValsT1 = (CSF>0.9).*t1w_t2pdcorr;
    
    fun = @(x,xdataT1)x(1).*(1-exp(-TR./xdataT1));
    x0 = [-1];
    opts = optimset('Display','off');
    
    converteqT1pre = '-TR./log(1-imgT1slicePre./a)';
    for i = 1:sz(3) %slice-by-slice fit
        if max(max(wmValsT1(:,:,i))) > 0.9 && max(max(gmValsT1(:,:,i))) > 0.9
            avgWMT1 = mean(nonzeros(wmValsT1(:,:,i)), 'omitnan');
            avgGMT1 = mean(nonzeros(gmValsT1(:,:,i)), 'omitnan');
            avgCSFT1 = mean(nonzeros(csfValsT1(:,:,i)), 'omitnan');
            ydataT1 = [avgWMT1 avgGMT1 avgCSFT1]; %Signal intensities
            xdataT1 = [wmT1 gmT1 csfT1];
            xpre = lsqcurvefit(fun, x0, xdataT1, ydataT1, [],[], opts);
            a = xpre(1);
            imgT1slicePre = double(t1w_t2pdcorr(:,:,i));
            quantT1preNoLim(:,:,i) = abs(eval(converteqT1pre));
            stepQuantT1pre = double(quantT1preNoLim(:,:,i) > 0).*quantT1preNoLim(:,:,i);
            stepQuantT1pre(quantT1preNoLim(:,:,i) > 3000) = 0;
            stepQuantT1pre(isnan(stepQuantT1pre)) = 0;
            quantT1pre(:,:,i) = stepQuantT1pre;
        else
            quantT1pre(:,:,i) = zeros(sz(1),sz(2));
            quantT1preNoLim(:,:,i) = zeros(sz(1),sz(2));
        end
    end
    
    
    t1map_require = quantT1pre.*brainFinal;
    
    %% MPRAGE
elseif sum(strcmp(infoT1.SeriesDescription, {'3DT1TFE', 'MPRAGE'})) > 0
    
    %mprage_solver based off of MPRAGE: A Three-dimensional, T1-weighted, Gradient-Echo Sequence -- Initial Experience in the Brain
    %by Michael Brant-Zawadzki, Gary D. Gillan, Wolfgang R. Nitz (Radiology 1992; 182:769-775)
    %
    %July 14, 2020
    %Adam Hasse
    %Graduate Program in Medical Physics, University of Chicago
    %ahasse@uchicago.edu
    %
    %The purpose of this code is to generate a lookup table that will take
    %signal intensities and put out estimated T1 values to convert a
    %T1-weighted MPRAGE image into a t1map
    %Important notes: There are quite a few assumptions that will be made. Most
    %importantly is the lack of any T2* weighting. This is considered
    %appropriate due to the small flip angle in MPRAGE (i.e. 9 degrees in our
    %case) and the small TE (on the order of 2-3 ms, T2* is on the order of
    %tens of milliseconds).
    
    %INPUTS
    %mprage -- MPRAGE t1-weighted image
    %gm -- gray matter segmentation direct from spm12
    %csf -- caf segmentation direct from spm12
    %wm -- wm segmentation direct from spm12
    %n -- phase encoding steps (echo train length)
    %trec -- time between GRE readout and following inversion pulse
    %TI -- inversion time between inversion pulse and GRE readout
    %TR -- repetition time
    %a -- flip angle for GRE readout
    %TE -- echo time for T2 correction
    
    mprage = T1w;
    try
        n = infoT1.GradientEchoTrainLength;
    catch
        n = infoT1.EchoTrainLength;
    end
    TR = 8;
    try
        a = infoT1.FlipAngle;
    catch
        a = 15;
    end
    TI = 358; %From scanner, min inversion delay
    trec = 400;
    n = double(n);
    mprage = mprage.*brainFinal;
    
    
    %Initialize M0 ranges
    M0_step = linspace(0,(10*max(max(max(mprage)))),100);
    M0_step_prev = M0_step;
    t1vals = [1331, 3000, 832]; %ms, GM, CSF, wm
    
    %Get average tissue signal intensities
    gm = GM > 0.9;
    wm = WM > 0.9;
    csf = CSF > 0.9;
    wmValsT1 = wm.*mprage;
    gmValsT1 = gm.*mprage;
    csfValsT1 = csf.*mprage;
    wmMean = mean(nonzeros(wmValsT1));
    gmMean = mean(nonzeros(gmValsT1));
    csfMean = mean(nonzeros(csfValsT1));
    
    %----------------------------------------------------
    
    %Iterate to find the best value of M0
    
    while M0_step_prev(2) - M0_step_prev(1) >= (0.001*max(max(max(mprage))))
        count = 1;
        for j = M0_step
            [SN] = mprage_equation(j,n,trec,TI,TR,a,t1vals);
            rmse(count) = sum(([gmMean, csfMean, wmMean]-SN).^2/length(M0_step));
            %-------------------------------------------------
            count = count + 1;
        end
        [~, ind] = min(rmse);
        M0_step_prev = M0_step;
        M0_step = linspace(M0_step_prev(ind-1), M0_step_prev(ind+1), 100);
    end
    M0 = M0_step(ind);
    rmse_final = rmse(ind);
    
    %Now that we have estimted M0, we will calculate the signal intensity for
    %T1s from 100ms to 3500 ms.
    lookup_table = 100:4000;
    lookup_table = lookup_table';
    [SN] = mprage_equation(M0, n, trec, TI, TR, a, lookup_table);
    lookup_table(:,2) = SN;
    
    %If you want to see how well it ended up fitting, uncomment the following 3
    %lines
%     figure(); plot(t1vals, [gmMean, csfMean, wmMean], 'o')
%     hold on;
%     plot(lookup_table(:,1), lookup_table(:,2))
    %-----------------------------
    
    %Finally, now that we have a lookup table, we can convert each signal
    %intensity into a T1 estimate
    sz = size(mprage);
    t1map_require = zeros(sz);
    for j = 1:sz(1)
        for k = 1:sz(2)
            for m = 1:sz(3)
                if mprage(j,k,m) >= min(lookup_table(:,2)) && mprage(j,k,m) <= max(lookup_table(:,2))
                    try
                        t1map_require(j,k,m) = lookup_table((lookup_table(:,2) == mprage(j,k,m)),1);
                    catch
                        diff = lookup_table(:,2) - mprage(j,k,m);
                        [~, ind] = min(abs(diff));
                        t1map_require(j,k,m) = lookup_table(ind,1);
                    end
                end
            end
        end
    end
    
    
end
end %t1REQUIRE_brain


function [SN] = mprage_equation(M0, n, trec, TI, TR, a, t1vals)
%MPRAGE equation -- easier to do this way

%Begin with X
sumX = 0;
for q = double(0:(n/2-2))
    sumX = sumX + (cosd(a).*exp(-TR./t1vals)).^q;
end
X = (1-exp(-TR./t1vals)).*sumX+(1-2.*exp(-TI./t1vals)).*(cosd(a).*exp(-TR./t1vals)).^(n/2-1);

%Move on to Meq
Meq = 1./(1-X.*exp(-trec./t1vals)).*M0.*(1-exp(-trec./t1vals));

%Now onto MN
sumMN = 0;
for q = double(0:(n/2-2))
    sumMN = sumMN + (cosd(a).*exp(-TR./t1vals)).^q;
end
MN = M0.*(1-exp(-TR./t1vals)).*sumMN + M0.*(1-exp(-TI./t1vals)) + Meq.*exp(TI./t1vals).*(cosd(a).*exp(-TR./t1vals)).^(n/2-1);

%Finally, solve for SN assuming no T2* weighting
SN = MN.*sind(a);

end