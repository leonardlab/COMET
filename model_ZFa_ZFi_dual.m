%Donahue PS, Draut JW, Muldoon JJ, Edelstein HI, Bagheri N, & Leonard JN.
%The COMET toolkit for composing customizable genetic programs in mammalian cells.

%Simulation of ZFa-inducible and ZFi-inhibitable gene expression.
%    The ZFi effect is represented as both competitive and
%    affecting cooperative RNAPII recruitment.


function sim = model_ZFa_ZFi_dual(doseZFa, doseZFi, m, wA, wI, b, l, u, Z)


%Input arguments:
%    doseZFa: ZFa plasmid (ng)
%    doseZFi: ZFi plasmid (ng)
%    m:       maximum activation parameter
%    wA:      steepness parameter for ZFa
%    wI:      steepness parameter for ZFi
%    b:       background (TF-independent transcription) parameter
%    l:       weight-normalized ratio of inhibitor to activator at which ramp-down begins
%    u:       weight-normalized ratio of inhibitor to activator at which ramp-down ends
%    Z:       one row from the population matrix (not the entire 2-D matrix)

%Output argument:
%    sim: simulation (dimensions: number of time points x number of state variables)

%Notes:
%    Output is in model-specific a.u.
%    The dual mechanism (competitive inhibition and loss of cooperativity)
%        corresponds to three of the four types of promoter-inihibitor
%        pairings.

%Example for ZF1a and the c'th cell in a heterogeneous population:
%    sim = model_ZFa_ZFi_dual(50, 50, 32.7, 0.036, 0.036, 0.08, 0, 1.5, Z(c, :));

%Example for a homogeneous simulation:
%    sim = model_ZFa_ZFi_dual(50, 50, 32.7, 0.036, 0.036, 0.08, 0, 1.5, [1, 1, 1]);


%*****************%
%**** Specify ****%
%*****************%


%timecourse (h)
T = 0:42; %the first number is zero, and the second number corresponds to
%          the time of measurement post-transfection

%general parameters
ktxCMV   = 1;     %transcription                         defined as 1 a.u.
ktxZF    = 9;     %transcription of regulated promoters  fitted
kdegR    = 2.7;   %degradation of RNA                    assumed half-life of 15 min.
ktl      = 1;     %translation                           defined as 1 a.u.
kdegZFP  = 0.35;  %degradation of ZF-TF                  assumed half-life of 2 h.
kdegRepP = 0.029; %degradation of reporter protein       half-life ~1 day.

%initial values (a.u.)
IV = zeros(6, 1);

%intercellular variation for the ZFa
za = Z(1);

%intercellular variation for the ZFi
zi = Z(2);

%intercellular variation for the reporter
zrep = Z(3);


%***************%
%**** Model ****%
%***************%


%simulate ODEs
[~, sim] = ode15s(@(t, y, options)[ 
    
    %Y1 ZFa RNA
      ktxCMV   * za * doseZFa                                          ... %constitutive transcription
    - kdegR    * y(1)                                                      %RNA degradation
    
    %Y2 ZFa protein
      ktl      * y(1)                                                  ... %translation
    - kdegZFP  * y(2)                                                      %protein degradation
    
    %Y3 ZFi RNA
      ktxCMV   * zi * doseZFi                                          ... %constitutive transcription
    - kdegR    * y(3)                                                      %RNA degradation
    
    %Y4 ZFi protein
      ktl      * y(3)                                                  ... %translation
    - kdegZFP  * y(4)                                                      %protein degradation    
    
    %Y5 Reporter RNA
      ktxZF    * zrep * (b + max(min(((wI * y(4)) ./ (wA * y(2)) - l)  ...
          * (1 - m) / (u - l) + m, m), 1) * wA .* y(2))                ...
          ./ (1 + wA * y(2) + wI * y(4))                               ... %inducible and inhibitable transcription
    - kdegR    * y(5)                                                      %RNA degradation
    
    %Y6 Reporter protein
      ktl      * y(5)                                                  ... %translation
    - kdegRepP * y(6)                                                      %protein degradation
    
    ], T, IV);


end

