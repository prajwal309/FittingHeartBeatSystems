J/MNRAS/489/1644     Modelling Kepler eclipsing binaries      (Windemuth+, 2019)
================================================================================
Modelling Kepler eclipsing binaries:
homogeneous inference of orbital and stellar properties.
    Windemuth D., Agol E., Ali A., Kiefer F.
   <Mon. Not. R. Astron. Soc., 489, 1644-1666 (2019)>
   =2019MNRAS.489.1644W    (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Binaries, eclipsing ; Binaries, orbits ; Stars, masses ;
              Stars, diameters ; Models ; Optical
Keywords: binaries: eclipsing - binaries: close - methods: data analysis

Abstract:
    We report on the properties of eclipsing binaries (EBs) from the
    Kepler mission with a newly developed photometric modelling code,
    which uses the light curve, spectral energy distribution of each
    binary, and stellar evolution models to infer stellar masses without
    the need for radial velocity (RV) measurements. We present solutions
    and posteriors to orbital and stellar parameters for 728 systems,
    forming the largest homogeneous catalogue of full Kepler binary
    parameter estimates to date. Using comparisons to published RV
    measurements, we demonstrate that the inferred properties (e.g.
    masses) are reliable for well-detached main-sequence (MS) binaries,
    which make up the majority of our sample. The fidelity of our inferred
    parameters degrades for a subset of systems not well described by
    input isochrones, such as short-period binaries that have undergone
    interactions, or binaries with post-MS components. Additionally, we
    identify 35 new systems which show evidence of eclipse timing
    variations, perhaps from apsidal motion due to binary tides or
    tertiary companions. We plan to subsequently use these models to
    search for and constrain the presence of circumbinary planets in
    Kepler EB systems.

Description:
    We developed a python code, dubbed 'KEBLAT', to simultaneously fit the
    LCs and SEDs of 728 Kepler EBs to obtain their orbital and stellar
    parameters. The open-source code is available on Github
    (https://github.com/savvytruffle/keblat/).

    We selected our binary systems from the VKEBC (Prsa et al.
    2011AJ....141...83P, Cat. J/AJ/141/83; Kirk et al.
    2016AJ....151...68K, Cat. J/AJ/151/68). This sample, totaling 2877
    targets, has been compiled from the entire Kepler prime mission and
    includes binaries which are not eclipsing (e.g. ellipsoidal
    variables). We made selection cuts to reduce grazing (i.e.
    non-constraining) geometries and the effects of close binarity that
    are not included in our model.

    In total, we investigated 728 ~detached EB systems. The large
    reduction in the number of systems from the full VKEBC is related to
    the fact that short-period binaries have the highest eclipse
    probability, but also have the largest tidal distortions; our sample
    emphasizes the longer period EBs.

    We extracted LCs from Kepler long-cadence (~30min) simple aperture
    photometry (SAP) based on Data Release (DR) 24 (Thompson et al. 2015,
    Kepler Data Release 24 Notes (KSCI-19064-002)), spanning ~4yr or 17
    quarters. The LC module couples a Keplerian orbit solver to the
    analytic Mandel & Agol (2002ApJ...580L.171M) transit model with a
    quadratic limb-darkening law to fit the observed LC.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
table2.dat       328      728   EB orbital and distance parameter posteriors
table3.dat       384      728   EB stellar parameter posteriors
table4.dat       779      728   ML parameter solutions
table5.dat        38       84   ETV candidates identified in our sample
table6.dat        42        3   Mass ratios from APOGEE RVs
table7.dat       128       55   RV-derived EB mass values from literature
--------------------------------------------------------------------------------

See also:
  J/AJ/142/160 : Kepler Mission. II. Eclipsing binaries in DR2 (Slawson+, 2011)
   J/AJ/151/68 : Kepler Mission. VII. Eclipsing binaries in DR3 (Kirk+, 2016)

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8 I8     ---     ID        KIC name (NNNNNNNN)
  10- 26 F17.12 d       P         Orbital period
  28- 44 E17.15 d     E_P         Upper error on P (1{sigma})
  46- 62 E17.15 d     e_P         Lower error on P (1{sigma})
  64- 77 F14.10 d       tPE       Time of primary eclipse (BJD-2454833)
  79- 95 E17.15 d     E_tPE       Upper error on tPE (1{sigma})
  97-113 E17.15 d     e_tPE       Lower error on tPE (1{sigma})
 115-132 E18.15 ---     esinw     Transformation of eccentricity and longitude
                                   of periastron (radial component)
 134-150 E17.15 ---   E_esinw     Upper error on esinw (1{sigma})
 152-168 E17.15 ---   e_esinw     Lower error on esinw (1{sigma})
 170-187 E18.15 ---     ecosw     Transformation of eccentricity and longitude
                                   of periastron (tangential component)
 189-205 E17.15 ---   E_ecosw     Upper error on ecosw (1{sigma})
 207-223 E17.15 ---   e_ecosw     Lower error on ecosw (1{sigma})
 225-237 F13.11 rad     i         Orbit inclination
 239-255 E17.15 rad   E_i         Upper error on i (1{sigma})
 257-273 E17.15 rad   e_i         Lower error on i (1{sigma})
 275-290 F16.10 pc      d         Distance to system
 292-310 F19.13 pc    E_d         Upper error on d (1{sigma})
 312-328 F17.12 pc    e_d         Lower error on d (1{sigma})
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8 I8     ---     ID        KIC name (NNNNNNNN)
  10- 25 F16.14 ---     Z         Binary metallicity defined as 1-X-Y
  27- 44 E18.15 ---   E_Z         Upper error on Z (1{sigma})
  46- 62 E17.15 ---   e_Z         Lower error on Z (1{sigma})
  64- 77 F14.11 [yr]    tau       Logarithm of binary age
  79- 96 E18.15 [yr]  E_tau       Upper error on tau (1{sigma})
  98-114 E17.15 [yr]  e_tau       Lower error on tau (1{sigma})
 116-132 F17.14 Msun    M1        Mass of the first component
 134-151 E18.15 Msun  E_M1        Upper error on M1 (1{sigma})
 153-169 E17.15 Msun  e_M1        Lower error on M1 (1{sigma})
 171-185 F15.13 Msun    M2        Mass of the second component
 187-204 E18.15 Msun  E_M2        Upper error on M2 (1{sigma})
 206-222 E17.15 Msun  e_M2        Lower error on M2 (1{sigma})
 224-238 F15.12 Rsun    R1        Radius of the first component
 240-257 E18.15 Rsun  E_R1        Upper error on R1 (1{sigma})
 259-275 E17.15 Rsun  e_R1        Lower error on R1 (1{sigma})
 277-291 F15.12 Rsun    R2        Radius of the second component
 293-310 E18.15 Rsun  E_R2        Upper error on R2 (1{sigma})
 312-328 E17.15 Rsun  e_R2        Lower error on R2 (1{sigma})
 330-347 F18.15 ---     F2/F1     Kepler flux ratio
 349-366 E18.15 ---   E_F2/F1     Upper error on F2/F1 (1{sigma})
 368-384 E17.15 ---   e_F2/F1     Lower error on F2/F1 (1{sigma})
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table4.dat
--------------------------------------------------------------------------------
   Bytes Format Units Label       Explanations
--------------------------------------------------------------------------------
   1-  8 I8     ---   ID          KIC identifier (NNNNNNNN)
  10- 28 F19.16 Msun  Msum        ?=0.0 Sum of masses M1+M2
  30- 48 F19.17 ---   Q           ?=0.0 Ratio of masses M2/M1
  50- 70 F21.19 ---   Z           ?=0.0 Binary metallicity defined as 1-X-Y
  72- 90 F19.16 [yr]  tau         ?=0.0 Logarithm of binary age
  92-111 F20.14 pc    d           ?=0.0 Distance to system
 113-133 F21.19 ---   E(B-V)      ?=0.0 Reddening assuming R_V_=3.1
 135-139 F5.1   pc    h0          ?=0.0 Dust vertical scale height (fixed at
                                   119)
 141-161 F21.16 d     P           ?=0.0 Orbital period
 163-180 F18.14 d     tPE         ?=0.0 Time of primary eclipse (BJD-2454833)
 182-204 E23.20 ---   esinw       ?=0.0 Transformation of eccentricity and
                                   longitude of periastron (radial component)
 206-228 E23.20 ---   ecosw       ?=0.0 Transformation of eccentricity and
                                   longitude of periastron (tangential
                                   component)
 230-251 E22.19 ---   b           ?=0.0 Impact parameter acosi/R1
 253-273 E21.19 ---   q11         ?=0.0 Transformed quadratic limb-darkening q11
 275-295 E21.19 ---   q12         ?=0.0 Transformed quadratic limb-darkening q12
 297-318 E22.20 ---   q21         ?=0.0 Transformed quadratic limb-darkening q21
 320-340 E21.19 ---   q22         ?=0.0 Transformed quadratic limb-darkening q22
 342-360 F19.16 ---   lnsigLC     ?=0.0 Systematic LC error
 362-384 F23.19 ---   lnsigSED    ?=0.0 Systematic SED error
 386-406 F21.17 ---   lnsigE(B-V) ?=0.0 Uncertainty in reddening prior
 408-426 F19.16 ---   lnsigd      ?=0.0 Uncertainty in Gaia distance prior
 428-445 F18.16 ---   c0          ?=0.0 Crowding parameter for quarter 0
 447-464 F18.16 ---   c1          ?=0.0 Crowding parameter for quarter 1
 466-483 F18.16 ---   c2          ?=0.0 Crowding parameter for quarter 2
 485-503 F19.17 ---   c3          ?=0.0 Crowding parameter for quarter 3
 505-522 F18.16 ---   c4          ?=0.0 Crowding parameter for quarter 4
 524-541 F18.16 ---   c5          ?=0.0 Crowding parameter for quarter 5
 543-560 F18.16 ---   c6          ?=0.0 Crowding parameter for quarter 6
 562-580 F19.17 ---   c7          ?=0.0 Crowding parameter for quarter 7
 582-600 F19.17 ---   c8          ?=0.0 Crowding parameter for quarter 8
 602-619 F18.16 ---   c9          ?=0.0 Crowding parameter for quarter 9
 621-638 F18.16 ---   c10         ?=0.0 Crowding parameter for quarter 10
 640-658 F19.17 ---   c11         ?=0.0 Crowding parameter for quarter 11
 660-677 F18.16 ---   c12         ?=0.0 Crowding parameter for quarter 12
 679-696 F18.16 ---   c13         ?=0.0 Crowding parameter for quarter 13
 698-715 F18.16 ---   c14         ?=0.0 Crowding parameter for quarter 14
 717-735 F19.17 ---   c15         ?=0.0 Crowding parameter for quarter 15
 737-754 F18.16 ---   c16         ?=0.0 Crowding parameter for quarter 16
 756-774 F19.17 ---   c17         ?=0.0 Crowding parameter for quarter 17
 776-779 F4.2   ---   morph       ?=0.0 Morphology value taken directly from
                                   the VKEBC (Prsa et al. 2011AJ....141...83P,
                                   Cat. J/AJ/141/83;
                                   Kirk et al.  2016AJ....151...68K,
                                   Cat. J/AJ/151/68)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8  I8    ---     ID        KIC identifier (NNNNNNNN)
  10- 21  F12.8 d       P         Orbital period
  23- 26  F4.2  ---     Q         Ratio of masses M2/M1
  28- 34  E7.5  ---     e         Eccentricity
  36- 38  A3    ---     refs      [brco] References (1)
--------------------------------------------------------------------------------
Note (1): References as follows:
   b = Borkovits et al. (2016MNRAS.455.4136B, Cat. J/MNRAS/455/4136)
   r = Rappaport et al. (2013ApJ...768...33R)
   c = Conroy et al. (2014AJ....147...45C, Cat. J/AJ/147/45)
   o = Orosz (2015ASPC..496...55O)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table6.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  7  I7    ---     ID        KIC identifier (NNNNNNN)
   9- 17  F9.6  d       P         Orbital period
  19- 26  F8.6  d     e_P         Error on P
  28- 34  F7.5  ---     Q         Ratio of masses M2/M1
  36- 42  F7.5  ---   e_Q         Error on Q
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table7.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8  I8    ---     ID        KIC identifier (NNNNNNNN)
  10- 22  F13.9 d       P         Orbital period
  24- 29  F6.4  Msun    M1        ? Mass of the first component (1)
  31- 36  F6.4  Msun  e_M1        ? Error on M1
  38- 44  F7.5  Msun    M2        ? Mass of the second component (1)
  46- 52  F7.5  Msun  e_M2        ? Error on M2
  54-128  A75   ---     ref       References
--------------------------------------------------------------------------------
Note (1): Note for binary studies which contain multiple mass estimates (e.g.
          from asteroseismology and RV), we adopted the RV-only derived values
          for consistency
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                                           Ana Fiallos [CDS]    05-Jan-2023