%% PV_LIB Toolbox Release Notes
% This toolbox implements functions that enable simulation of the
% performance of photovoltaic (PV) energy systems.
%
% The PV_LIB Toolbox requires Matlab software.
% 
%% References
% Additional information and documentation is available on the PV
% Performance Modeling Collaborative website (<http://pvpmc.org>)
%
% The latest version of this fully functional toolbox is available on the
% PVPMC website. (Check back for periodic updates)
%
%% Bug Reporting
% Report bugs and problmes to Joshua Stein (jsstein@sandia.gov)
%% Credits for Non-Sandia contributions
% * *Rob Andrews, Queen's University:* Multiple bug finds and fixes in
% PV_LIB version 1.0 functions pvl_perez, pvl_haydavies1980,
% pvl_klucher1979, and pvl_reindl1990
% * *Martin Herrerias Azcue:* Found incorrect documentation notes for
% pvl_calcparams_desoto 
% * *John McKeen, DOW Solar Solutions:* Found incorrect error checking in
% pvl_calcparams_desoto.
% * *Mark Campanelli, NREL:* Found bug in pvl_clearsky_ineichen.
%
%% Versions 
%%
% * *Version 1.0:* June-2012  Initial Release
% * *Version 1.1:* December-2012  Update with new functions and bug fixes
%%
% * pvl_clearsky_haurwitz - Clear sky global horizontal irradiance model (simple) 
% * pvl_clearsky_ineichen - Clear sky irradiance (GHI, DNI, and DHI) model
% * pvl_physicaliam - Incident angle modifer model based on Snell’s Law
% * pvl_ashraeiam - Incident angle modifier model from ASHRAE (used in PVsyst)
% * pvl_calcparams_desoto – Create module performance coefficient structure for the single diode model form described by DeSoto et al., 2006
% * pvl_singlediode – Solves the single-diode model to obtain a photovoltaic IV curve 
%
% Significant Changes in Version 1.1
%%
% * Fixed numerous text typos in documentation and help files
% * Made numerical fix to pvl_spa
% * Fixed angle of incidence calculation in pvl_perez, pvl_haydavies1980, pvl_klucher 1979, and pvl_reindl1990
% * Fixed numerical errors in final calculation of pvl_perez
% * Fixed pvl_perez to accept scalar input values with vector input values
%
% * *Version 1.2:* December-2014  Update with new functions and bug fixes
%%
% * pvl_louche - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the LOUCHE model
% * pvl_erbs - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Erbs model
% * pvl_orgill_hollands - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Orgill and Hollands model
% * pvl_reindl_1 - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Reindl 1 model
% * pvl_reindl_2 - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Reindl 2 model
% * pvl_martinruiziam - Determine the incidence angle modifier using the Martin and Ruiz incident angle model
%
% Signifficant Changes in Version 1.2
%%
% * Corrected pvl_klucher1979. Statement  GHI(GHI<DHI) = DHI did not work for vectors, changed to GHI(GHI<DHI) = DHI(GHI<DHI)
% * Corrected pvl_calcparams_desoto.m as fillows:
%   a.	Removed line M = max(M, 0); which set values of M which were less than 0 to 0
%   b.	Removed line S(S==0) = 1E-10; 
%   c.	Added line IL(isnan(M) | M<0 | S <=0) = 0; This sets IL to 0 when
%       i.	Irradiance is <0 or 
%       ii.	Airmass modifier is <0 (most likely due to evaluating the polynomial at very high airmass or 
%       iii.	Airmass modifier is NaN, which is returned by pvl_relativeairmass for sun zenith angles > 90 degrees
%   d.	Added line I0(IL==0) = 0; According the circuit diagram, IL is the only source, and therefore, there can be no reverse saturation current (I0) when IL is 0
%   e.	Added line Rsh(S <= 0) = inf; This is due to the fact that Rsh is determined by dividing by S, negative values would give a negative Rsh, and at S=0 Rsh is undefined
% * 2.	Modified pvl_singlediode.m 
%   a.	Pre-allocate memory to Imax, negPmp, Vmax, Ix, Ixx, Voc, and Isc. This is necessary in order to use the filter (see below)
%   b.	Added a filter u = IL > 0; Thus we only compute IV points and IV curves when there is a photocurrent (IL)
%   c.	Changed the Isc and Voc generation lines to only generate Isc and Voc under conditions which satisfy filter u.
%   d.	Changed Imax and negPmp finding line to only find Imax and negPmp when conditions satisfy filter u
%   e.	Changed Vmax, Ixx, and Ix generation lines to only generate the values under conditions which satisfy filter u
%   f.	Pre-allocate memory for Result.I and Result.V (necessary to implement filter u).
%   g.	Changed Result.I finding line in order to only find I-V curve currents under conditions which satisfy filter u
%   h.	Added line Result.I(:,end) = 0; in order to ensure that the current at Voc is exactly 0 (prevents numerical error which may result in a negative current at Voc). 
% * Corrected errors in wapr_vec.m.  From line 134 through line 154 (inclusive), the Lambert W function is evaluated piecewise over the domain (-exp(-1), inf) using four approximations.  The domain is partitioned into four disjoint sets and each approximation is applied to one (and only one) set.  Filters were implemented incorrectly which allowed later approximations to overwrite values obtained from early approximations.
% * Corrected pvl_calcparams_desoto.m documentation. The documentation used to specify that alpha_isc should be in units of 1/C, it now specifies that alpha_isc should be in units of A/C or A/K. Documentation was changed rather than code due to the fact that the SAM CEC module database lists alpha_isc in units of A/C. 
%   a.	Thanks to Martin Herrerias Azcue for the bug find.
% * Corrected pvl_calcparams_desoto.m. Prior input checking for EgRef was p.addRequired('EgRef', @(x) (isnumeric(x) & isvector(x) & x>0)); this has been changed to p.addRequired('EgRef', @(x) (isnumeric(x) & isvector(x) & all(x>0))); in order to ensure that all EgRef values are positive.
%   a.	Thanks to John McKeen at DOW Solar Solutions for the bug find
% * Corrected pvl_clearsky_ineichen to use system-dependent file separators
% when using the default Linke Turbidity index. Prior code was:  load('Required Data\LinkeTurbidities.mat'); this has been changed to load(['Required Data' filesep 'LinkeTurbidities.mat']); in order to allow for system-specific file separators.
%   a.	Thanks to Mark Campanelli at NREL for finding this bug
% * Modified pvl_singlediode to not require pvl_fminbnd in order to find the maximum power point. Prior to this change, maximum power point was found by finding the maximum value of the power-current curve using a modified version of MATLAB’s fminbnd function. pvl_singlediode now finds where dP/dV as a function of current is equal to 0. Uses bisection techniques.
%   a.	New subfunctions “calc_phi_exact”, “calc_Imp_bisect”, “g”, and “calc_Pmp_bisect” were added.
%   b.	Prior code was left in the function and commented out to allow for reversion later if necessary.
% * 1.	Modified pvl_getaoi to avoid complex values. Changed the line: 
%   AOI = acosd(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz)); 
%   to the following:
%   AOI = acosd(max(min(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz), 1),-1)); 
%   a.	Roundoff errors previously could cause the argument of the arcos function to be greater than 1 or less than -1 (resulting in a complex output). The min and max functions prevent such an occurrence.
% * Modified pvl_getaoi to avoid complex values. Changed the line: 
%   AOI = acosd(max(min(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz), 1),-1));
%   to the following:
%   temp = cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz); 
%   temp(temp>1) = 1; temp(temp<-1) = -1; 
%   AOI = acosd(temp); 
%   a.	Roundoff errors previously could cause the argument of the arcos function to be greater than 1 or less than -1 (resulting in a complex output). The min and max functions did not prevent all occurrences.

%%
% Copyright 2014 Sandia National Laboratories