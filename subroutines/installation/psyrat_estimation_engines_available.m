function avail = psyrat_estimation_engines_available(minver)
%Report which estimation engines PsyRAT can use on this machine.
%
%avail = psyrat_estimation_engines_available()
%avail = psyrat_estimation_engines_available(minver)
%
%Purpose:
% PsyRAT can estimate variance components with three engines: CmdStan (the
% preferred Bayesian engine), a native Hamiltonian Monte Carlo engine built on
% MATLAB's hmcSampler, and a native frequentist engine built on fitlme. Which
% of these can actually run depends on what is installed on the machine
% (CmdStan, the Statistics and Machine Learning Toolbox, and - for the HMC
% engine's automatic-differentiation gradients - the Deep Learning Toolbox).
%
% This function is the SINGLE SOURCE OF TRUTH for that detection. It is reused
% by (a) the processing-preferences GUI dropdown so only runnable engines are
% offered, (b) the automatic engine cascade (CmdStan -> native HMC -> fitlme)
% so it only considers installed engines, and (c) command-line/forced-engine
% validation so a forced choice whose toolbox is missing fails with a clear,
% actionable message. Keeping detection in one place prevents these callers
% from drifting apart.
%
% Detection is conservative and non-destructive: it never launches an
% estimation run, only checks for the presence of installs/toolboxes.
%
%Optional Inputs:
% minver - minimum acceptable CmdStan version as 'major.minor.patch'
%  (default '2.26.0', matching psyrat_require_min_cmdstan; the generated Stan
%  models need the array[N] T syntax introduced in Stan 2.26).
%
%Output:
% avail - struct of logical flags plus human-readable detail strings:
%  .cmdstan - CmdStan is installed and meets the minimum version (or its
%             version could not be read but a runnable install was located,
%             mirroring the fail-open behavior of psyrat_require_min_cmdstan).
%  .cmdstan_status - char reason behind .cmdstan: 'ok', 'absent', 'broken'
%             (stanc present but won't run), 'unreadable' (ran, version
%             unparseable, assumed usable), or 'tooold' (below the minimum). The
%             auto cascade uses 'tooold' to error with an upgrade message instead
%             of silently falling back to an approximate native engine.
%  .hmc     - the native HMC engine can run (Statistics and Machine Learning
%             Toolbox: hmcSampler).
%  .fitlme  - the native fitlme engine can run (Statistics and Machine Learning
%             Toolbox: fitlme).
%  .ad      - automatic differentiation is available for HMC gradients (Deep
%             Learning Toolbox: dlarray). Informational; the HMC engine can
%             also fall back to numerical/analytic gradients.
%  .native  - convenience flag, true if either native engine (.hmc or .fitlme)
%             is available.
%  .details - struct of one-line human-readable strings (.cmdstan, .hmc,
%             .fitlme, .ad) suitable for warnings and dialog text.
%
%See also: psyrat_detect_cmdstan_version, psyrat_require_min_cmdstan

% Copyright (C) 2016-2025 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see
%     <http://www.gnu.org/licenses/>.
%

if nargin < 1 || isempty(minver)
    minver = '2.26.0';
end

avail = struct();

%--- CmdStan (preferred engine) -----------------------------------------------
%Locate the install and read its version. home is non-empty when an install is
%found even if the version string could not be read; ranok is true only when the
%stanc compiler actually ran. cmdstan_status records WHY CmdStan is/ isn't usable
%('ok' | 'absent' | 'broken' | 'unreadable' | 'tooold') so callers (the auto
%cascade) can distinguish a too-old CmdStan from an absent one.
[verstr, home, ranok] = psyrat_detect_cmdstan_version();
if isempty(home)
    avail.cmdstan = false;
    avail.cmdstan_status = 'absent';
    cmdstan_detail = ['CmdStan: no installation found on the MATLAB path or via ' ...
        'the CMDSTAN / STAN_HOME environment variables.'];
elseif ~ranok
    %An install directory was located but stanc could not be run (wrong-arch
    %binary, missing shared libraries, no execute permission). It cannot compile
    %a model, so treat CmdStan as unavailable and let the cascade fall back to a
    %native engine rather than choosing a broken CmdStan that dies at compile time.
    avail.cmdstan = false;
    avail.cmdstan_status = 'broken';
    cmdstan_detail = sprintf(['CmdStan: install found at %s, but its compiler ' ...
        '(stanc) could not be run; treating CmdStan as unavailable.'],home);
elseif isempty(verstr)
    %stanc ran but its version banner could not be read. Mirror the fail-open
    %preflight (psyrat_require_min_cmdstan) and treat CmdStan as usable; the run
    %itself will surface any real problem.
    avail.cmdstan = true;
    avail.cmdstan_status = 'unreadable';
    cmdstan_detail = sprintf(['CmdStan: install found at %s (version could not ' ...
        'be read; assumed usable).'],home);
else
    avail.cmdstan = psyrat_version_ge(verstr,minver);
    if avail.cmdstan
        avail.cmdstan_status = 'ok';
        cmdstan_detail = sprintf('CmdStan %s found at %s.',verstr,home);
    else
        avail.cmdstan_status = 'tooold';
        cmdstan_detail = sprintf(['CmdStan %s found at %s, but %s or newer is ' ...
            'required.'],verstr,home,minver);
    end
end

%--- Statistics and Machine Learning Toolbox (fitlme + hmcSampler) ------------
%Require both a checkout-able license AND the function being present on the
%path, so a licensed-but-not-installed toolbox is reported unavailable.
stats_lic = license('test','Statistics_Toolbox');
avail.fitlme = stats_lic && (exist('fitlme','file') > 0);
avail.hmc    = stats_lic && (exist('hmcSampler','file') > 0);

%--- Deep Learning Toolbox (automatic differentiation for HMC gradients) ------
%The license feature name remains 'Neural_Network_Toolbox' after the product's
%rename to Deep Learning Toolbox.
dl_lic = license('test','Neural_Network_Toolbox');
avail.ad = dl_lic && (exist('dlarray','file') > 0);

%--- convenience flag ---------------------------------------------------------
avail.native = avail.fitlme || avail.hmc;

%--- human-readable details for warnings/dialogs ------------------------------
avail.details = struct();
avail.details.cmdstan = cmdstan_detail;
avail.details.fitlme = local_toolbox_detail(...
    'Statistics and Machine Learning Toolbox (fitlme)',avail.fitlme);
avail.details.hmc = local_toolbox_detail(...
    'Statistics and Machine Learning Toolbox (hmcSampler)',avail.hmc);
avail.details.ad = local_toolbox_detail(...
    'Deep Learning Toolbox (dlarray, automatic-differentiation gradients)',avail.ad);

end

function s = local_toolbox_detail(name,present)
if present
    s = sprintf('%s: available.',name);
else
    s = sprintf('%s: NOT available.',name);
end
end
