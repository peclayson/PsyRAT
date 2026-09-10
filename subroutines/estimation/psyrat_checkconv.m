function RELout = psyrat_checkconv(REL)
%Check for convergence by examining the potential scale reduction factor
%(r_hat) and the effective sample size (n_eff)
%
%[REL, rerun] = psyrat_checkconv(REL)
%
%Lasted Updated 1/19/17
%
%Required Input:
% REL - structure array created by psyrat_computevarcomp
%
%Outputs:
% RELout - structure array with the following added fields
%  out.conv.converged - 1 data converged, 0 data did not converge

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

%nargin is a scalar, so the length(nargin) this guard used to test was always 1
%and the guard never fired. The signature takes exactly one positional input, so
%MATLAB rejects two or more at the call itself; nargin ~= 1 therefore means the
%caller passed none, which previously surfaced as an opaque "Undefined variable
%REL" from the size() call below instead of the message written here.
if nargin ~= 1
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Input not specified \n\n',... 
    'See help psyrat_checkconv for more information on  inputs'));
end

%pull the dimensions to figure out if there are nested cell arrays
%if there are there will be 1 row and multiple columns (a = 1)
[a,nloops] = size(REL.out.conv.data);

%FAIL CLOSED (audit 2026-09-05). Two edge states used to certify convergence:
%(1) a NaN R-hat or n_eff, which stansummary prints for a parameter with zero
%variance across draws (a chain stuck at a constant is the textbook pathology
%this hides) and str2double returns as NaN -- and NaN >= 1.1 and NaN < floor
%are both false, so the parameter was never flagged; (2) a convergence table
%with no parameter rows at all (no summary row matched the design's prefix
%list), for which find([] ...) is empty and the run read as converged. Both
%now read as NOT converged: the comparisons below are spelled as the negation
%of the pass condition so NaN lands on the failing side, a table with fewer
%than two rows (header only) fails, and a result with no tables at all fails.
%Finite values behave exactly as before, so every ordinary run is unchanged.
%REFINED 2026-09-07 (owner ruling, at the tutorial refit): a row whose n_eff AND
%R-hat are both NaN is a parameter that is constant by construction, not a stuck
%chain, and is excluded from the gate; see local_varying_rows below.
result = 0;

%if there is more than 1 row then that means there is only one set of
%convergence statistics to examine
if a == 1

    for i=1:nloops
        %pull data from REL
        %first check r_hats
        data = REL.out.conv.data{i};
        if size(data,1) < 2
            %header only: nothing was monitored, so nothing was shown to converge
            result = 0;
            break;
        end
        rhatv = [data{2:end,3}];
        neffv = [data{2:end,2}];
        varying = local_varying_rows(neffv, rhatv);
        if ~any(varying)
            %every monitored row is a constant by construction: nothing was
            %shown to converge (the same reading as the header-only table)
            result = 0;
            break;
        end

        %see if any of the rhats didn't equal 1 (i.e., did not converge);
        %spelled as ~(passes) so a NaN R-hat on a varying row is flagged
        indbad = find(~(rhatv(varying) < 1.1));

        %specify whether convergenece between chains was reached
        if ~isempty(indbad)
            result = 0;
            break;
        else
            result = 1;
        end

        if result == 1
            %check whether neff is (5*2*nchains)
            %see if any n_eff fell below the floor (i.e., did not converge);
            %spelled as ~(passes) so a NaN n_eff on a varying row is flagged
            indbad = find(~(neffv(varying) >= (5*2*REL.nchains)));

            %specify whethere convergenece between chains was reached
            if ~isempty(indbad)
                result = 0;
                break;
            else
                result = 1;
            end
        end
    end
    
elseif a ~= 1
    
    %pull data from REL
    %first check r_hats
    rhatv = [REL.out.conv.data{2:end,3}];
    neffv = [REL.out.conv.data{2:end,2}];
    varying = local_varying_rows(neffv, rhatv);

    %see if any of the rhats didn't equal 1 (i.e., did not converge);
    %spelled as ~(passes) so a NaN R-hat on a varying row is flagged; a
    %header-only table (no monitored parameter) and an all-constant table (no
    %varying row) fail rather than pass
    indbad = find(~(rhatv(varying) < 1.1));

    %specify whether convergenece between chains was reached
    if ~isempty(indbad) || size(REL.out.conv.data,1) < 2 || ~any(varying)
        result = 0;
    else
        result = 1;
    end

    if result == 1
        %check whether neff is (5*2*nchains)
        %see if any n_eff fell below the floor (i.e., did not converge);
        %spelled as ~(passes) so a NaN n_eff on a varying row is flagged
        indbad = find(~(neffv(varying) >= (5*2*REL.nchains)));

        %specify whethere convergenece between chains was reached
        if ~isempty(indbad)
            result = 0;
        else
            result = 1;
        end
    end
end
    
    
REL.out.conv.converged = result;
RELout = REL;

end

function varying = local_varying_rows(neffv, rhatv)
%Which monitored rows can carry a convergence verdict. stansummary prints NaN for
%BOTH n_eff and R-hat when a parameter never moves in any chain, and the monitored
%set contains such parameters by design: correlation-matrix diagonals fixed at 1,
%and disabled covariance / residual-correlation placeholders fixed at 0. The
%2026-09-05 fail-closed rule flagged them, so every design that monitors one (most
%designs beyond the plain one-facet and test-retest models) could never converge,
%and the tutorial refit broke on chapter 15's subject-level difference fit after
%four automatic reruns (owner ruling 2026-09-07). Rows with NaN in both columns
%are excluded here. A NaN in only one column is not the constant signature and is
%still flagged by the callers; chains stuck at DIFFERENT constants give an Inf
%R-hat, which ~(Inf < 1.1) still flags.
varying = ~(isnan(neffv) & isnan(rhatv));
end

