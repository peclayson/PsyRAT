function [template,titleText] = psyrat_formula_template(varargin)
%Return a runnable command template for a catalog formula.
%
%[template,titleText] = psyrat_formula_template('formula_id','t2_global_d')

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

formula_id = '';

if ~isempty(varargin)
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n'));
    end
    
    ind = find(strcmpi('formula_id',varargin),1);
    if ~isempty(ind)
        formula_id = char(string(varargin{ind+1}));
    end
end

if isempty(formula_id)
    error('varargin:formula_id',... %Error code and associated error
        'WARNING: formula_id is required.');
end

switch lower(formula_id)
    case 't2_global_d'
        titleText = 'Table 2 | Global Coefficient | Absolute';
        template = ['[ll,pt,ul] = psyrat_rel_sing(''',...
            'gcoeff'',1,''metric'',''global'',''bp'',0.50,''wp'',1.00,''obs'',20,''CI'',0.95);'];
    case 't2_global_g'
        titleText = 'Table 2 | Global Coefficient | Relative';
        template = ['[ll,pt,ul] = psyrat_rel_sing(''',...
            'gcoeff'',2,''metric'',''global'',''bp'',0.50,''wp'',1.00,''i'',0.20,''obs'',20,''CI'',0.95);'];
    case 't2_icc_abs'
        titleText = 'Table 2 | ICC | Absolute';
        template = ['[ll,pt,ul] = psyrat_rel_sing(''',...
            'gcoeff'',1,''metric'',''icc'',''bp'',0.50,''wp'',1.00,''CI'',0.95);'];
    case 't2_icc_rel'
        titleText = 'Table 2 | ICC | Relative';
        template = ['[ll,pt,ul] = psyrat_rel_sing(''',...
            'gcoeff'',2,''metric'',''icc'',''bp'',0.50,''wp'',1.00,''i'',0.20,''CI'',0.95);'];
    case 't2_cut_d'
        titleText = 'Table 2 | Cut-Score | Absolute';
        template = ['[ll,pt,ul] = psyrat_cutscore_sing(''',...
            'gcoeff'',1,''bp'',0.50,''wp'',1.00,''mu'',0.10,''cut'',0,''obs'',20,''CI'',0.95);'];
    case 't3_rel_d'
        titleText = 'Table 3 | CE/CS/CES | Absolute';
        template = ['[ll,pt,ul] = psyrat_rel_trt(''',...
            'gcoeff'',1,''reltype'',1,''bp'',0.50,''bo'',0.10,''bt'',0.20,''txp'',0.15,''oxp'',0.12,''txo'',0.09,''err'',0.70,''obs'',20,''nocc'',2,''CI'',0.95);'];
    case 't3_rel_g'
        titleText = 'Table 3 | CE/CS/CES | Relative';
        template = ['[ll,pt,ul] = psyrat_rel_trt(''',...
            'gcoeff'',2,''reltype'',1,''bp'',0.50,''bo'',0.10,''bt'',0.20,''txp'',0.15,''oxp'',0.12,''txo'',0.09,''err'',0.70,''obs'',20,''nocc'',2,''CI'',0.95);'];
    case 't3_icc_abs'
        titleText = 'Table 3 | ICC (CE/CS/CES) | Absolute';
        template = ['[ll,pt,ul] = psyrat_rel_trt(''',...
            'gcoeff'',1,''reltype'',1,''bp'',0.50,''bo'',0.10,''bt'',0.20,''txp'',0.15,''oxp'',0.12,''txo'',0.09,''err'',0.70,''obs'',1,''nocc'',1,''CI'',0.95);'];
    case 't3_icc_rel'
        titleText = 'Table 3 | ICC (CE/CS/CES) | Relative';
        template = ['[ll,pt,ul] = psyrat_rel_trt(''',...
            'gcoeff'',2,''reltype'',1,''bp'',0.50,''bo'',0.10,''bt'',0.20,''txp'',0.15,''oxp'',0.12,''txo'',0.09,''err'',0.70,''obs'',1,''nocc'',1,''CI'',0.95);'];
    case 't3_cut_d'
        titleText = 'Table 3 | Cut-Score | Absolute';
        template = ['[ll,pt,ul] = psyrat_cutscore_trt(''',...
            'gcoeff'',1,''reltype'',1,''bp'',0.50,''bo'',0.10,''bt'',0.20,''txp'',0.15,''oxp'',0.12,''txo'',0.09,''err'',0.70,''mu'',0.10,''cut'',0,''obs'',20,''nocc'',2,''CI'',0.95);'];
    case 't6_diff1_d'
        titleText = 'Table 6 | One-Facet Difference Reliability | Absolute';
        template = ['bp = reshape([0.50 0.10 0.10 0.45],[1 2 2]); bt = reshape([0.25 0.04 0.04 0.20],[1 2 2]); ',...
            'er_var = log([sqrt(0.55) sqrt(0.62)]); out = psyrat_diffrel(''',...
            'bp'',bp,''bt'',bt,''er_var'',er_var,''wp_cov'',0.02,''obs'',[20 20],''est'',''dep'',''CI'',0.95);'];
    case 't6_diff1_g'
        titleText = 'Table 6 | One-Facet Difference Reliability | Relative';
        template = ['bp = reshape([0.50 0.10 0.10 0.45],[1 2 2]); bt = reshape([0.25 0.04 0.04 0.20],[1 2 2]); ',...
            'er_var = log([sqrt(0.55) sqrt(0.62)]); out = psyrat_diffrel(''',...
            'bp'',bp,''bt'',bt,''er_var'',er_var,''wp_cov'',0.02,''obs'',[20 20],''est'',''gen'',''CI'',0.95);'];
    case 't6_diff2_d'
        titleText = 'Table 6 | Two-Facet Difference Reliability | Absolute';
        template = ['bp=reshape([0.50 0.10 0.10 0.45],[1 2 2]); bpi=reshape([0.25 0.06 0.06 0.22],[1 2 2]); ',...
            'bpo=reshape([0.20 0.05 0.05 0.18],[1 2 2]); bt=reshape([0.18 0.04 0.04 0.16],[1 2 2]); ',...
            'bo=reshape([0.10 0.02 0.02 0.09],[1 2 2]); boi=reshape([0.08 0.02 0.02 0.07],[1 2 2]); ',...
            'er_var=log([sqrt(0.50) sqrt(0.55)]); out=psyrat_diffrel_trt(''',...
            'bp'',bp,''bpi'',bpi,''bpo'',bpo,''bt'',bt,''bo'',bo,''boi'',boi,''er_var'',er_var,''er_cov'',0.03,''obs'',[20 20],''nocc'',[2 2],''reltype'',3,''est'',''dep'',''CI'',0.95);'];
    case 't6_diff2_g'
        titleText = 'Table 6 | Two-Facet Difference Reliability | Relative';
        template = ['bp=reshape([0.50 0.10 0.10 0.45],[1 2 2]); bpi=reshape([0.25 0.06 0.06 0.22],[1 2 2]); ',...
            'bpo=reshape([0.20 0.05 0.05 0.18],[1 2 2]); bt=reshape([0.18 0.04 0.04 0.16],[1 2 2]); ',...
            'bo=reshape([0.10 0.02 0.02 0.09],[1 2 2]); boi=reshape([0.08 0.02 0.02 0.07],[1 2 2]); ',...
            'er_var=log([sqrt(0.50) sqrt(0.55)]); out=psyrat_diffrel_trt(''',...
            'bp'',bp,''bpi'',bpi,''bpo'',bpo,''bt'',bt,''bo'',bo,''boi'',boi,''er_var'',er_var,''er_cov'',0.03,''obs'',[20 20],''nocc'',[2 2],''reltype'',3,''est'',''gen'',''CI'',0.95);'];
    case 'dyn_sigma_z'
        titleText = 'Rast & Clayson | Conditional residual SD | sigma_pi,e(z)';
        template = ['sige0 = 3; b_sigma = 0.3; z1 = [-1.5 0 1.5]; ',...
            'sigma_e_z = exp(log(sige0) + b_sigma .* z1);'];
    case 'dyn_g'
        titleText = 'Rast & Clayson | Dynamic Generalizability | G(z)';
        template = ['out = psyrat_rel_dynrel(''',...
            'sig_p'',4,''sig_log0'',log(3),''b_sigma'',0.3,''sig_i'',2,''ndim'',1,',...
            '''z1'',[-1.5 0 1.5],''obs'',20,''CI'',0.95); Gz = out.G.pt;'];
    case 'dyn_d'
        titleText = 'Rast & Clayson | Dynamic Dependability | D(z)';
        template = ['out = psyrat_rel_dynrel(''',...
            'sig_p'',4,''sig_log0'',log(3),''b_sigma'',0.3,''sig_i'',2,''ndim'',1,',...
            '''z1'',[-1.5 0 1.5],''obs'',20,''CI'',0.95); Dz = out.D.pt;'];
    case 'dyn_icc_rel'
        titleText = 'Rast & Clayson | Conditional ICC | Relative';
        template = ['out = psyrat_rel_dynrel(''',...
            'sig_p'',4,''sig_log0'',log(3),''b_sigma'',0.3,''sig_i'',2,''ndim'',1,',...
            '''z1'',[-1.5 0 1.5],''obs'',1,''CI'',0.95); ICCg = out.ICCg.pt;'];
    case 'dyn_icc_abs'
        titleText = 'Rast & Clayson | Conditional ICC | Absolute';
        template = ['out = psyrat_rel_dynrel(''',...
            'sig_p'',4,''sig_log0'',log(3),''b_sigma'',0.3,''sig_i'',2,''ndim'',1,',...
            '''z1'',[-1.5 0 1.5],''obs'',1,''CI'',0.95); ICCd = out.ICCd.pt;'];
    case 'dyn_diff_g'
        titleText = 'Rast & Clayson | Dynamic difference reliability | G(z)';
        template = ['% Dynamic difference reliability is estimated from single-trial data ',...
            'via Specify Inputs (difference scores + Dimension 1). The dimension-conditioned ',...
            'residual sigma_pi,e(z) feeds psyrat_diffrel per z; see psyrat_dynrel_summary.'];
    case 'dyn_diff_d'
        titleText = 'Rast & Clayson | Dynamic difference reliability | D(z)';
        template = ['% As dyn_diff_g, but the absolute (dependability) form adds the ',...
            'difference trial main effect (and, with occasions, the occasion terms) to the ',...
            'error. Estimated via Specify Inputs (difference scores + Dimension 1).'];
    case 'dyn_rescor'
        titleText = 'Rast & Clayson | Concurrent residual correlation | rho_e,12 / SD(rho_e,s)';
        template = ['% Concurrent residual correlation for co-occurring events. Population ',...
            '(shared) estimates one rho_e,12 across participants (cases 17-20); Per-subject ',...
            'adds the between-person SD(rho_e,s) (cases 21-22). Under the gamma modular cut ',...
            'designs (19/20) the reported quantity is rho_e,12(copula), the LATENT-normal ',...
            'copula correlation - a different estimand from the Gaussian observed-scale ',...
            'rescor. Estimated via Specify Inputs ',...
            '(difference scores + subject-specific error variances + within-person residual covariance).'];
    case 't4_d'
        titleText = 'Rocha Table 4 | Splits one-facet | Absolute (D)';
        template = ['% Nonparallel splits, single occasion (analysis 23, p x (i:s)). Estimated ',...
            'data-driven via Specify Inputs: one row per split, the required n_i (weight) column, ',...
            'and Nonparallel. psyrat_splits_summary reuses psyrat_rel_sing with an SD-scaled residual ',...
            '(nbar_i = harmonic mean of n_i). Dependability (absolute) is exact.'];
    case 't4_g'
        titleText = 'Rocha Table 4 | Splits one-facet | Relative (G)';
        template = ['% As t4_d but the relative coefficient. Generalizability is a downward-biased ',...
            'LOWER BOUND: split means lump the item main effect sigma_(i:s) into the residual, so G ',...
            'is understated and gbias-flagged in the output. Estimated via Specify Inputs (n_i column, Nonparallel).'];
    case 't4_sem'
        titleText = 'Rocha Table 4 | Splits one-facet | SEM';
        template = ['% Standard error of measurement for the single-occasion splits design. ',...
            'Absolute SEM = sqrt((sig_ps^2 + sig_s^2 + sig_err^2/nbar_i)/ns) is exact; the relative ',...
            'SEM (drops sig_s^2) is an upper bound. Reported per stratum by psyrat_splits_summary.'];
    case 't5_d'
        titleText = 'Rocha Table 5 | Splits two-facet (test-retest) | Absolute (D)';
        template = ['% Nonparallel splits, multiple occasions (analysis 24, p x (i:s) x o). Estimated ',...
            'via Specify Inputs (n_i column, Nonparallel, plus a Time variable). psyrat_splits_summary ',...
            'reuses psyrat_rel_trt for all three reltypes (CE/CS/CES). Dependability (absolute) is exact ',...
            'under trials-nested-in-occasions.'];
    case 't5_g'
        titleText = 'Rocha Table 5 | Splits two-facet (test-retest) | Relative (G)';
        template = ['% As t5_d but the relative coefficient (CE/CS/CES). Generalizability is a ',...
            'downward-biased LOWER BOUND (the residual lumps sigma_(i:s)); gbias-flagged in the output.'];
    case 't5_sem'
        titleText = 'Rocha Table 5 | Splits two-facet (test-retest) | SEM';
        template = ['% SEM for the multi-occasion splits design, per reltype. The error term uses ',...
            'err = sqrt(sig_pos^2 + sig_err^2/nbar_i). Absolute SEM is exact; the relative SEM is an ',...
            'upper bound. Reported per stratum/reltype by psyrat_splits_summary.'];
    otherwise
        error('varargin:formula_id',... %Error code and associated error
            strcat('WARNING: Unknown formula_id: ',formula_id,'\n'));
end

end
