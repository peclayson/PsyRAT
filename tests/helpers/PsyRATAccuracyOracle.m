classdef PsyRATAccuracyOracle
    % REGRESSION MIRROR of the production calculation formulas.
    %
    % This class was previously described as holding "independent closed-form oracles".
    % It does not, and the distinction decides what a green test run is allowed to prove.
    % The expressions here are the SAME ALGEBRA as the production functions -- in several
    % cases character-identical -- so agreement demonstrates that the formulas are
    % UNCHANGED, never that they are CORRECT. Both sides would have to be wrong in the
    % same way for a mirror to notice, which is exactly what duplicating a formula
    % guarantees. A comment in this file's own range-layout helper says as much.
    %
    % Two further limits worth knowing before citing a green run:
    %   * dep() mirrors psyrat_dep, which has NO production call sites. The live one-facet
    %     coefficient is psyrat_rel_sing, and this class has no entry for it -- so two
    %     planted one-facet defects have been shown to pass this suite 8/8.
    %   * cutscoreDep() likewise targets psyrat_cutscore_dep, which is also unreachable
    %     from production; the shipped cut-score is psyrat_cutscore_sing.
    %
    % For EXTERNAL corroboration -- formulas transcribed from the published tables by a
    % reader with no access to this source -- see PsyRATIndependentOracle and
    % TestIndependentOracleAgreement. Keep both: a mutation sweep found each catches
    % defects the other misses.

    methods (Static)
        function [ll,pt,ul] = dep(bp,wp,obs,ci)
            [ll,pt,ul] = PsyRATAccuracyOracle.rangeStats(@(n) ...
                (bp(:).^2) ./ ((bp(:).^2) + (wp(:).^2 ./ n)), obs, ci);
        end

        function [ll,pt,ul] = icc(bp,wp,ci)
            icc = (bp(:).^2) ./ ((bp(:).^2) + (wp(:).^2));
            [ll,pt,ul] = PsyRATAccuracyOracle.interval(icc,ci);
        end

        function [ll,pt,ul] = relTrt(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,obs,nocc,ci)
            bps = bp(:).^2;
            bos = bo(:).^2;
            bts = bt(:).^2;
            txps = txp(:).^2;
            oxps = oxp(:).^2;
            txos = txo(:).^2;
            errs = err(:).^2;

            [ll,pt,ul] = PsyRATAccuracyOracle.rangeStats(@(n) ...
                PsyRATAccuracyOracle.relTrtDraws(gcoeff,reltype,bps,bos,bts,txps,oxps,txos,errs,n,nocc), ...
                obs, ci);
        end

        function [ll,pt,ul] = cutscoreDep(bp,wp,mu,cut,obs,ci,est,i)
            if nargin < 8 || isempty(i)
                i = 0;
            end
            if nargin < 7 || isempty(est)
                est = 'dep';
            end

            bps = bp(:).^2;
            wps = wp(:).^2;
            iv = i(:).^2;
            muv = mu(:);

            ndraw = length(bps);
            wps = PsyRATAccuracyOracle.expandDraw(wps,ndraw,'wp');
            iv = PsyRATAccuracyOracle.expandDraw(iv,ndraw,'i');
            muv = PsyRATAccuracyOracle.expandDraw(muv,ndraw,'mu');

            [ll,pt,ul] = PsyRATAccuracyOracle.rangeStats(@(n) ...
                PsyRATAccuracyOracle.cutscoreDepDraws(bps,wps,iv,muv,cut,n,est), obs, ci);
        end

        function [ll,pt,ul] = cutscoreTrt(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,mu,cut,obs,nocc,ci)
            bps = bp(:).^2;
            bos = bo(:).^2;
            bts = bt(:).^2;
            txps = txp(:).^2;
            oxps = oxp(:).^2;
            txos = txo(:).^2;
            errs = err(:).^2;
            muv = mu(:);

            ndraw = length(bps);
            bos = PsyRATAccuracyOracle.expandDraw(bos,ndraw,'bo');
            bts = PsyRATAccuracyOracle.expandDraw(bts,ndraw,'bt');
            txps = PsyRATAccuracyOracle.expandDraw(txps,ndraw,'txp');
            oxps = PsyRATAccuracyOracle.expandDraw(oxps,ndraw,'oxp');
            txos = PsyRATAccuracyOracle.expandDraw(txos,ndraw,'txo');
            errs = PsyRATAccuracyOracle.expandDraw(errs,ndraw,'err');
            muv = PsyRATAccuracyOracle.expandDraw(muv,ndraw,'mu');

            [ll,pt,ul] = PsyRATAccuracyOracle.rangeStats(@(n) ...
                PsyRATAccuracyOracle.cutscoreTrtDraws(gcoeff,reltype,bps,bos,bts,txps,oxps,txos,errs,muv,cut,n,nocc), ...
                obs, ci);
        end

        function out = diffrel(bp,bt,er_var,wp_cov,obs,est,ci)
            bp1 = bp(:,1,1);
            bp2 = bp(:,2,2);
            bp_cov = bp(:,2,1);

            bt1 = bt(:,1,1);
            bt2 = bt(:,2,2);
            bt_cov = bt(:,2,1);

            wp1 = exp(er_var(:,1)).^2;
            wp2 = exp(er_var(:,2)).^2;

            if isempty(wp_cov)
                wp_cov = zeros(size(wp1));
            else
                wp_cov = wp_cov(:);
                wp_cov = PsyRATAccuracyOracle.expandDraw(wp_cov,length(wp1),'wp_cov');
            end

            obs1 = obs(1);
            obs2 = obs(2);

            uni = bp1 + bp2 - (2 .* bp_cov);
            rel_err = (wp1 ./ obs1) + (wp2 ./ obs2) - ((2 .* wp_cov) ./ psyrat_harmmean(obs));
            abs_err = rel_err + (bt1 ./ obs1) + (bt2 ./ obs2) - ((2 .* bt_cov) ./ psyrat_harmmean(obs));

            denom_rel_err = wp1 + wp2 - (2 .* wp_cov);
            denom_abs_err = denom_rel_err + bt1 + bt2 - (2 .* bt_cov);

            if strcmpi(est,'dep')
                rel = uni ./ (uni + abs_err);
                icc = uni ./ (uni + denom_abs_err);
                out.est_type = 'dependability';
            else
                rel = uni ./ (uni + rel_err);
                icc = uni ./ (uni + denom_rel_err);
                out.est_type = 'generalizability';
            end

            [out.ll,out.pt,out.ul] = PsyRATAccuracyOracle.interval(rel,ci);
            [out.icc_ll,out.icc_pt,out.icc_ul] = PsyRATAccuracyOracle.interval(icc,ci);
            out.bp_cov_pt = mean(bp_cov);
            out.bp_cov_ll = quantile(bp_cov,(1-ci)/2);
            out.bp_cov_ul = quantile(bp_cov,1-((1-ci)/2));
            out.bt_cov_pt = mean(bt_cov);
            out.bt_cov_ll = quantile(bt_cov,(1-ci)/2);
            out.bt_cov_ul = quantile(bt_cov,1-((1-ci)/2));
            out.wp_cov = mean(wp_cov);
            out.wp_cov_ll = quantile(wp_cov,(1-ci)/2);
            out.wp_cov_ul = quantile(wp_cov,1-((1-ci)/2));
        end

        function out = diffrelTrt(bp,bpi,bpo,bt,bo,boi,er_var,er_cov,obs,nocc,reltype,est,ci)
            [bp1,bp2,bp_cov] = PsyRATAccuracyOracle.parseCov(bp,'bp');
            [bpi1,bpi2,bpi_cov] = PsyRATAccuracyOracle.parseCov(bpi,'bpi');
            [bpo1,bpo2,bpo_cov] = PsyRATAccuracyOracle.parseCov(bpo,'bpo');
            [bt1,bt2,bt_cov] = PsyRATAccuracyOracle.parseCov(bt,'bt');
            [bo1,bo2,bo_cov] = PsyRATAccuracyOracle.parseCov(bo,'bo');
            [boi1,boi2,boi_cov] = PsyRATAccuracyOracle.parseCov(boi,'boi');

            wp1 = exp(er_var(:,1)).^2;
            wp2 = exp(er_var(:,2)).^2;

            ndraw = max([length(bp1),length(bpi1),length(bpo1),length(bt1),...
                length(bo1),length(boi1),length(wp1),length(er_cov)]);

            bp1 = PsyRATAccuracyOracle.expandDraw(bp1,ndraw,'bp1');
            bp2 = PsyRATAccuracyOracle.expandDraw(bp2,ndraw,'bp2');
            bp_cov = PsyRATAccuracyOracle.expandDraw(bp_cov,ndraw,'bp_cov');

            bpi1 = PsyRATAccuracyOracle.expandDraw(bpi1,ndraw,'bpi1');
            bpi2 = PsyRATAccuracyOracle.expandDraw(bpi2,ndraw,'bpi2');
            bpi_cov = PsyRATAccuracyOracle.expandDraw(bpi_cov,ndraw,'bpi_cov');

            bpo1 = PsyRATAccuracyOracle.expandDraw(bpo1,ndraw,'bpo1');
            bpo2 = PsyRATAccuracyOracle.expandDraw(bpo2,ndraw,'bpo2');
            bpo_cov = PsyRATAccuracyOracle.expandDraw(bpo_cov,ndraw,'bpo_cov');

            bt1 = PsyRATAccuracyOracle.expandDraw(bt1,ndraw,'bt1');
            bt2 = PsyRATAccuracyOracle.expandDraw(bt2,ndraw,'bt2');
            bt_cov = PsyRATAccuracyOracle.expandDraw(bt_cov,ndraw,'bt_cov');

            bo1 = PsyRATAccuracyOracle.expandDraw(bo1,ndraw,'bo1');
            bo2 = PsyRATAccuracyOracle.expandDraw(bo2,ndraw,'bo2');
            bo_cov = PsyRATAccuracyOracle.expandDraw(bo_cov,ndraw,'bo_cov');

            boi1 = PsyRATAccuracyOracle.expandDraw(boi1,ndraw,'boi1');
            boi2 = PsyRATAccuracyOracle.expandDraw(boi2,ndraw,'boi2');
            boi_cov = PsyRATAccuracyOracle.expandDraw(boi_cov,ndraw,'boi_cov');

            wp1 = PsyRATAccuracyOracle.expandDraw(wp1,ndraw,'wp1');
            wp2 = PsyRATAccuracyOracle.expandDraw(wp2,ndraw,'wp2');
            er_cov = PsyRATAccuracyOracle.expandDraw(er_cov(:),ndraw,'er_cov');

            obs1 = obs(1);
            obs2 = obs(2);
            nocc1 = nocc(1);
            nocc2 = nocc(2);
            hmi = psyrat_harmmean(obs);
            hmo = psyrat_harmmean(nocc);

            bp_raw = bp1 + bp2 - (2 .* bp_cov);
            bpi_raw = bpi1 + bpi2 - (2 .* bpi_cov);
            bpo_raw = bpo1 + bpo2 - (2 .* bpo_cov);
            bt_raw = bt1 + bt2 - (2 .* bt_cov);
            bo_raw = bo1 + bo2 - (2 .* bo_cov);
            boi_raw = boi1 + boi2 - (2 .* boi_cov);
            bpoi_raw = wp1 + wp2 - (2 .* er_cov);

            bpiAdj = (bpi1 ./ obs1) + (bpi2 ./ obs2) - ((2 .* bpi_cov) ./ hmi);
            bpoAdj = (bpo1 ./ nocc1) + (bpo2 ./ nocc2) - ((2 .* bpo_cov) ./ hmo);
            btAdj = (bt1 ./ obs1) + (bt2 ./ obs2) - ((2 .* bt_cov) ./ hmi);
            boAdj = (bo1 ./ nocc1) + (bo2 ./ nocc2) - ((2 .* bo_cov) ./ hmo);
            boiAdj = (boi1 ./ (obs1*nocc1)) + (boi2 ./ (obs2*nocc2)) - ((2 .* boi_cov) ./ (hmi*hmo));
            bpoiAdj = (wp1 ./ (obs1*nocc1)) + (wp2 ./ (obs2*nocc2)) - ((2 .* er_cov) ./ (hmi*hmo));

            if reltype == 1
                uni = bp_raw + bpoAdj;
                rel_err = bpiAdj + bpoiAdj;
                abs_err = rel_err + btAdj + boiAdj;

                uni_icc = bp_raw + bpo_raw;
                rel_err_icc = bpi_raw + bpoi_raw;
                abs_err_icc = rel_err_icc + bt_raw + boi_raw;
            elseif reltype == 2
                uni = bp_raw + bpiAdj;
                rel_err = bpoAdj + bpoiAdj;
                abs_err = rel_err + boAdj + boiAdj;

                uni_icc = bp_raw + bpi_raw;
                rel_err_icc = bpo_raw + bpoi_raw;
                abs_err_icc = rel_err_icc + bo_raw + boi_raw;
            else
                uni = bp_raw;
                rel_err = bpiAdj + bpoAdj + bpoiAdj;
                abs_err = rel_err + btAdj + boAdj + boiAdj;

                uni_icc = bp_raw;
                rel_err_icc = bpi_raw + bpo_raw + bpoi_raw;
                abs_err_icc = rel_err_icc + bt_raw + bo_raw + boi_raw;
            end

            if strcmpi(est,'dep')
                rel = uni ./ (uni + abs_err);
                icc = uni_icc ./ (uni_icc + abs_err_icc);
                out.est_type = 'dependability';
            else
                rel = uni ./ (uni + rel_err);
                icc = uni_icc ./ (uni_icc + rel_err_icc);
                out.est_type = 'generalizability';
            end

            [out.ll,out.pt,out.ul] = PsyRATAccuracyOracle.interval(rel,ci);
            [out.icc_ll,out.icc_pt,out.icc_ul] = PsyRATAccuracyOracle.interval(icc,ci);
            out.wp_cov = mean(er_cov);
            out.bp_cov_pt = mean(bp_cov);
            out.bpi_cov_pt = mean(bpi_cov);
            out.bpo_cov_pt = mean(bpo_cov);
            out.bt_cov_pt = mean(bt_cov);
            out.bo_cov_pt = mean(bo_cov);
            out.boi_cov_pt = mean(boi_cov);
        end

        function out = ssrel(bp,wp_pop,wp_ss,idtable,ci,i,gcoeff)
            % Subject-level reliability with the relative/absolute split. i is
            % the trial main-effect SD draws (sigma_i); gcoeff selects the
            % decision type (1 = absolute/dependability phi_s, keeps sigma_i^2;
            % 2 = relative/generalizability G_s, excludes it). Defaults
            % reproduce the pre-trial-effect absolute-only behavior.
            if nargin < 6 || isempty(i)
                i = 0;
            end
            if nargin < 7 || isempty(gcoeff)
                gcoeff = 1;
            end

            bpVar = cell2mat(bp).^2;
            % After the trial main effect is modeled on the mean, this is the
            % per-subject RELATIVE residual sigma_pi,e^2(s).
            wpVar = exp(wp_pop + cell2mat(wp_ss)).^2;
            i2 = i(:).^2;
            i2 = PsyRATAccuracyOracle.expandDraw(i2,size(wpVar,1),'i');
            ciedge = (1-ci)/2;

            if gcoeff == 1
                errwp = wpVar + i2; % absolute keeps sigma_i^2
            else
                errwp = wpVar;      % relative excludes sigma_i^2
            end

            dep = bpVar ./ (bpVar + (errwp ./ idtable.trls(:)'));
            icc = bpVar ./ (bpVar + errwp);
            sem = sqrt(errwp ./ idtable.trls(:)');

            out.dep_ll = quantile(dep,ciedge)';
            out.dep_pt = mean(dep)';
            out.dep_ul = quantile(dep,1-ciedge)';

            out.icc_ll = quantile(icc,ciedge)';
            out.icc_pt = mean(icc)';
            out.icc_ul = quantile(icc,1-ciedge)';

            out.sem_ll = quantile(sem,ciedge)';
            out.sem_pt = mean(sem)';
            out.sem_ul = quantile(sem,1-ciedge)';

            out.bp_var = repmat(mean(bpVar), height(idtable), 1);
            out.ss_errvar = mean(wpVar,1)';
            out.trl_var = repmat(mean(i2), height(idtable), 1);
            out.pop_errvar = repmat(mean(exp(wp_pop).^2), height(idtable), 1);
        end

        function out = sscutscore(bp,wp_pop,wp_ss,mu,cut,idtable,ci,est,i)
            if nargin < 9 || isempty(i)
                i = 0;
            end
            if nargin < 8 || isempty(est)
                est = 'dep';
            end
            % Cut-scores are absolute-error only (Rocha 2026, Table 2; audit
            % F7). A relative-error (est = 'gen') cut-score is not a valid
            % quantity.
            if ~strcmpi(est,'dep')
                error('oracle:est', ...
                    ['Subject-level cut-scores are absolute-error only ', ...
                    '(est = ''dep''). A relative-error cut-score is not a ', ...
                    'valid quantity. How to fix: request est = ''dep''.']);
            end

            bpVar = cell2mat(bp).^2;
            % Per-subject RELATIVE residual sigma_pi,e^2(s).
            wpVar = exp(wp_pop + cell2mat(wp_ss)).^2;
            mu = PsyRATAccuracyOracle.coerceMu(mu,size(wpVar,1),size(wpVar,2));
            i = i(:).^2;
            i = PsyRATAccuracyOracle.expandDraw(i,size(wpVar,1),'i');
            ciedge = (1-ci)/2;

            offsq = (mu - cut).^2;
            % Absolute error keeps the trial main effect sigma_i^2.
            err = (wpVar ./ idtable.trls(:)') + (i ./ idtable.trls(:)');

            rel = (bpVar + offsq) ./ (bpVar + offsq + err);
            out.ll = quantile(rel,ciedge)';
            out.pt = mean(rel)';
            out.ul = quantile(rel,1-ciedge)';
        end

        function out = ssrelTrt(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,idtable,ci,nocc,txp_ss,oxp_ss,err_ss)
            if nargin < 12 || isempty(nocc)
                nocc = 1;
            end
            if nargin < 13
                txp_ss = [];
            end
            if nargin < 14
                oxp_ss = [];
            end
            if nargin < 15
                err_ss = [];
            end

            ciedge = (1-ci)/2;
            out.rel_pt = zeros(height(idtable),1);
            out.rel_ll = zeros(height(idtable),1);
            out.rel_ul = zeros(height(idtable),1);
            out.icc_pt = zeros(height(idtable),1);
            out.icc_ll = zeros(height(idtable),1);
            out.icc_ul = zeros(height(idtable),1);
            out.sem_pt = zeros(height(idtable),1);
            out.sem_ll = zeros(height(idtable),1);
            out.sem_ul = zeros(height(idtable),1);

            btv = bt.^2;
            bov = bo.^2;
            txov = txo.^2;
            bpv = bp.^2;

            for r = 1:height(idtable)
                txp_draw = txp;
                oxp_draw = oxp;
                err_draw = err;

                has_ss = ~isempty(txp_ss) && ~isempty(oxp_ss) && ~isempty(err_ss) && ...
                    any(strcmp(idtable.Properties.VariableNames,'id2'));
                if has_ss
                    sid = round(idtable.id2(r));
                    txp_draw = txp_ss(:,sid);
                    oxp_draw = oxp_ss(:,sid);
                    err_draw = err_ss(:,sid);
                end

                obs = max(1,idtable.trls(r));
                [rll,rpt,rul] = PsyRATAccuracyOracle.relTrt(gcoeff,reltype,bp,bo,bt,txp_draw,oxp_draw,txo,err_draw,obs,nocc,ci);
                %The subject-level ICC is the same two-facet partition at
                %obs = nocc = 1, NOT the one-facet PsyRATAccuracyOracle.icc.
                %Mirrors psyrat_ssrel_trt; see the pinning note there.
                [ill,ipt,iul] = PsyRATAccuracyOracle.relTrt(gcoeff,reltype,bp,bo,bt,txp_draw,oxp_draw,txo,err_draw,1,1,ci);

                txpv = txp_draw.^2;
                oxpv = oxp_draw.^2;
                errv = err_draw.^2;

                if gcoeff == 1 && reltype == 1
                    err_term = (txpv ./ obs) + (errv ./ (obs*nocc)) + (btv ./ obs) + (txov ./ (obs*nocc));
                elseif gcoeff == 1 && reltype == 2
                    err_term = (oxpv ./ nocc) + (errv ./ (obs*nocc)) + (bov ./ nocc) + (txov ./ (obs*nocc));
                elseif gcoeff == 1 && reltype == 3
                    err_term = (txpv ./ obs) + (oxpv ./ nocc) + (errv ./ (obs*nocc)) + (btv ./ obs) + (bov ./ nocc) + (txov ./ (obs*nocc));
                elseif gcoeff == 2 && reltype == 1
                    err_term = (txpv ./ obs) + (errv ./ (obs*nocc));
                elseif gcoeff == 2 && reltype == 2
                    err_term = (oxpv ./ nocc) + (errv ./ (obs*nocc));
                else
                    err_term = (txpv ./ obs) + (oxpv ./ nocc) + (errv ./ (obs*nocc));
                end

                sem_draw = sqrt(err_term);
                out.rel_ll(r) = rll;
                out.rel_pt(r) = rpt;
                out.rel_ul(r) = rul;
                out.icc_ll(r) = ill;
                out.icc_pt(r) = ipt;
                out.icc_ul(r) = iul;
                out.sem_ll(r) = quantile(sem_draw,ciedge);
                out.sem_pt(r) = mean(sem_draw);
                out.sem_ul(r) = quantile(sem_draw,1-ciedge);
            end

            out.dep_ll = out.rel_ll;
            out.dep_pt = out.rel_pt;
            out.dep_ul = out.rel_ul;
            out.bp_var = repmat(mean(bpv), height(idtable), 1);

            %Group reference lines for the caterpillar plots, from the POPULATION
            %draws. Mirrors psyrat_ssrel_trt's gro_icc / gro_rel columns.
            gro_obs = max(1,mean(idtable.trls));
            [~,gro_icc,~] = PsyRATAccuracyOracle.relTrt( ...
                gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,1,1,ci);
            [~,gro_rel,~] = PsyRATAccuracyOracle.relTrt( ...
                gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,gro_obs,nocc,ci);
            out.gro_icc = repmat(gro_icc, height(idtable), 1);
            out.gro_rel = repmat(gro_rel, height(idtable), 1);
        end

        function out = ssrelDiff(bp,bt,er_var,idtable,ci,est,wp_cov,er_var_ss1,er_var_ss2,wp_cov_ss)
            if nargin < 6 || isempty(est)
                est = 'dep';
            end
            if nargin < 7 || isempty(wp_cov)
                wp_cov = 0;
            end
            if nargin < 8
                er_var_ss1 = [];
            end
            if nargin < 9
                er_var_ss2 = [];
            end
            if nargin < 10
                wp_cov_ss = [];
            end

            ciedge = (1-ci)/2;
            out.rel_pt = zeros(height(idtable),1);
            out.rel_ll = zeros(height(idtable),1);
            out.rel_ul = zeros(height(idtable),1);
            out.icc_pt = zeros(height(idtable),1);
            out.icc_ll = zeros(height(idtable),1);
            out.icc_ul = zeros(height(idtable),1);
            out.sem_pt = zeros(height(idtable),1);
            out.sem_ll = zeros(height(idtable),1);
            out.sem_ul = zeros(height(idtable),1);

            bt1 = bt(:,1,1);
            bt2 = bt(:,2,2);
            bt_cov = bt(:,2,1);

            for r = 1:height(idtable)
                obs = [max(1,idtable.trls1(r)) max(1,idtable.trls2(r))];
                er_draw = er_var;
                wp_cov_draw = wp_cov;
                wp1 = exp(er_draw(:,1)).^2;
                wp2 = exp(er_draw(:,2)).^2;

                has_ss = ~isempty(er_var_ss1) && ~isempty(er_var_ss2) && ...
                    ~isempty(wp_cov_ss) && any(strcmp(idtable.Properties.VariableNames,'id2'));
                if has_ss
                    sid = round(idtable.id2(r));
                    er_draw = [er_var_ss1(:,sid) er_var_ss2(:,sid)];
                    wp_cov_draw = wp_cov_ss(:,sid);
                    wp1 = exp(er_draw(:,1)).^2;
                    wp2 = exp(er_draw(:,2)).^2;
                else
                    wp_cov_draw = PsyRATAccuracyOracle.expandDraw(wp_cov_draw(:),size(er_var,1),'wp_cov');
                end

                diffscore = PsyRATAccuracyOracle.diffrel(bp,bt,er_draw,wp_cov_draw,obs,est,ci);

                rel_err = (wp1 ./ obs(1)) + (wp2 ./ obs(2)) - ((2 .* wp_cov_draw) ./ psyrat_harmmean(obs));
                abs_err = rel_err + (bt1 ./ obs(1)) + (bt2 ./ obs(2)) - ((2 .* bt_cov) ./ psyrat_harmmean(obs));
                %Mirrors psyrat_ssrel_diff: a negative error variance has no
                %real square root, so the SEM is undefined for that draw and is
                %reported as NaN rather than floored at zero. This duplicated
                %the old max(...,0) and had to move with production, or the
                %suite would have stayed green over a stale expectation.
                sem_draw = nan(size(abs_err));
                sem_draw(abs_err >= 0) = sqrt(abs_err(abs_err >= 0));

                out.rel_pt(r) = diffscore.pt;
                out.rel_ll(r) = diffscore.ll;
                out.rel_ul(r) = diffscore.ul;
                out.icc_pt(r) = diffscore.icc_pt;
                out.icc_ll(r) = diffscore.icc_ll;
                out.icc_ul(r) = diffscore.icc_ul;
                out.sem_pt(r) = mean(sem_draw);
                out.sem_ll(r) = quantile(sem_draw,ciedge);
                out.sem_ul(r) = quantile(sem_draw,1-ciedge);
            end

            out.dep_pt = out.rel_pt;
            out.dep_ll = out.rel_ll;
            out.dep_ul = out.rel_ul;
        end

        function [ll,pt,ul] = splitsSingle(gcoeff,sig_p,sig_s,sig_ps,sig_e,nsplit,nbar_i,ci)
            % Rocha Table 4 (one-facet data splits, observed design). From split
            % means the per-item residual sig_e^2 lumps the item main effect with
            % the person-involving item residual: dependability (gcoeff 1) is
            % exact, generalizability (gcoeff 2) excludes the split main effect
            % sig_s^2 and is a downward-biased lower bound. resid is the per-item
            % residual scaled to one split (variance / nbar_i).
            sp = sig_p(:).^2; ss = sig_s(:).^2; sps = sig_ps(:).^2;
            resid = sig_e(:).^2 ./ nbar_i;
            if gcoeff == 1
                errterm = (sps + ss + resid) ./ nsplit;
            else
                errterm = (sps + resid) ./ nsplit;
            end
            draws = sp ./ (sp + errterm);
            [ll,pt,ul] = PsyRATAccuracyOracle.interval(draws,ci);
        end

        function [ll,pt,ul] = splitsSingleSem(kind,sig_s,sig_ps,sig_e,nsplit,nbar_i,ci)
            % SEM = sqrt(error variance) for Table 4. kind = 'rel' (lower-bound
            % relative, an upper bound on the true relative SEM) or 'abs' (exact).
            ss = sig_s(:).^2; sps = sig_ps(:).^2; resid = sig_e(:).^2 ./ nbar_i;
            relvar = (sps + resid) ./ nsplit;
            absvar = (sps + ss + resid) ./ nsplit;
            if strcmpi(kind,'rel'); v = relvar; else; v = absvar; end
            [ll,pt,ul] = PsyRATAccuracyOracle.interval(sqrt(v),ci);
        end

        function [ll,pt,ul] = splitsTrt(gcoeff,reltype,sig_p,sig_o,sig_s,sig_ps,...
                sig_po,sig_os,sig_pos,sig_e,nsplit,nocc,nbar_i,ci)
            % Rocha Table 5 (two-facet data splits, observed design, nested-items
            % component set). The split facet maps onto the test-retest "trial"
            % slot; the effective highest-order residual folds the p:o:s 3-way and
            % the n_i-scaled per-item residual, so it lands at /(n_s n_o) and the
            % item part at /(n_s n_o n_i). Reuses relTrtDraws (inputs are
            % variances).
            bp = sig_p(:).^2; bo = sig_o(:).^2; bt = sig_s(:).^2;
            txp = sig_ps(:).^2; oxp = sig_po(:).^2; txo = sig_os(:).^2;
            erreff = sig_pos(:).^2 + sig_e(:).^2 ./ nbar_i;
            draws = PsyRATAccuracyOracle.relTrtDraws(gcoeff,reltype,bp,bo,bt,...
                txp,oxp,txo,erreff,nsplit,nocc);
            [ll,pt,ul] = PsyRATAccuracyOracle.interval(draws,ci);
        end

        function [ll,pt,ul] = splitsTrtSem(kind,reltype,sig_o,sig_s,sig_ps,...
                sig_po,sig_os,sig_pos,sig_e,nsplit,nocc,nbar_i,ci)
            % SEM = sqrt(error variance) for Table 5, per reltype (matching
            % relTrtDraws's err_term). kind = 'rel' or 'abs'.
            bo = sig_o(:).^2; bt = sig_s(:).^2; txp = sig_ps(:).^2;
            oxp = sig_po(:).^2; txo = sig_os(:).^2;
            er = sig_pos(:).^2 + sig_e(:).^2 ./ nbar_i;
            switch reltype
                case 1
                    relv = txp ./ nsplit + er ./ (nsplit*nocc);
                    absv = relv + bt ./ nsplit + txo ./ (nsplit*nocc);
                case 2
                    relv = oxp ./ nocc + er ./ (nsplit*nocc);
                    absv = relv + bo ./ nocc + txo ./ (nsplit*nocc);
                case 3
                    relv = txp ./ nsplit + oxp ./ nocc + er ./ (nsplit*nocc);
                    absv = relv + bt ./ nsplit + bo ./ nocc + txo ./ (nsplit*nocc);
            end
            if strcmpi(kind,'rel'); v = relv; else; v = absv; end
            [ll,pt,ul] = PsyRATAccuracyOracle.interval(sqrt(v),ci);
        end

        function comp = dodComposite(varcov, cvec, rule, n, tol)
            % Regression mirror of PSYRAT_DOD_COMPOSITE (analyses 28/29):
            % the signed four-cell quadratic form, per-cell count divisors on
            % the variance terms and pairwise harmonic-mean divisors on the
            % covariance terms under 'harmonic', the reference's magnitude
            % scale, the tiny-negative floor, and the material-negative flag.
            if nargin < 5 || isempty(tol), tol = 1e-10; end
            if strcmpi(rule,'unscaled'), n = ones(1,4); end
            cvec = cvec(:); n = n(:)';
            nd = size(varcov,1);
            vars = zeros(nd,4);
            for q = 1:4
                vars(:,q) = varcov(:,q,q);
            end
            value = (vars ./ n) * cvec.^2;
            abscov = zeros(nd,1);
            for a = 1:3
                for b = (a+1):4
                    nharm = 2 / ((1/n(a)) + (1/n(b)));
                    value = value + 2 .* cvec(a) .* cvec(b) .* varcov(:,a,b) ./ nharm;
                    abscov = abscov + abs(varcov(:,a,b));
                end
            end
            scale = sum(abs(vars ./ n),2) + 2 .* abscov;
            threshold = tol .* max(1, abs(scale));
            tiny = value < 0 & value >= -threshold;
            value(tiny) = 0;
            comp = struct('value', value, ...
                'invalid_negative', double(value < -threshold), ...
                'scale', scale);
        end

        function mats = dodComponentsLs(facet, alpha, log_sigma, sd_p, ...
                cor_p, sd_trl, cor_trl, extra)
            % Regression mirror of PSYRAT_GAMMA_DOD_COMPONENTS_LS: the
            % observed-scale 4x4 component matrices assembled from the
            % location-scale converter formulas, with the six residual
            % cross-covariances at exactly zero (the nonconcurrent estimand).
            % One-facet crosses use the product-identity grouping
            % (PSYRAT_GAMMA_CROSSCOV_RESCOR); two-facet diagonals and crosses
            % use the subtraction groupings of PSYRAT_GAMMA_VARCOMPS_TRT_LS /
            % PSYRAT_GAMMA_CROSSCOV_TRT, including their max(.,0) clips.
            two = strcmpi(facet,'trial_occasion');
            nd = size(alpha,1);
            z44 = zeros(nd,4,4);
            if two
                keys = {'p','i','o','pi','po','io','pio_e'};
            else
                keys = {'p','i','pi_e'};
            end
            for k = 1:numel(keys), mats.(keys{k}) = z44; end
            if two
                vtot = sd_p.^2 + extra.sd_occ.^2 + sd_trl.^2 + ...
                    extra.sd_oid.^2 + extra.sd_tid.^2 + extra.sd_to.^2;
            else
                vtot = sd_p.^2 + sd_trl.^2;
            end
            for q = 1:4
                if two
                    vp = sd_p(:,q).^2; vo = extra.sd_occ(:,q).^2;
                    vt = sd_trl(:,q).^2; vpo = extra.sd_oid(:,q).^2;
                    vpt = extra.sd_tid(:,q).^2; vot = extra.sd_to(:,q).^2;
                    B = exp(2.*alpha(:,q) + vtot(:,q));
                    p = B .* (exp(vp) - 1);
                    o = B .* (exp(vo) - 1);
                    t = B .* (exp(vt) - 1);
                    po = max(B .* (exp(vp+vo+vpo) - exp(vp) - exp(vo) + 1), 0);
                    pt = max(B .* (exp(vp+vt+vpt) - exp(vp) - exp(vt) + 1), 0);
                    ot = max(B .* (exp(vo+vt+vot) - exp(vo) - exp(vt) + 1), 0);
                    var_mu = B .* (exp(vtot(:,q)) - 1);
                    res = max(var_mu + exp(2.*log_sigma(:,q)) ...
                        - p - o - t - po - pt - ot, 0);
                    mats.p(:,q,q) = p; mats.i(:,q,q) = t; mats.o(:,q,q) = o;
                    mats.pi(:,q,q) = pt; mats.po(:,q,q) = po;
                    mats.io(:,q,q) = ot; mats.pio_e(:,q,q) = res;
                else
                    sp2 = sd_p(:,q).^2; si2 = sd_trl(:,q).^2;
                    p = exp(2.*alpha(:,q) + si2) .* exp(sp2) .* (exp(sp2) - 1);
                    i = exp(2.*alpha(:,q) + sp2) .* exp(si2) .* (exp(si2) - 1);
                    var_mu = exp(2.*alpha(:,q)) .* ...
                        (exp(2.*(sp2+si2)) - exp(sp2+si2));
                    pi_e = max(var_mu + exp(2.*log_sigma(:,q)) - p - i, 0);
                    mats.p(:,q,q) = p; mats.i(:,q,q) = i;
                    mats.pi_e(:,q,q) = pi_e;
                end
            end
            for a = 1:3
                for b = (a+1):4
                    pref = exp(alpha(:,a) + alpha(:,b) + ...
                        0.5 .* (vtot(:,a) + vtot(:,b)));
                    cp = cor_p(:,a,b) .* sd_p(:,a) .* sd_p(:,b);
                    ci = cor_trl(:,a,b) .* sd_trl(:,a) .* sd_trl(:,b);
                    if two
                        co = extra.cor_occ(:,a,b) .* extra.sd_occ(:,a) .* extra.sd_occ(:,b);
                        cpi = extra.cor_tid(:,a,b) .* extra.sd_tid(:,a) .* extra.sd_tid(:,b);
                        cpo = extra.cor_oid(:,a,b) .* extra.sd_oid(:,a) .* extra.sd_oid(:,b);
                        cio = extra.cor_to(:,a,b) .* extra.sd_to(:,a) .* extra.sd_to(:,b);
                        cfun = @(c) pref .* (exp(c) - 1);
                        p = cfun(cp); o = cfun(co); t = cfun(ci);
                        po = cfun(cp + co + cpo) - p - o;
                        pt = cfun(cp + ci + cpi) - p - t;
                        ot = cfun(co + ci + cio) - o - t;
                        total = cfun(cp + co + ci + cpo + cpi + cio);
                        res = total - p - o - t - po - pt - ot;
                        vals = {p, t, o, pt, po, ot, res};
                    else
                        p = pref .* expm1(cp);
                        i = pref .* expm1(ci);
                        pim = pref .* expm1(cp) .* expm1(ci);
                        vals = {p, i, pim};
                    end
                    for k = 1:numel(keys)
                        mats.(keys{k})(:,a,b) = vals{k};
                        mats.(keys{k})(:,b,a) = vals{k};
                    end
                end
            end
        end

        function mats = dodComponentsGauss(facet, alpha, log_sigma, sd_p, ...
                cor_p, sd_trl, cor_trl, extra)
            % Regression mirror of PSYRAT_DOD_COMPONENTS (the GAUSSIAN arm of
            % analyses 28/29): identity link, so every signal component is its
            % model covariance block directly and the residual channel is
            % exactly diagonal (SCIENTIFIC_FORMULA_AUDIT.md section 26).
            % Assembled here PER DRAW via D(sd)*R*D(sd) matrix products - a
            % different evaluation route from production's vectorized
            % entrywise assembly, so agreement is not a byte-copy. alpha is
            % unused in every second moment by the identity-link theorem; it
            % is accepted to keep the mirror's signature parallel to
            % dodComponentsLs.
            two = strcmpi(facet,'trial_occasion');
            nd = size(sd_p,1);
            if two
                blocks = {'p',sd_p,cor_p; 'i',sd_trl,cor_trl; ...
                    'o',extra.sd_occ,extra.cor_occ; ...
                    'pi',extra.sd_tid,extra.cor_tid; ...
                    'po',extra.sd_oid,extra.cor_oid; ...
                    'io',extra.sd_to,extra.cor_to};
                reskey = 'pio_e';
            else
                blocks = {'p',sd_p,cor_p; 'i',sd_trl,cor_trl};
                reskey = 'pi_e';
            end
            for k = 1:size(blocks,1)
                key = blocks{k,1}; sd = blocks{k,2}; cor = blocks{k,3};
                M = zeros(nd,4,4);
                for d = 1:nd
                    D = diag(sd(d,:));
                    M(d,:,:) = D * squeeze(cor(d,:,:)) * D;
                end
                mats.(key) = M;
            end
            R = zeros(nd,4,4);
            for q = 1:4
                R(:,q,q) = exp(2 .* log_sigma(:,q));
            end
            mats.(reskey) = R;
        end
    end

    methods (Static, Access = private)
        function draws = relTrtDraws(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,obs,nocc)
            if gcoeff == 1 && reltype == 1
                uni = bp + (oxp ./ nocc);
                err_term = (txp ./ obs) + (err ./ (obs*nocc)) + (bt ./ obs) + (txo ./ (obs*nocc));
            elseif gcoeff == 1 && reltype == 2
                uni = bp + (txp ./ obs);
                err_term = (oxp ./ nocc) + (err ./ (obs*nocc)) + (bo ./ nocc) + (txo ./ (obs*nocc));
            elseif gcoeff == 1 && reltype == 3
                uni = bp;
                err_term = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc)) + (bt ./ obs) + (bo ./ nocc) + (txo ./ (obs*nocc));
            elseif gcoeff == 2 && reltype == 1
                uni = bp + (oxp ./ nocc);
                err_term = (txp ./ obs) + (err ./ (obs*nocc));
            elseif gcoeff == 2 && reltype == 2
                uni = bp + (txp ./ obs);
                err_term = (oxp ./ nocc) + (err ./ (obs*nocc));
            else
                uni = bp;
                err_term = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc));
            end
            draws = uni ./ (uni + err_term);
        end

        function draws = cutscoreDepDraws(bp,wp,i,mu,cut,obs,est) %#ok<INUSD>
            % Cut-scores are absolute-error only (Rocha 2026, Table 2): the
            % denominator error keeps the trial main effect sigma_i^2. The est
            % argument is retained for signature compatibility but is ignored;
            % a relative-error cut-score is not a valid quantity.
            offsq = (mu - cut).^2;
            err = (wp ./ obs) + (i ./ obs);
            draws = (bp + offsq) ./ (bp + offsq + err);
        end

        function draws = cutscoreTrtDraws(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,mu,cut,obs,nocc)
            % Cut-scores are absolute-error only (Rocha 2026, Table 3): the
            % only valid decision type is gcoeff = 1 (dependability). A
            % relative-error (gcoeff = 2) cut-score is not a valid quantity.
            if gcoeff ~= 1
                error('oracle:cutscoreRelative', ...
                    ['Cut-scores are absolute-error only (gcoeff = 1). ', ...
                    'How to fix: request the absolute cut-score; a relative ', ...
                    'cut-score is not a valid quantity.']);
            end
            offsq = (mu - cut).^2;
            if reltype == 1
                uni = bp + (oxp ./ nocc);
                err_term = (txp ./ obs) + (err ./ (obs*nocc)) + (bt ./ obs) + (txo ./ (obs*nocc));
            elseif reltype == 2
                uni = bp + (txp ./ obs);
                err_term = (oxp ./ nocc) + (err ./ (obs*nocc)) + (bo ./ nocc) + (txo ./ (obs*nocc));
            else % reltype == 3
                uni = bp;
                err_term = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc)) + (bt ./ obs) + (bo ./ nocc) + (txo ./ (obs*nocc));
            end
            draws = (uni + offsq) ./ (uni + offsq + err_term);
        end

        function [ll,pt,ul] = rangeStats(drawFn,obs,ci)
            if length(obs) == 1
                draws = drawFn(obs);
                [ll,pt,ul] = PsyRATAccuracyOracle.interval(draws,ci);
            else
                % Mirror the production range layout (psyrat_dep /
                % psyrat_rel_sing): one column per trial count, aligned to
                % obs(1):obs(2) with no leading zeros when obs(1) > 1.
                nobs = obs(2)-obs(1)+1;
                ll = zeros(1,nobs);
                pt = zeros(1,nobs);
                ul = zeros(1,nobs);
                for i = obs(1):obs(2)
                    draws = drawFn(i);
                    [lli,pti,uli] = PsyRATAccuracyOracle.interval(draws,ci);
                    idx = i - obs(1) + 1;
                    ll(idx) = lli;
                    pt(idx) = pti;
                    ul(idx) = uli;
                end
            end
        end

        function [ll,pt,ul] = interval(draws,ci)
            ciedge = (1-ci)/2;
            ll = quantile(draws,ciedge);
            pt = mean(draws);
            ul = quantile(draws,1-ciedge);
        end

        function out = expandDraw(v,ndraw,name)
            v = v(:);
            if length(v) == 1
                out = repmat(v,ndraw,1);
            elseif length(v) == ndraw
                out = v;
            else
                error('oracle:drawsize', ...
                    ['Draw-size mismatch for %s. Expected a scalar or %d draws. ',...
                    'How to fix: pass oracle inputs with matching posterior draw counts.'], ...
                    name, ndraw);
            end
        end

        function [v1,v2,cv] = parseCov(inp,name)
            if ndims(inp) ~= 3 || size(inp,2) ~= 2 || size(inp,3) ~= 2
                error('oracle:shape', ...
                    ['%s should be [draw x 2 x 2]. How to fix: provide ',...
                    'a covariance matrix with diagonal variances and off-diagonal covariance.'], ...
                    name);
            end
            v1 = inp(:,1,1);
            v2 = inp(:,2,2);
            cv = inp(:,2,1);
        end

        function mu = coerceMu(mu,nrow,ncol)
            if iscell(mu)
                mu = cell2mat(mu);
            end
            if isvector(mu)
                if length(mu) == nrow
                    mu = repmat(mu(:),1,ncol);
                elseif length(mu) == ncol
                    mu = repmat(mu(:)',nrow,1);
                else
                    error('oracle:mu', ...
                        ['mu shape is incompatible with draw/subject dimensions. ',...
                        'How to fix: pass mu as [draw x subject], [subject x draw], ',...
                        'or a vector of draw or subject length.']);
                end
            elseif size(mu,1) == nrow && size(mu,2) == ncol
                % expected shape
            elseif size(mu,2) == nrow && size(mu,1) == ncol
                mu = mu';
            else
                error('oracle:mu', ...
                    ['mu shape is incompatible with draw/subject dimensions. ',...
                    'How to fix: pass mu as [draw x subject], [subject x draw], ',...
                    'or a vector of draw or subject length.']);
            end
        end
    end
end
