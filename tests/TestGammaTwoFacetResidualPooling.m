classdef TestGammaTwoFacetResidualPooling < PsyRATTestBase
    %TESTGAMMATWOFACETRESIDUALPOOLING Q8/PBA-08 (owner ruling 2026-09-13,
    % gamma_target_version 2): every two-facet gamma subject-level residual
    % (analyses 25 and 27, log-nu and location-scale) is
    %
    %   E_p = R + C_p,
    %
    % the POPULATION highest-order (person x occasion x trial) mean-surface
    % remainder plus the participant's own conditional observation variance,
    % with the six population mean components P, O, I, PO, PI, OI retained at
    % the participant's own z and the existing CE/CS/CES allocations, ICC
    % (both counts 1) and draw-level SEM unchanged. No CmdStan.
    %
    % INDEPENDENT REFERENCES. The anchors are the frozen expectations in
    % the 2C Q8/Q9 expectation set (local development evidence, not shipped)
    % (70-digit Decimal, support-union product enumeration, produced without
    % running production code), transcribed here to double precision. The
    % in-test remainder is the memo's INCLUSION-EXCLUSION closed form (eq. 2,
    % localClosedFormRemainder), a different algorithm from the production
    % helper's support-union product enumeration. The in-test coefficients
    % are the contract's own U/(U+A) and U/(U+A+B) arithmetic
    % (localContractMetrics), not psyrat_rel_trt.
    %
    % WHAT IS EXECUTED. The production helper psyrat_gamma_mean_remainder_trt,
    % the production converter psyrat_gamma_varcomps_trt (the six components),
    % the two production dynamic calculators psyrat_ssrel_dynrel_trt_gamma /
    % _ls, and the production static consumer psyrat_ssrel_trt fed exactly as
    % psyrat_relsummary's local_ssrel_trt_call feeds it. The two static
    % extractors are local functions of psyrat_computevarcomp and cannot be
    % called here (the repo's standing limitation); their ENCODING is the same
    % two-line pop_sdlog/ind_sdlog construction exercised through
    % psyrat_ssrel_trt below, and their WIRING (helper call, encoding form,
    % old form absent, v2 stamp) is pinned by source in the last two tests.
    % The live check of the extractors is the CmdStan recovery lane.

    methods (Static, Access = private)

        function R = localClosedFormRemainder(alpha, v)
            %Memo eq. 2, the inclusion-exclusion closed form over shared-index
            %sets: R = B [e^V - e^(vp+vo+vpo) - e^(vp+vt+vpt) - e^(vo+vt+vot)
            %+ e^vp + e^vo + e^vt - 1]. v is a struct with fields P, O, I, PO,
            %PI, OI (log variances). Deliberately NOT the production helper's
            %support-union enumeration; adequate at the fixtures' variance
            %scale (it cancels badly only in the additive limit).
            V = v.P + v.O + v.I + v.PO + v.PI + v.OI;
            B = exp(2*alpha + V);
            R = B * (exp(V) - exp(v.P + v.O + v.PO) - exp(v.P + v.I + v.PI) ...
                - exp(v.O + v.I + v.OI) + exp(v.P) + exp(v.O) + exp(v.I) - 1);
        end

        function m = localContractMetrics(c, E, n_i, n_o, reltype)
            %Contract table: c has fields P,O,I,PO,PI,OI (observed-scale
            %variances); E is the residual variance in the E/(n_i n_o) slot.
            switch reltype
                case 1  % CE
                    U = c.P + c.PO/n_o;  A = c.PI/n_i + E/(n_i*n_o);
                    Bx = c.I/n_i + c.OI/(n_i*n_o);
                    U1 = c.P + c.PO;     A1 = c.PI + E;  B1 = c.I + c.OI;
                case 2  % CS
                    U = c.P + c.PI/n_i;  A = c.PO/n_o + E/(n_i*n_o);
                    Bx = c.O/n_o + c.OI/(n_i*n_o);
                    U1 = c.P + c.PI;     A1 = c.PO + E;  B1 = c.O + c.OI;
                otherwise  % CES
                    U = c.P;  A = c.PI/n_i + c.PO/n_o + E/(n_i*n_o);
                    Bx = c.I/n_i + c.O/n_o + c.OI/(n_i*n_o);
                    U1 = c.P; A1 = c.PI + c.PO + E;  B1 = c.I + c.O + c.OI;
            end
            m.G = U/(U+A);          m.D = U/(U+A+Bx);
            m.ICC_G = U1/(U1+A1);   m.ICC_D = U1/(U1+A1+B1);
            m.SEM_rel = sqrt(A);    m.SEM_abs = sqrt(A+Bx);
        end

        function c = localComponents(vc)
            %Six observed-scale components from the production converter's
            %struct, under the contract's names (varcomps: sigma_po2 = person x
            %occasion, sigma_pt2 = person x trial).
            c = struct('P', vc.sigma_p2, 'O', vc.sigma_o2, 'I', vc.sigma_t2, ...
                'PO', vc.sigma_po2, 'PI', vc.sigma_pt2, 'OI', vc.sigma_ot2);
        end

        function f = localMainFixture()
            %Memo/contract main fixture: alpha = log 5, s_p = .6, s_o = s_i = .5,
            %zero log two-way SDs, n_i = 3, n_o = 2; LS sigma = .2; log-nu nu = 18.
            f.alpha = log(5); f.s_p = .6; f.s_o = .5; f.s_t = .5;
            f.s_po = 0; f.s_pt = 0; f.s_ot = 0;
            f.n_i = 3; f.n_o = 2;
            f.logsig = log(.2); f.lognu = log(18);
            %frozen (Decimal -> double)
            f.R  = 2.0652177285562118;
            f.P  = 25.600675997888492; f.O = 16.779942518247016;
            f.I  = 16.779942518247016; f.OI = 4.7659301657414588;
            f.PI = 7.2712426677881410; f.PO = 7.2712426677881410;
            f.C_ls = .04;                          f.E_ls = 2.1052177285562118;
            f.C_nu = 7.5507828568306812;           f.E_nu = 9.6160005853868930;
            %frozen person-target metrics, reltype 1 (CE), 2 (CS), 3 (CES)
            f.ls(3) = struct('G',.79974834803861693,'D',.54715718947985399, ...
                'ICC_G',.60595640747083997,'ICC_D',.31772798017597068, ...
                'SEM_rel',2.5318448829360550,'SEM_abs',4.6030257046221502);
            f.ls(1) = struct('G',.91332277693711746,'D',.76139065265133028, ...
                'ICC_G',.77806342859735215,'ICC_D',.51528026061898091, ...
                'SEM_rel',1.6657182166727407,'SEM_abs',3.0269213806101182);
            f.ls(2) = struct('G',.87546463397095062,'D',.68028358831044403, ...
                'ICC_G',.77806342859735215,'ICC_D',.51528026061898091, ...
                'SEM_rel',1.9966198825315012,'SEM_abs',3.6291574654089972);
            f.nu(3) = struct('G',.76965090115835429,'D',.53289979436097916, ...
                'ICC_G',.51449170385221994,'ICC_D',.29063612019443668, ...
                'SEM_rel',2.7680382320435171,'SEM_abs',4.7370500082031387);
            f.nu(1) = struct('G',.87895111011115295,'D',.73735286475826926, ...
                'ICC_G',.66062042442123280,'ICC_D',.46100417521383989, ...
                'SEM_rel',2.0065927140704287,'SEM_abs',3.2270807531265433);
            f.nu(2) = struct('G',.84251770712688673,'D',.66022146472689312, ...
                'ICC_G',.66062042442123280,'ICC_D',.46100417521383989, ...
                'SEM_rel',2.2887306740036538,'SEM_abs',3.7977073414810376);
        end

        function f = localRetainedFixture()
            %Retained log two-way effects: (v_po, v_pi, v_oi) = (.09, .04, .0625).
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            f.s_po = .3; f.s_pt = .2; f.s_ot = .25;
            f.R  = 13.023130380021398;
            f.P  = 31.035097991267339; f.O = 20.341930048432374;
            f.I  = 20.341930048432374; f.OI = 13.393233403241014;
            f.PI = 14.194101887589867; f.PO = 21.228046894297880;
            f.C_ls = .04;                 f.E_ls = 13.063130380021398;
            f.C_nu = 11.096738624490647;  f.E_nu = 24.119869004512044;
            f.ls = struct('G',{.85772474919048112,.73657694108134588,.63913885146466610}, ...
                'D',{NaN,NaN,.45814016639152272});
            f.nu = struct('G',{NaN,NaN,.61577004934892007}, ...
                'D',{NaN,NaN,.44600731106885274});
        end

        function [pop_sdlog, ind_sdlog] = localStaticEncoding(R, C_typ, C_ss)
            %The two-line encoding both static extractors emit (mirrors
            %psyrat_gamma_extract_trt_sserr / _ls): pop_sdlog = .5 log(R+C_typ),
            %ind_sdlog = .5 log(R+C_p) - pop_sdlog. R, C_typ [nd x 1]; C_ss [nd x NSUB].
            pop_sdlog = 0.5 .* log(R + C_typ);
            ind_sdlog = 0.5 .* log(R + C_ss) - pop_sdlog;
        end

        function T = localStaticConsumer(vc, pop_sdlog, ind_sdlog, idtable, gcoeff, reltype, nocc)
            %psyrat_ssrel_trt fed exactly as psyrat_relsummary's
            %local_ssrel_trt_call feeds a gamma analysis-25 chunk.
            err_ss = exp(pop_sdlog + ind_sdlog);
            nsub = size(err_ss,2);
            T = psyrat_ssrel_trt('gcoeff',gcoeff,'reltype',reltype, ...
                'bp',sqrt(vc.sigma_p2),'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2), ...
                'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2),'txo',sqrt(vc.sigma_ot2), ...
                'err',exp(pop_sdlog), ...
                'txp_ss',repmat(sqrt(vc.sigma_pt2),1,nsub), ...
                'oxp_ss',repmat(sqrt(vc.sigma_po2),1,nsub), ...
                'err_ss',err_ss,'idtable',idtable,'CI',.95,'nocc',nocc);
        end

        function d = rep(v, n)
            d = repmat(v, n, 1);
        end

        function src = localSource(relpath)
            here = fileparts(mfilename('fullpath'));
            src = fileread(fullfile(fileparts(here), relpath));
        end

        function body = localFunctionBody(testCase, src, fname)
            starts = regexp(src, ['(?m)^function [^\n]*[ =]' fname '\s*\('], 'start');
            testCase.assertEqual(numel(starts), 1, ...
                sprintf('expected exactly one definition of %s', fname));
            rest = src(starts(1):end);
            nxt = regexp(rest(2:end), '(?m)^function ', 'start', 'once');
            if isempty(nxt); body = rest; else; body = rest(1:nxt); end
            body = regexprep(body, '(?m)%[^\n]*', '');      % comments cannot satisfy a pin
            body = regexprep(body, '\.\.\.\s*\n\s*', ' ');  % join continuations
        end

    end

    methods (Test)

        % ---------- the remainder itself ----------

        function testRemainderMatchesFrozenMainAnchor(testCase)
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            R = psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot);
            testCase.verifyEqual(R, f.R, 'RelTol', 1e-12, ...
                'R does not match the frozen main-fixture remainder (memo eq. 3).');
            % inclusion-exclusion closed form, a different algorithm
            v = struct('P',f.s_p^2,'O',f.s_o^2,'I',f.s_t^2,'PO',0,'PI',0,'OI',0);
            testCase.verifyEqual(R, ...
                TestGammaTwoFacetResidualPooling.localClosedFormRemainder(f.alpha, v), ...
                'RelTol', 1e-12);
            % with zero two-way effects R = B q_p q_o q_i (contract line 17)
            B = exp(2*f.alpha + f.s_p^2 + f.s_o^2 + f.s_t^2);
            testCase.verifyEqual(R, B*expm1(f.s_p^2)*expm1(f.s_o^2)*expm1(f.s_t^2), ...
                'RelTol', 1e-12);
            % and it is exactly what the production converter leaves inside
            % sigma_res2 by subtraction (up to that subtraction's rounding)
            vc = psyrat_gamma_varcomps_trt(f.alpha, f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, exp(f.lognu), 0, 0);
            testCase.verifyEqual(R, vc.sigma_res2 - vc.e_cond_var, 'RelTol', 1e-10, ...
                'R must be the subtraction remainder of psyrat_gamma_varcomps_trt.');
            testCase.verifyGreaterThan(R, 0);
        end

        function testRemainderWithRetainedTwoWayEffectsAndSeparatePiPo(testCase)
            f = TestGammaTwoFacetResidualPooling.localRetainedFixture();
            R = psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot);
            testCase.verifyEqual(R, f.R, 'RelTol', 1e-12, ...
                'R does not match the frozen retained-two-way remainder.');
            v = struct('P',f.s_p^2,'O',f.s_o^2,'I',f.s_t^2, ...
                'PO',f.s_po^2,'PI',f.s_pt^2,'OI',f.s_ot^2);
            testCase.verifyEqual(R, ...
                TestGammaTwoFacetResidualPooling.localClosedFormRemainder(f.alpha, v), ...
                'RelTol', 1e-12);
            % the six components are retained as their own orthogonal supports:
            % PI and PO do not absorb R, and match the frozen values
            vc = psyrat_gamma_varcomps_trt(f.alpha, f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, exp(f.lognu), 0, 0);
            c = TestGammaTwoFacetResidualPooling.localComponents(vc);
            for nm = {'P','O','I','PO','PI','OI'}
                testCase.verifyEqual(c.(nm{1}), f.(nm{1}), 'RelTol', 1e-12, ...
                    sprintf('%s component differs from the frozen value.', nm{1}));
            end
        end

        function testRemainderVanishesOnlyWhenAWholeFacetIsRemoved(testCase)
            f = TestGammaTwoFacetResidualPooling.localRetainedFixture();
            % occasion facet removed entirely: v_o = v_po = v_oi = 0 -> R = 0
            % EXACTLY (every product whose union covers the occasion carries a
            % zero factor); PI survives as its own component (frozen
            % complete_occasion)
            R = psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, 0, f.s_t, 0, f.s_pt, 0);
            testCase.verifyEqual(R, 0, ...
                'R must vanish exactly when the occasion facet is removed.');
            vc = psyrat_gamma_varcomps_trt(f.alpha, f.s_p, 0, f.s_t, 0, f.s_pt, 0, ...
                exp(f.lognu), 0, 0);
            testCase.verifyEqual(vc.sigma_pt2, 9.4908342614818533, 'RelTol', 1e-12, ...
                'PI must survive on its own when R vanishes.');
            testCase.verifyEqual(vc.sigma_p2, 20.751504650075511, 'RelTol', 1e-12);
            testCase.verifyEqual(vc.sigma_t2, 13.601557053576340, 'RelTol', 1e-12);
            % trial facet removed entirely: the symmetric identity
            R = psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, f.s_o, 0, f.s_po, 0, 0);
            testCase.verifyEqual(R, 0, ...
                'R must vanish exactly when the trial facet is removed.');
            % NOT a zero-facet control: only the occasion MAIN effect zeroed
            % while v_po and v_oi stay (frozen only_main_occasion_zero)
            R = psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, 0, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot);
            testCase.verifyEqual(R, 5.4537081777759397, 'RelTol', 1e-12, ...
                'Zeroing only the occasion main effect must not zero R.');
            vc = psyrat_gamma_varcomps_trt(f.alpha, f.s_p, 0, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, exp(f.lognu), 0, 0);
            testCase.verifyEqual(vc.sigma_po2, 7.5290411278020280, 'RelTol', 1e-12);
            % and with no facets at all, exactly zero
            testCase.verifyEqual(psyrat_gamma_mean_remainder_trt(f.alpha, f.s_p, 0, 0, 0, 0, 0), 0);
        end

        function testRemainderIsAccurateWhereSubtractionIsNot(testCase)
            % additive limit: log variances v_A / L^2 with L = 1000. R is ~1e-8
            % against var_mu ~ .4, so the converter's six-way subtraction (and
            % the inclusion-exclusion bracket) cancel to a few digits at best;
            % the enumeration must agree with the exact zero-two-way product
            % B q_p q_o q_i to machine precision, and the subtraction must not.
            L = 1000;
            alpha = log(L);
            sp = sqrt(.36)/L; so = sqrt(.25)/L; st = sqrt(.25)/L;
            R = psyrat_gamma_mean_remainder_trt(alpha, sp, so, st, 0, 0, 0);
            B = exp(2*alpha + sp^2 + so^2 + st^2);
            Rref = B*expm1(sp^2)*expm1(so^2)*expm1(st^2);
            testCase.verifyEqual(R, Rref, 'RelTol', 1e-12);
            testCase.verifyGreaterThan(R, 0);
            vc = psyrat_gamma_varcomps_trt(alpha, sp, so, st, 0, 0, 0, 18, 0, 0);
            Rsub = vc.var_mu - vc.sigma_p2 - vc.sigma_o2 - vc.sigma_t2 ...
                - vc.sigma_po2 - vc.sigma_pt2 - vc.sigma_ot2;
            testCase.verifyGreaterThan(abs(Rsub - Rref)/Rref, 1e-6, ...
                ['positive control: the subtraction remainder is NOT accurate ', ...
                'here, which is why the helper enumerates instead']);
            % draw-vector / scalar broadcasting and shape contract
            Rv = psyrat_gamma_mean_remainder_trt([alpha; alpha], sp, so, st, 0, 0, 0);
            testCase.verifyEqual(Rv, [R; R]);
            testCase.verifyError(@() psyrat_gamma_mean_remainder_trt([1;2], [1;2;3], 0, 0, 0, 0, 0), ...
                'psyrat_gamma_mean_remainder_trt:draws');
        end

        % ---------- static 25 encoding through the production consumer ----------

        function testStaticLogNuEncodingReproducesFrozenPersonMetrics(testCase)
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            vc = psyrat_gamma_varcomps_trt(rp(f.alpha), f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, exp(f.lognu), 0, 0);
            R = psyrat_gamma_mean_remainder_trt(rp(f.alpha), f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot);
            V_w = f.s_o^2 + f.s_t^2 + f.s_po^2 + f.s_pt^2 + f.s_ot^2;
            % subject 1: zero person effects (the frozen "person" target);
            % subject 2: own (u_p, v_p) = (.3, -.2)
            u = [0 .3]; v = [0 -.2];
            C_ss  = 2 .* exp(2.*(rp(f.alpha) + u) + 2.*V_w - (f.lognu + v));
            C_typ = 2 .* exp(2.*rp(f.alpha) + 2.*V_w - f.lognu);
            testCase.verifyEqual(C_ss(1,1), f.C_nu, 'RelTol', 1e-12);
            [pop_sdlog, ind_sdlog] = ...
                TestGammaTwoFacetResidualPooling.localStaticEncoding(R, C_typ, C_ss);
            % the encoding reproduces E_p = R + C_p once, per subject
            testCase.verifyEqual(exp(pop_sdlog + ind_sdlog).^2, R + C_ss, 'RelTol', 1e-13);
            testCase.verifyEqual(exp(2*pop_sdlog), R + C_typ, 'RelTol', 1e-13);

            idtable = table({'a';'b'}, [1;2], [f.n_i; f.n_i], ...
                'VariableNames', {'id','id2','trls'});
            c = TestGammaTwoFacetResidualPooling.localComponents(structfun(@(x) x(1), vc, 'UniformOutput', false));
            for reltype = 1:3
                Tg = TestGammaTwoFacetResidualPooling.localStaticConsumer(vc, pop_sdlog, ind_sdlog, idtable, 2, reltype, f.n_o);
                Td = TestGammaTwoFacetResidualPooling.localStaticConsumer(vc, pop_sdlog, ind_sdlog, idtable, 1, reltype, f.n_o);
                % subject 1 against the frozen person-target metrics
                fr = f.nu(reltype);
                testCase.verifyEqual(Tg.dep_pt(1), fr.G,     'RelTol', 1e-12, sprintf('reltype %d G', reltype));
                testCase.verifyEqual(Td.dep_pt(1), fr.D,     'RelTol', 1e-12, sprintf('reltype %d D', reltype));
                testCase.verifyEqual(Tg.icc_pt(1), fr.ICC_G, 'RelTol', 1e-12, sprintf('reltype %d ICC_G', reltype));
                testCase.verifyEqual(Td.icc_pt(1), fr.ICC_D, 'RelTol', 1e-12, sprintf('reltype %d ICC_D', reltype));
                testCase.verifyEqual(Tg.sem_pt(1), fr.SEM_rel, 'RelTol', 1e-12, sprintf('reltype %d SEM rel', reltype));
                testCase.verifyEqual(Td.sem_pt(1), fr.SEM_abs, 'RelTol', 1e-12, sprintf('reltype %d SEM abs', reltype));
                testCase.verifyEqual(Tg.ss_errvar(1), f.E_nu, 'RelTol', 1e-12, 'ss_errvar must be R + C_p');
                % subject 2 against the contract arithmetic with its own C_p
                m = TestGammaTwoFacetResidualPooling.localContractMetrics(c, R(1) + C_ss(1,2), f.n_i, f.n_o, reltype);
                testCase.verifyEqual(Tg.dep_pt(2), m.G, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.dep_pt(2), m.D, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.icc_pt(2), m.ICC_G, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.icc_pt(2), m.ICC_D, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.sem_pt(2), m.SEM_rel, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.sem_pt(2), m.SEM_abs, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.ss_errvar(2), R(1) + C_ss(1,2), 'RelTol', 1e-12);
                % the reference line (err = exp(pop_sdlog)) is R + C_typ, i.e.
                % the typical person, not the marginal group residual
                testCase.verifyEqual(Tg.gro_rel(1), Tg.dep_pt(1), 'RelTol', 1e-12, ...
                    'gro_rel at equal counts must sit on the typical-person point.');
            end
            % the person coefficient depends on the own effects (subject 2 differs)
            testCase.verifyGreaterThan(abs(Tg.dep_pt(2) - Tg.dep_pt(1)), 1e-3);
        end

        function testStaticLocationScaleEncodingReproducesFrozenPersonMetrics(testCase)
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            vc = psyrat_gamma_varcomps_trt(rp(f.alpha), f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, NaN, NaN, NaN);   % the LS extractor's NaN sentinel
            R = psyrat_gamma_mean_remainder_trt(rp(f.alpha), f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot);
            w = [0 .45];
            C_ss  = exp(2 .* (rp(f.logsig) + w));
            C_typ = exp(2 .* rp(f.logsig));
            [pop_sdlog, ind_sdlog] = ...
                TestGammaTwoFacetResidualPooling.localStaticEncoding(R, C_typ, C_ss);
            testCase.verifyEqual(exp(pop_sdlog + ind_sdlog).^2, R + C_ss, 'RelTol', 1e-13);
            % NO own-location and NO exp(2*V_w) factor on the LS observation
            % term: the frozen residual is R + .04, which an exp(2*V_w) = e
            % factor (.109) or an amplitude factor would miss by a wide margin
            testCase.verifyEqual(R(1) + C_ss(1,1), f.E_ls, 'RelTol', 1e-12);

            idtable = table({'a';'b'}, [1;2], [f.n_i; f.n_i], ...
                'VariableNames', {'id','id2','trls'});
            c = TestGammaTwoFacetResidualPooling.localComponents(structfun(@(x) x(1), vc, 'UniformOutput', false));
            for reltype = 1:3
                Tg = TestGammaTwoFacetResidualPooling.localStaticConsumer(vc, pop_sdlog, ind_sdlog, idtable, 2, reltype, f.n_o);
                Td = TestGammaTwoFacetResidualPooling.localStaticConsumer(vc, pop_sdlog, ind_sdlog, idtable, 1, reltype, f.n_o);
                fr = f.ls(reltype);
                testCase.verifyEqual(Tg.dep_pt(1), fr.G,     'RelTol', 1e-12, sprintf('reltype %d G', reltype));
                testCase.verifyEqual(Td.dep_pt(1), fr.D,     'RelTol', 1e-12, sprintf('reltype %d D', reltype));
                testCase.verifyEqual(Tg.icc_pt(1), fr.ICC_G, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.icc_pt(1), fr.ICC_D, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.sem_pt(1), fr.SEM_rel, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.sem_pt(1), fr.SEM_abs, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.ss_errvar(1), f.E_ls, 'RelTol', 1e-12);
                m = TestGammaTwoFacetResidualPooling.localContractMetrics(c, R(1) + C_ss(1,2), f.n_i, f.n_o, reltype);
                testCase.verifyEqual(Tg.dep_pt(2), m.G, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.dep_pt(2), m.D, 'RelTol', 1e-12);
                testCase.verifyEqual(Tg.sem_pt(2), m.SEM_rel, 'RelTol', 1e-12);
                testCase.verifyEqual(Td.sem_pt(2), m.SEM_abs, 'RelTol', 1e-12);
            end
            % the previous (observation-only) LS encoding is the refuted
            % baseline: G .808441253962432 under C alone vs .799748348038617
            [p0, i0] = TestGammaTwoFacetResidualPooling.localStaticEncoding(zeros(nd,1), C_typ, C_ss);
            T0 = TestGammaTwoFacetResidualPooling.localStaticConsumer(vc, p0, i0, idtable, 2, 3, f.n_o);
            testCase.verifyEqual(T0.dep_pt(1), .808441253962432, 'RelTol', 1e-12, ...
                'positive control: the observation-only encoding reproduces the pre-Q8 value');
        end

        % ---------- dynamic 27 calculators ----------

        function testDynamicCalculatorsMatchFrozenAtZeroSlopes(testCase)
            % zero slopes and zero person effects: every participant sits on the
            % frozen person target whatever their z (own-z invariance at b = 0)
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            z = [-1.1; 0; .7];
            n = numel(z);
            ssinfo = table((1:n)', (1:n)', z, repmat(f.n_i, n, 1), repmat(f.n_o, n, 1), ...
                'VariableNames', {'id','id2','z1','trls','occs'});
            for reltype = 1:3
                [gen, dep] = psyrat_ssrel_dynrel_trt_gamma( ...
                    'alpha0', rp(f.alpha), 'b', rp(0), 'lognu0', rp(f.lognu), 'b_nu', rp(0), ...
                    'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                    'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                    'ind_bs', zeros(nd, n), 'ind_nu', zeros(nd, n), ...
                    'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, 'nocc', f.n_o, 'CI', .95);
                fr = f.nu(reltype);
                for k = 1:n
                    testCase.verifyEqual(gen.dep_pt(k), fr.G, 'RelTol', 1e-12, sprintf('lognu reltype %d G', reltype));
                    testCase.verifyEqual(dep.dep_pt(k), fr.D, 'RelTol', 1e-12, sprintf('lognu reltype %d D', reltype));
                    testCase.verifyEqual(gen.icc_pt(k), fr.ICC_G, 'RelTol', 1e-12);
                    testCase.verifyEqual(dep.icc_pt(k), fr.ICC_D, 'RelTol', 1e-12);
                    testCase.verifyEqual(gen.sem_pt(k), fr.SEM_rel, 'RelTol', 1e-12);
                    testCase.verifyEqual(dep.sem_pt(k), fr.SEM_abs, 'RelTol', 1e-12);
                    testCase.verifyEqual(gen.ss_errvar(k), f.E_nu, 'RelTol', 1e-12, 'lognu ss_errvar must be R + C_typ');
                end
                [gen, dep] = psyrat_ssrel_dynrel_trt_gamma_ls( ...
                    'alpha0', rp(f.alpha), 'b', rp(0), 'logsig0', rp(f.logsig), 'b_sigma', rp(0), ...
                    'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                    'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                    'ind_sd', zeros(nd, n), ...
                    'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, 'nocc', f.n_o, 'CI', .95);
                fr = f.ls(reltype);
                for k = 1:n
                    testCase.verifyEqual(gen.dep_pt(k), fr.G, 'RelTol', 1e-12, sprintf('LS reltype %d G', reltype));
                    testCase.verifyEqual(dep.dep_pt(k), fr.D, 'RelTol', 1e-12, sprintf('LS reltype %d D', reltype));
                    testCase.verifyEqual(gen.icc_pt(k), fr.ICC_G, 'RelTol', 1e-12);
                    testCase.verifyEqual(dep.icc_pt(k), fr.ICC_D, 'RelTol', 1e-12);
                    testCase.verifyEqual(gen.sem_pt(k), fr.SEM_rel, 'RelTol', 1e-12);
                    testCase.verifyEqual(dep.sem_pt(k), fr.SEM_abs, 'RelTol', 1e-12);
                    testCase.verifyEqual(gen.ss_errvar(k), f.E_ls, 'RelTol', 1e-12, 'LS ss_errvar must be R + C_typ');
                end
            end
        end

        function testDynamicCalculatorsRetainedTwoWayAnchor(testCase)
            f = TestGammaTwoFacetResidualPooling.localRetainedFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            ssinfo = table(1, 1, .3, f.n_i, f.n_o, 'VariableNames', {'id','id2','z1','trls','occs'});
            for reltype = 1:3
                [gen, dep] = psyrat_ssrel_dynrel_trt_gamma_ls( ...
                    'alpha0', rp(f.alpha), 'b', rp(0), 'logsig0', rp(f.logsig), 'b_sigma', rp(0), ...
                    'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                    'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                    'ind_sd', zeros(nd, 1), 'idtable', ssinfo, 'ndim', 1, ...
                    'reltype', reltype, 'nocc', f.n_o, 'CI', .95);
                testCase.verifyEqual(gen.dep_pt, f.ls(reltype).G, 'RelTol', 1e-12, sprintf('LS retained reltype %d', reltype));
                testCase.verifyEqual(gen.ss_errvar, f.E_ls, 'RelTol', 1e-12);
                if reltype == 3
                    testCase.verifyEqual(dep.dep_pt, f.ls(3).D, 'RelTol', 1e-12);
                end
            end
            [gen, dep] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', rp(f.alpha), 'b', rp(0), 'lognu0', rp(f.lognu), 'b_nu', rp(0), ...
                'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                'ind_bs', zeros(nd, 1), 'ind_nu', zeros(nd, 1), 'idtable', ssinfo, 'ndim', 1, ...
                'reltype', 3, 'nocc', f.n_o, 'CI', .95);
            testCase.verifyEqual(gen.dep_pt, f.nu(3).G, 'RelTol', 1e-12, 'lognu retained CES G');
            testCase.verifyEqual(dep.dep_pt, f.nu(3).D, 'RelTol', 1e-12, 'lognu retained CES D');
            testCase.verifyEqual(gen.ss_errvar, f.E_nu, 'RelTol', 1e-12);
        end

        function testDynamicOwnZAndIdPermutation(testCase)
            % nonzero slopes, permuted id2 order, distinct z and own effects: row
            % k must use alpha(z_k), R(z_k), the components at z_k and the draw
            % column id2(k). Expected values from the contract arithmetic with
            % the union-enumeration remainder, not from psyrat_rel_trt.
            f = TestGammaTwoFacetResidualPooling.localRetainedFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            b = .2; b_nu = -.25; b_sigma = .15;
            id2  = [3; 1; 2];
            z    = [.5; -.8; 1.1];
            trls = [9; 4; 6];
            n_o  = 2;
            u = [.25 -.10 .05];   % columns indexed by id2
            v = [-.20 .30 .00];
            w = [.30 -.15 .05];
            ssinfo = table((101:103)', id2, z, trls, repmat(n_o, 3, 1), ...
                'VariableNames', {'id','id2','z1','trls','occs'});
            V_w = f.s_o^2 + f.s_t^2 + f.s_po^2 + f.s_pt^2 + f.s_ot^2;
            vlog = struct('P',f.s_p^2,'O',f.s_o^2,'I',f.s_t^2, ...
                'PO',f.s_po^2,'PI',f.s_pt^2,'OI',f.s_ot^2);
            for reltype = [1 3]
                [genN, depN] = psyrat_ssrel_dynrel_trt_gamma( ...
                    'alpha0', rp(f.alpha), 'b', rp(b), 'lognu0', rp(f.lognu), 'b_nu', rp(b_nu), ...
                    'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                    'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                    'ind_bs', repmat(u, nd, 1), 'ind_nu', repmat(v, nd, 1), ...
                    'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, 'nocc', n_o, 'CI', .95);
                [genL, depL] = psyrat_ssrel_dynrel_trt_gamma_ls( ...
                    'alpha0', rp(f.alpha), 'b', rp(b), 'logsig0', rp(f.logsig), 'b_sigma', rp(b_sigma), ...
                    'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                    'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                    'ind_sd', repmat(w, nd, 1), ...
                    'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, 'nocc', n_o, 'CI', .95);
                for k = 1:3
                    col = id2(k);
                    alpha_z  = f.alpha + b * z(k);
                    lognu_z  = f.lognu + b_nu * z(k);
                    logsig_z = f.logsig + b_sigma * z(k);
                    vc = psyrat_gamma_varcomps_trt(alpha_z, f.s_p, f.s_o, f.s_t, ...
                        f.s_po, f.s_pt, f.s_ot, exp(lognu_z), 0, 0);
                    c = TestGammaTwoFacetResidualPooling.localComponents(vc);
                    R_z = TestGammaTwoFacetResidualPooling.localClosedFormRemainder(alpha_z, vlog);
                    E_nu = R_z + 2*exp(2*(alpha_z + u(col)) + 2*V_w - (lognu_z + v(col)));
                    E_ls = R_z + exp(2*(logsig_z + w(col)));
                    mN = TestGammaTwoFacetResidualPooling.localContractMetrics(c, E_nu, trls(k), n_o, reltype);
                    mL = TestGammaTwoFacetResidualPooling.localContractMetrics(c, E_ls, trls(k), n_o, reltype);
                    testCase.verifyEqual(genN.dep_pt(k), mN.G, 'RelTol', 1e-12, sprintf('lognu row %d G', k));
                    testCase.verifyEqual(depN.dep_pt(k), mN.D, 'RelTol', 1e-12, sprintf('lognu row %d D', k));
                    testCase.verifyEqual(genN.icc_pt(k), mN.ICC_G, 'RelTol', 1e-12);
                    testCase.verifyEqual(depN.icc_pt(k), mN.ICC_D, 'RelTol', 1e-12);
                    testCase.verifyEqual(genN.sem_pt(k), mN.SEM_rel, 'RelTol', 1e-12);
                    testCase.verifyEqual(depN.sem_pt(k), mN.SEM_abs, 'RelTol', 1e-12);
                    testCase.verifyEqual(genN.ss_errvar(k), E_nu, 'RelTol', 1e-12);
                    testCase.verifyEqual(genL.dep_pt(k), mL.G, 'RelTol', 1e-12, sprintf('LS row %d G', k));
                    testCase.verifyEqual(depL.dep_pt(k), mL.D, 'RelTol', 1e-12, sprintf('LS row %d D', k));
                    testCase.verifyEqual(genL.icc_pt(k), mL.ICC_G, 'RelTol', 1e-12);
                    testCase.verifyEqual(depL.icc_pt(k), mL.ICC_D, 'RelTol', 1e-12);
                    testCase.verifyEqual(genL.sem_pt(k), mL.SEM_rel, 'RelTol', 1e-12);
                    testCase.verifyEqual(depL.sem_pt(k), mL.SEM_abs, 'RelTol', 1e-12);
                    testCase.verifyEqual(genL.ss_errvar(k), E_ls, 'RelTol', 1e-12);
                end
                % row order is preserved (id column) and columns were permuted
                testCase.verifyEqual(genN.id, (101:103)');
                testCase.verifyEqual(genN.id2, id2);
            end
        end

        function testResidualIsNeitherGroupResidualNorDoublePooled(testCase)
            % at zero person effects with s_v = 0 the residual is R + C_typ:
            % strictly between the observation-only C_typ and the marginal group
            % residual sigma_res2 = R + C_pop (C_pop carries exp(2 v_p) over
            % C_typ), and not 2R + C.
            f = TestGammaTwoFacetResidualPooling.localMainFixture();
            nd = 4;
            rp = @(x) TestGammaTwoFacetResidualPooling.rep(x, nd);
            ssinfo = table(1, 1, 0, f.n_i, f.n_o, 'VariableNames', {'id','id2','z1','trls','occs'});
            [gen, ~] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', rp(f.alpha), 'b', rp(0), 'lognu0', rp(f.lognu), 'b_nu', rp(0), ...
                'sig_p', rp(f.s_p), 'sig_occ', rp(f.s_o), 'sig_trl', rp(f.s_t), ...
                'sig_occxid', rp(f.s_po), 'sig_trlxid', rp(f.s_pt), 'sig_trlxocc', rp(f.s_ot), ...
                'ind_bs', zeros(nd,1), 'ind_nu', zeros(nd,1), 'idtable', ssinfo, 'ndim', 1, ...
                'reltype', 3, 'nocc', f.n_o, 'CI', .95);
            vc = psyrat_gamma_varcomps_trt(f.alpha, f.s_p, f.s_o, f.s_t, ...
                f.s_po, f.s_pt, f.s_ot, exp(f.lognu), 0, 0);
            E = gen.ss_errvar;
            testCase.verifyEqual(E, f.R + f.C_nu, 'RelTol', 1e-12);
            testCase.verifyGreaterThan(E, f.C_nu * (1 + 1e-3), 'the remainder must be present');
            testCase.verifyLessThan(E, vc.sigma_res2 * (1 - 1e-3), ...
                'the entire marginal group residual must not be used');
            testCase.verifyLessThan(E, 2*f.R + f.C_nu - 1e-3, 'R must be pooled once');
        end

        % ---------- static extractor wiring and v2 stamp (source pins) ----------

        function testStaticExtractorsAreWiredToTheHelperAndEncodeRPlusC(testCase)
            % The two extractors cannot be executed here; pin their wiring to
            % the helper and the encoding form. The numerical content of that
            % encoding is exercised above through psyrat_ssrel_trt.
            src = TestGammaTwoFacetResidualPooling.localSource('subroutines/estimation/psyrat_computevarcomp.m');
            nu = TestGammaTwoFacetResidualPooling.localFunctionBody(testCase, src, 'psyrat_gamma_extract_trt_sserr');
            ls = TestGammaTwoFacetResidualPooling.localFunctionBody(testCase, src, 'psyrat_gamma_extract_trt_sserr_ls');
            for body = {nu, ls}
                testCase.verifyTrue(~isempty(regexp(body{1}, ...
                    'R\s*=\s*psyrat_gamma_mean_remainder_trt\(alpha,\s*s_p,\s*s_o,\s*s_t,\s*s_po,\s*s_pt,\s*s_ot\)', 'once')), ...
                    'extractor must build R from alpha and the six log-mean SDs via the helper');
                testCase.verifyTrue(isempty(regexp(body{1}, 'psyrat_gamma_mean_remainder_trt\([^)]*ind_', 'once')), ...
                    'R must never be built from a person effect');
            end
            % log-nu: E_p = R + C_p, reference R + C_typ, old single-exponential
            % residual gone
            testCase.verifyTrue(~isempty(regexp(nu, 'sig_pot_e2\s*=\s*R\s*\+\s*sig_pot_c2\s*;', 'once')));
            testCase.verifyTrue(~isempty(regexp(nu, 'sig_pot_e2_pop\s*=\s*R\s*\+\s*sig_pot_c2_typ\s*;', 'once')));
            testCase.verifyTrue(isempty(regexp(nu, 'sig_pot_e2\s*=\s*2\s*\.\*\s*exp\(', 'once')), ...
                'the observation-only log-nu residual must be gone');
            testCase.verifyTrue(~isempty(regexp(nu, 'pop_sdlog\s*=\s*0\.5\s*\.\*\s*log\(sig_pot_e2_pop\)', 'once')));
            testCase.verifyTrue(~isempty(regexp(nu, 'ind_sdlog\s*=\s*0\.5\s*\.\*\s*log\(sig_pot_e2\)\s*-\s*pop_sdlog', 'once')));
            % location-scale: composite encoding, pass-through gone
            testCase.verifyTrue(~isempty(regexp(ls, 'pop_sdlog\s*=\s*0\.5\.\*log\(R\s*\+\s*exp\(2\.\*logsig\(:\)\)\)', 'once')));
            testCase.verifyTrue(~isempty(regexp(ls, ...
                'ind_sdlog\s*=\s*0\.5\.\*log\(R\s*\+\s*exp\(2\.\*\(logsig\(:\)\s*\+\s*ind_sd\)\)\)\s*-\s*pop_sdlog', 'once')));
            testCase.verifyTrue(isempty(regexp(ls, 'ind_sdlog\s*=\s*ind_sd\s*;', 'once')), ...
                'the raw log-sigma pass-through must be gone');
            testCase.verifyTrue(isempty(regexp(ls, '(pop_sdlog|ind_sdlog|R)\s*=[^;]*ind_bs', 'once')), ...
                'no LS residual-channel line may reference ind_bs');
            % raw scale pieces preserved on the provenance struct
            for body = {nu, ls}
                for fld = {'mean_remainder', 'obs_var_typ', 'obs_var_ss', 'scale0', 'ind_scale', 'parameterization'}
                    testCase.verifyTrue(~isempty(regexp(body{1}, ['''' fld{1} ''''], 'once')), ...
                        sprintf('provenance field %s missing from an extractor', fld{1}));
                end
            end
            % positive control: the executable dynamic calculators carry the same call
            dyn = TestGammaTwoFacetResidualPooling.localSource('subroutines/calculation/psyrat_ssrel_dynrel_trt_gamma.m');
            testCase.verifyTrue(~isempty(regexp(dyn, 'psyrat_gamma_mean_remainder_trt\(alpha_z,', 'once')));
        end

        function testCase25StoresTheVersionedTargetAndPerStratumRawPieces(testCase)
            src = TestGammaTwoFacetResidualPooling.localSource('subroutines/estimation/psyrat_computevarcomp.m');
            code = regexprep(src, '(?m)%[^\n]*', '');
            code = regexprep(code, '\.\.\.\s*\n\s*', ' ');
            % the case-25 block runs from its case line to the next 'case 7'
            % (\s, not \b: MATLAB's \b is a backspace, not a word boundary)
            s = regexp(code, '(?m)^\s*case 25\s', 'start', 'once');
            e = regexp(code(s:end), '(?m)^\s*case 7\s', 'start', 'once');
            testCase.assertTrue(~isempty(s) && ~isempty(e), 'case 25 block not located');
            blk = code(s:s+e);
            testCase.verifyTrue(~isempty(regexp(blk, 'REL\.out\.gamma_target\s*=\s*struct\(', 'once')));
            testCase.verifyTrue(~isempty(regexp(blk, '''gamma_target_version'',\s*2\s*,', 'once')), ...
                'gamma_target_version must be 2');
            testCase.verifyTrue(~isempty(regexp(blk, ...
                '''remainder_scheme'',\s*''population_highest_order_plus_person_observation_v2''', 'once')));
            for fld = {'subject_target','reference_target','parameterization','count_source','analysis'}
                testCase.verifyTrue(~isempty(regexp(blk, ['''' fld{1} ''','], 'once')), ...
                    sprintf('gamma_target field %s missing', fld{1}));
            end
            testCase.verifyTrue(~isempty(regexp(blk, 'REL\.out\.gamma_resid\{end\+1\}\s*=\s*struct\(', 'once')), ...
                'raw residual pieces must be stored per stratum');
            for fld = {'mean_remainder','obs_var_typ','obs_var_ss','scale0','ind_scale','parameterization'}
                testCase.verifyTrue(~isempty(regexp(blk, ['''' fld{1} ''',\s*gout\.gamma\.' fld{1}], 'once')), ...
                    sprintf('gamma_resid field %s missing', fld{1}));
            end
            % the stamp is gamma-only (a Gaussian analysis-25 result gets none)
            testCase.verifyTrue(~isempty(regexp(blk, ...
                'if strcmp\(family,''gamma''\)\s+if gammascale == 2\s+gparam = ''ls'';', 'once')));
        end

    end
end
