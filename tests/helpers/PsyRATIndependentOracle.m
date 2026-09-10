classdef PsyRATIndependentOracle
    %PsyRATIndependentOracle  Paper-sourced reliability formulas, derived WITHOUT reading
    %                         PsyRAT's implementation.
    %
    %WHY THIS CLASS EXISTS, AND HOW IT DIFFERS FROM PsyRATAccuracyOracle
    %   PsyRATAccuracyOracle is a REGRESSION MIRROR: its expressions are the same algebra
    %   as the production functions, so a green run there proves "the formulas are
    %   unchanged", never "the formulas are correct". Both would have to be wrong in the
    %   same way for it to notice, which is exactly what a duplicated formula guarantees.
    %
    %   This class is different in provenance, which is the only thing that matters for
    %   this purpose. Every expression below was transcribed from the published tables of
    %   Rocha et al. (2026) by a reader working from the page ALONE, with no access to the
    %   PsyRAT source, and with the paper's own symbols rather than PsyRAT's argument
    %   names. Tables 2 and 3 were each transcribed TWICE from independent renderings of
    %   the page (a text extraction and a 300-dpi image), and the two transcriptions agreed
    %   before either was compared against any code.
    %
    %   PROVENANCE OF THE TABLE 6 BLOCK IS WEAKER, AND IS RECORDED HONESTLY RATHER THAN
    %   ASSIMILATED TO THE ABOVE. Table 6 was also transcribed twice from independent
    %   renderings (the text layer and a 400-dpi image of p.13) and the two agreed before
    %   either met any code, so the double-rendering check is the same. What differs is the
    %   reader: Tables 2 and 3 were transcribed by someone who never read the source, while
    %   Table 6 was transcribed by a reader who went on to read psyrat_diffrel in the same
    %   session. The protection there is ORDERING, not ignorance -- the transcription was
    %   written down and frozen before the production algebra was opened, so it could not be
    %   bent to match. That is strictly weaker than a source-blind reader. Treat the Table 6
    %   block as independent-by-ordering and, if the difference-score coefficients ever
    %   become contested, re-transcribe it blind.
    %
    %   Consequence: agreement between this class and production is EXTERNAL corroboration.
    %   Disagreement is a finding to adjudicate -- possibly a deliberate estimand choice,
    %   possibly a transcription error here, possibly a defect there -- and must never be
    %   "fixed" by editing this file to match the code. Doing that would silently convert
    %   this class into a second mirror and destroy the only independent check in the repo.
    %
    %CONVENTIONS (deliberately NOT PsyRAT's)
    %   * Inputs are VARIANCES, named for the paper's components (sigma2_p, sigma2_i, ...).
    %     PsyRAT's own functions take STANDARD DEVIATIONS and square them internally. The
    %     difference is intentional: a caller must consciously convert, so a unit error
    %     cannot pass silently from one to the other.
    %   * Replicate counts are nprime_i and nprime_o (the paper's n'_i and n'_o).
    %   * No posterior draws, no credible intervals, no (ll, pt, ul) triple. These are the
    %     closed-form population quantities only.
    %
    %NOTATION
    %   sigma2_p     person (universe-score) variance
    %   sigma2_i     trial main effect
    %   sigma2_o     occasion main effect
    %   sigma2_pi    person x trial interaction        (two-facet)
    %   sigma2_po    person x occasion interaction     (two-facet)
    %   sigma2_oi    occasion x trial interaction, NO person term (two-facet)
    %   sigma2_pi_e  person x trial + undifferentiated error (one-facet; not separately
    %                identified in that design, hence the single symbol)
    %   sigma2_poi_e three-way interaction + undifferentiated error (two-facet)
    %   mu           grand mean of the universe score
    %   C            cut score (criterion)
    %
    %   Difference scores (Table 6) index the two constituent measures X and Y:
    %   sigma2_pX    person variance of measure X          (likewise sigma2_pY)
    %   sigmaXY_p    person COVARIANCE between X and Y     (a covariance, not a variance)
    %   nprime_iX    trial replicates for measure X        (likewise nprime_iY, nprime_oX,
    %                                                       nprime_oY)
    %
    %THE DISTINCTION THIS CLASS EXISTS TO PIN DOWN
    %   Relative error contains only components carrying a PERSON index. Absolute error
    %   adds every component that does NOT carry one -- the facet main effects and the
    %   facet-only interactions. The reason is not algebraic convenience: a component with
    %   no person index shifts every person's score by the same amount, which is invisible
    %   to a comparison AMONG persons (relative) but decides where a score falls on the
    %   scale itself (absolute). Hence:
    %       generalizability (G) uses RELATIVE error
    %       dependability    (D) uses ABSOLUTE error
    %       cut-scores       are ALWAYS absolute -- see the note on cutScoreOneFacet.
    %
    %SOURCE
    %   Rocha, H. A., et al. (2026). International Journal of Psychophysiology, 222,
    %   113321. Table 2 (p. 7), Table 3 (p. 8), Table 6 (p. 13).
    %
    %   Copyright (C) 2016-2025 Peter E. Clayson
    %
    %   This program is free software: you can redistribute it and/or modify it under the
    %   terms of the GNU General Public License as published by the Free Software
    %   Foundation, either version 3 of the License, or any later version.
    %
    %   This program is distributed in the hope that it will be useful, but WITHOUT ANY
    %   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    %   PARTICULAR PURPOSE. See the GNU General Public License for more details.
    %
    %   You should have received a copy of the GNU General Public License along with this
    %   program (gpl.txt). If not, see <http://www.gnu.org/licenses/>.

    methods (Static)

        %% ================================================================ TABLE 2
        %  One-facet design, Persons x Trials. Rocha 2026, p. 7.

        function v = errVarRelativeOneFacet(sigma2_pi_e, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Error Variances / Relative".
            % Only the person-carrying residual. The trial main effect is absent.
            v = sigma2_pi_e ./ nprime_i;
        end

        function v = errVarAbsoluteOneFacet(sigma2_pi_e, sigma2_i, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Error Variances / Absolute".
            % Adds sigma2_i/nprime_i -- the trial main effect, divided by the SAME
            % replicate count as the residual, not by 1.
            v = sigma2_pi_e ./ nprime_i + sigma2_i ./ nprime_i;
        end

        function v = semRelativeOneFacet(sigma2_pi_e, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Standard Error of Measurement / Relative".
            v = sqrt(PsyRATIndependentOracle.errVarRelativeOneFacet(sigma2_pi_e, nprime_i));
        end

        function v = semAbsoluteOneFacet(sigma2_pi_e, sigma2_i, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Standard Error of Measurement / Absolute".
            v = sqrt(PsyRATIndependentOracle.errVarAbsoluteOneFacet(sigma2_pi_e, sigma2_i, nprime_i));
        end

        function v = gCoefficientOneFacet(sigma2_p, sigma2_pi_e, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Global coefficients / G-coefficient".
            % GENERALIZABILITY: relative error. sigma2_i does not appear.
            v = sigma2_p ./ (sigma2_p + sigma2_pi_e ./ nprime_i);
        end

        function v = dCoefficientOneFacet(sigma2_p, sigma2_pi_e, sigma2_i, nprime_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Global coefficients / D-Coefficient".
            % DEPENDABILITY: absolute error. Same numerator as G, strictly larger
            % denominator, so D <= G always, with equality only when sigma2_i = 0.
            v = sigma2_p ./ (sigma2_p + sigma2_pi_e ./ nprime_i + sigma2_i ./ nprime_i);
        end

        function v = iccRelativeOneFacet(sigma2_p, sigma2_pi_e)
            % SOURCE: Rocha 2026 p.7 Table 2, row "ICC / Relative".
            % Note the printed cell carries NO divisor: this is the single-observation
            % (n' = 1) form, not the G-coefficient with nprime_i substituted.
            v = sigma2_p ./ (sigma2_p + sigma2_pi_e);
        end

        function v = iccAbsoluteOneFacet(sigma2_p, sigma2_pi_e, sigma2_i)
            % SOURCE: Rocha 2026 p.7 Table 2, row "ICC / Absolute".
            v = sigma2_p ./ (sigma2_p + sigma2_pi_e + sigma2_i);
        end

        function v = cutScoreOneFacet(sigma2_p, sigma2_pi_e, sigma2_i, nprime_i, mu, C)
            % SOURCE: Rocha 2026 p.7 Table 2, row "Cut-Score".
            %
            % Table 2 prints exactly ONE Cut-Score row and it is UNPAIRED: every other
            % coefficient row in the table is a two-cell pair (Relative | Absolute, or
            % G-coefficient | D-Coefficient), and the cut-score's relative cell is blank.
            % Its denominator error term is character-for-character the table's own
            % "Error Variances / Absolute" cell. So a relative-error cut score is not
            % merely unimplemented here -- the paper does not define one.
            %
            % (mu - C)^2 is added to BOTH numerator and denominator, which is why the
            % coefficient rises toward 1 as the cut score moves away from the mean, and
            % why at mu == C it reduces exactly to the D-coefficient.
            offsq = (mu - C) .^ 2;
            v = (sigma2_p + offsq) ./ ...
                (sigma2_p + offsq + ...
                 PsyRATIndependentOracle.errVarAbsoluteOneFacet(sigma2_pi_e, sigma2_i, nprime_i));
        end

        %% ================================================================ TABLE 3
        %  Two-facet design, Persons x Trials x Occasions. Rocha 2026, p. 8.
        %
        %  Three coefficient types, distinguished by which facet is treated as FIXED:
        %    CE  coefficient of equivalence                   -- occasion fixed
        %    CS  coefficient of stability                     -- trial fixed
        %    CES coefficient of trial equivalence and stability -- both random
        %
        %  A fixed facet moves its person-by-facet interaction OUT of the error term and
        %  INTO the universe score, and drops that facet's own main effect from the
        %  absolute error. That is why CE_D omits sigma2_o and CS_D omits sigma2_i, while
        %  CES_D carries the full absolute error. Those two omissions are printed in the
        %  table and are NOT transcription slips.

        function v = errVarRelativeTwoFacet(sigma2_pi, sigma2_po, sigma2_poi_e, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, row "Error Variances / Relative".
            v = sigma2_pi ./ nprime_i + sigma2_po ./ nprime_o + ...
                sigma2_poi_e ./ (nprime_o .* nprime_i);
        end

        function v = errVarAbsoluteTwoFacet(sigma2_pi, sigma2_po, sigma2_poi_e, ...
                                            sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, row "Error Variances / Absolute".
            % Adds exactly the three components with no person index.
            v = PsyRATIndependentOracle.errVarRelativeTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, nprime_i, nprime_o) + ...
                sigma2_o ./ nprime_o + sigma2_i ./ nprime_i + ...
                sigma2_oi ./ (nprime_o .* nprime_i);
        end

        function v = semRelativeTwoFacet(sigma2_pi, sigma2_po, sigma2_poi_e, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, row "Standard Error of Measurement / Relative".
            v = sqrt(PsyRATIndependentOracle.errVarRelativeTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, nprime_i, nprime_o));
        end

        function v = semAbsoluteTwoFacet(sigma2_pi, sigma2_po, sigma2_poi_e, ...
                                         sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, row "Standard Error of Measurement / Absolute".
            v = sqrt(PsyRATIndependentOracle.errVarAbsoluteTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o));
        end

        % ---- CE: coefficient of equivalence (occasion FIXED) --------------------
        function v = ceG(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CE", row "G-coefficient".
            % Occasion fixed => sigma2_po/nprime_o joins the UNIVERSE SCORE.
            v = (sigma2_p + sigma2_po ./ nprime_o) ./ ...
                (sigma2_p + sigma2_po ./ nprime_o + sigma2_pi ./ nprime_i + ...
                 sigma2_poi_e ./ (nprime_o .* nprime_i));
        end

        function v = ceD(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, ...
                         sigma2_i, sigma2_oi, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CE", row "D-Coefficient".
            % Absolute error here adds sigma2_i and sigma2_oi but NOT sigma2_o: with
            % occasion fixed, the occasion main effect is not a source of error. This
            % omission is printed in the table.
            v = (sigma2_p + sigma2_po ./ nprime_o) ./ ...
                (sigma2_p + sigma2_po ./ nprime_o + sigma2_pi ./ nprime_i + ...
                 sigma2_poi_e ./ (nprime_o .* nprime_i) + ...
                 sigma2_i ./ nprime_i + sigma2_oi ./ (nprime_o .* nprime_i));
        end

        function v = ceIccRelative(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CE", row "ICC / Relative".
            v = (sigma2_p + sigma2_po) ./ (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e);
        end

        function v = ceIccAbsolute(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, sigma2_i, sigma2_oi)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CE", row "ICC / Absolute".
            v = (sigma2_p + sigma2_po) ./ ...
                (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e + sigma2_i + sigma2_oi);
        end

        % ---- CS: coefficient of stability (trial FIXED) -------------------------
        function v = csG(sigma2_p, sigma2_pi, sigma2_po, sigma2_poi_e, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CS", row "G-coefficient".
            % Trial fixed => sigma2_pi/nprime_i joins the UNIVERSE SCORE.
            v = (sigma2_p + sigma2_pi ./ nprime_i) ./ ...
                (sigma2_p + sigma2_po ./ nprime_o + sigma2_pi ./ nprime_i + ...
                 sigma2_poi_e ./ (nprime_o .* nprime_i));
        end

        function v = csD(sigma2_p, sigma2_pi, sigma2_po, sigma2_poi_e, ...
                         sigma2_o, sigma2_oi, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CS", row "D-Coefficient".
            % Absolute error adds sigma2_o and sigma2_oi but NOT sigma2_i -- the mirror
            % image of CE. Also printed in the table.
            v = (sigma2_p + sigma2_pi ./ nprime_i) ./ ...
                (sigma2_p + sigma2_po ./ nprime_o + sigma2_pi ./ nprime_i + ...
                 sigma2_poi_e ./ (nprime_o .* nprime_i) + ...
                 sigma2_o ./ nprime_o + sigma2_oi ./ (nprime_o .* nprime_i));
        end

        function v = csIccRelative(sigma2_p, sigma2_pi, sigma2_po, sigma2_poi_e)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CS", row "ICC / Relative".
            % TRANSCRIPTION NOTE: the last denominator component is printed
            % "sigma2_pio,e" in this row where sibling rows print "sigma2_poi,e". Read as
            % the same component; recorded here because the paper's own spelling varies.
            v = (sigma2_p + sigma2_pi) ./ (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e);
        end

        function v = csIccAbsolute(sigma2_p, sigma2_pi, sigma2_po, sigma2_poi_e, sigma2_o, sigma2_oi)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CS", row "ICC / Absolute".
            v = (sigma2_p + sigma2_pi) ./ ...
                (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e + sigma2_o + sigma2_oi);
        end

        % ---- CES: trial equivalence and stability (both RANDOM) -----------------
        function v = cesG(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CES", row "G-coefficient".
            % Nothing joins the universe score: it is bare sigma2_p.
            v = sigma2_p ./ (sigma2_p + ...
                PsyRATIndependentOracle.errVarRelativeTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, nprime_i, nprime_o));
        end

        function v = cesD(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, ...
                          sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CES", row "D-Coefficient".
            % The full absolute error: all three facet-only components present.
            v = sigma2_p ./ (sigma2_p + ...
                PsyRATIndependentOracle.errVarAbsoluteTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o));
        end

        function v = cesIccRelative(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CES", row "ICC / Relative".
            v = sigma2_p ./ (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e);
        end

        function v = cesIccAbsolute(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, ...
                                    sigma2_o, sigma2_i, sigma2_oi)
            % SOURCE: Rocha 2026 p.8 Table 3, block "CES", row "ICC / Absolute".
            v = sigma2_p ./ (sigma2_p + sigma2_po + sigma2_pi + sigma2_poi_e + ...
                             sigma2_o + sigma2_i + sigma2_oi);
        end

        function v = cutScoreTwoFacet(sigma2_p, sigma2_po, sigma2_pi, sigma2_poi_e, ...
                                      sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o, mu, C)
            % SOURCE: Rocha 2026 p.8 Table 3, row "Cut-Score".
            %
            % As in Table 2 this row is unpaired and its relative cell is blank: ONE
            % cut-score, absolute error only. The printed error term matches the CES
            % D-Coefficient's, i.e. the paper's cut score is the both-facets-random form.
            % PsyRAT additionally implements CE and CS cut-scores, which are a defensible
            % extension beyond the printed table rather than something it specifies.
            offsq = (mu - C) .^ 2;
            v = (sigma2_p + offsq) ./ (sigma2_p + offsq + ...
                PsyRATIndependentOracle.errVarAbsoluteTwoFacet(sigma2_pi, sigma2_po, ...
                    sigma2_poi_e, sigma2_o, sigma2_i, sigma2_oi, nprime_i, nprime_o));
        end

        %% ================================================================ TABLE 6
        %  Difference scores, X minus Y. Rocha 2026, p. 13.
        %
        %  WHAT THE TABLE PRINTS, AND WHAT IT DOES NOT
        %    Table 6 prints exactly four coefficients: a G and a D for the one-facet design
        %    (p x i), and a G and a D for the two-facet design (p x i x o). It prints NO ICC
        %    row, NO SEM row and NO cut-score row -- consistent with the cross-table
        %    cut-score census in cutScoreSchema below. The four coefficients here are
        %    therefore the entire anchorable content of the table.
        %
        %  THE ONE PATTERN EVERY BLOCK FOLLOWS
        %    A difference-score variance component is built from three printed pieces:
        %    measure X's own component over X's OWN replicate count, measure Y's over Y's
        %    own count, and MINUS TWICE the X-Y covariance over the HARMONIC MEAN of the
        %    two counts. That asymmetry -- own count for the variances, harmonic mean for
        %    the covariance -- is the whole substance of the table and is what diffBlock
        %    below encodes once.
        %
        %  WHY THE TWO-FACET ENTRIES ARE CES ONLY
        %    The table's Note states that its two-facet formulas "represent coefficients of
        %    equivalence and stability" (CES, both facets random), and then gives ONE
        %    instruction each for CE and CS: add the po block (CE) or the pi block (CS) to
        %    the NUMERATOR. It says nothing about two further moves that a CE or CS
        %    coefficient requires -- whether the block added to the numerator also LEAVES
        %    the denominator's error, and whether the fixed facet's own main effect leaves
        %    the ABSOLUTE error. Table 3 prints those omissions for the univariate case, so
        %    carrying them over to difference scores is a defensible analogy; but it is an
        %    analogy, not a transcription. Encoding it here would import a code-derived
        %    assumption into the one file whose value rests on having been written without
        %    reading the code. CE and CS difference coefficients are therefore deliberately
        %    NOT provided: they are an unanchored PsyRAT extension beyond the printed table.
        %
        %  DIFFERENCE-SCORE ICCs ARE LIKEWISE NOT PROVIDED
        %    PsyRAT reports difference-score ICCs as the n' = 1 form of these coefficients.
        %    That is a reasonable derivation by analogy with Table 2's undivided ICC row,
        %    but Table 6 has no ICC row to transcribe, so there is nothing here to anchor
        %    them against. They are labeled "derived, not Rocha Table 6" in production for
        %    exactly this reason.

        function v = covDivisor(nprime_X, nprime_Y)
            %COVDIVISOR  The divisor Table 6 puts under every covariance term.
            %
            % SOURCE: Rocha 2026 p.13 Table 6, Note: "For covariance terms, the harmonic
            % mean (n-with-double-dot) of the number of included observations for X and Y."
            % The table writes this symbol as n-double-dot-prime-i and n-double-dot-prime-o.
            %
            % COMPUTED HERE RATHER THAN CALLED FROM PRODUCTION, DELIBERATELY. PsyRAT has a
            % psyrat_harmmean helper that returns the same number. Calling it would make
            % this oracle agree with production by construction and collapse the only
            % independent check in the repository into a second mirror -- and it would do so
            % for the single quantity most likely to change (see below). The two-element
            % harmonic mean is one line; write the line.
            %
            % THIS HARDCODING ENCODES ROCHA'S NOTE AND DOES NOT PREJUDGE THE OPEN SCALER
            % QUESTION. Findings RC-01 / S13 ask whether the harmonic mean is the right
            % scaler at all, on the grounds that HM(n1,n2) implies the two condition means
            % share AM(n1,n2) trial levels while positional pairing lets them share only
            % min(n1,n2). That question is about what the TOOLBOX should do. This file's job
            % is only to say what the PAPER prints, and the paper prints the harmonic mean.
            % If the toolbox moves to a different scaler, the Table 6 tests are SUPPOSED to
            % fail: that failure is the record of a deliberate, adjudicated deviation from
            % the published formula. Do not silence it by editing this function.
            v = 2 ./ ((1 ./ nprime_X) + (1 ./ nprime_Y));
        end

        function v = diffBlock(sigma2_X, sigma2_Y, sigmaXY, divX, divY, divCov)
            %DIFFBLOCK  One printed block of a Table 6 denominator.
            %
            % SOURCE: Rocha 2026 p.13 Table 6, the repeated three-term group.
            %
            % Each measure's variance is divided by ITS OWN replicate count; the covariance
            % is doubled, subtracted, and divided by the harmonic-mean count. Callers pass
            % divCov explicitly rather than having it inferred, so that every divisor in
            % every coefficient below is visible at the point of composition.
            v = sigma2_X ./ divX + sigma2_Y ./ divY - 2 .* sigmaXY ./ divCov;
        end

        function v = diffUniverseScore(sigma2_pX, sigma2_pY, sigmaXY_p)
            %DIFFUNIVERSESCORE  The numerator shared by all four Table 6 coefficients.
            %
            % SOURCE: Rocha 2026 p.13 Table 6, the numerator of every printed formula.
            %
            % Undivided: no replicate count appears, because a universe score is defined
            % for the mean over the whole universe of admissible observations. Note that a
            % POSITIVE person-level covariance between X and Y SHRINKS this, which is the
            % algebraic form of the standard result that highly correlated constituents make
            % for an unreliable difference.
            v = sigma2_pX + sigma2_pY - 2 .* sigmaXY_p;
        end

        % ---- One-facet difference score, p x i -----------------------------------
        function v = diffErrVarRelativeOneFacet(sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                                                nprime_iX, nprime_iY)
            % SOURCE: Rocha 2026 p.13 Table 6, one-facet G-Coefficient denominator, the
            % terms beyond the universe score.
            % Only the person-carrying residual appears. The trial main effect is absent,
            % exactly as in Table 2's relative error.
            v = PsyRATIndependentOracle.diffBlock(sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                    nprime_iX, nprime_iY, ...
                    PsyRATIndependentOracle.covDivisor(nprime_iX, nprime_iY));
        end

        function v = diffErrVarAbsoluteOneFacet(sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                                                sigma2_iX, sigma2_iY, sigmaXY_i, ...
                                                nprime_iX, nprime_iY)
            % SOURCE: Rocha 2026 p.13 Table 6, one-facet D-Coefficient denominator, the
            % terms beyond the universe score.
            % Adds the trial main effect block -- the components with no person index --
            % divided by the same counts as the residual, not by 1.
            v = PsyRATIndependentOracle.diffErrVarRelativeOneFacet(sigma2_piX_e, ...
                    sigma2_piY_e, sigmaXY_pi_e, nprime_iX, nprime_iY) + ...
                PsyRATIndependentOracle.diffBlock(sigma2_iX, sigma2_iY, sigmaXY_i, ...
                    nprime_iX, nprime_iY, ...
                    PsyRATIndependentOracle.covDivisor(nprime_iX, nprime_iY));
        end

        function v = diffGOneFacet(sigma2_pX, sigma2_pY, sigmaXY_p, ...
                                   sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                                   nprime_iX, nprime_iY)
            % SOURCE: Rocha 2026 p.13 Table 6, block "One-facet Design (p x i)", row
            % "G-Coefficient".
            % GENERALIZABILITY: relative error. No trial main effect anywhere.
            uni = PsyRATIndependentOracle.diffUniverseScore(sigma2_pX, sigma2_pY, sigmaXY_p);
            v = uni ./ (uni + PsyRATIndependentOracle.diffErrVarRelativeOneFacet( ...
                    sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, nprime_iX, nprime_iY));
        end

        function v = diffDOneFacet(sigma2_pX, sigma2_pY, sigmaXY_p, ...
                                   sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                                   sigma2_iX, sigma2_iY, sigmaXY_i, ...
                                   nprime_iX, nprime_iY)
            % SOURCE: Rocha 2026 p.13 Table 6, block "One-facet Design (p x i)", row
            % "D-Coefficient".
            % DEPENDABILITY: absolute error. Same numerator as G. Unlike the univariate
            % Table 2 case, D <= G is NOT guaranteed here: the added trial block carries a
            % subtracted covariance, so a sufficiently large positive trial-level covariance
            % between X and Y can make that block negative and push D ABOVE G. The
            % possibility is a property of the printed formula, not a defect in it.
            uni = PsyRATIndependentOracle.diffUniverseScore(sigma2_pX, sigma2_pY, sigmaXY_p);
            v = uni ./ (uni + PsyRATIndependentOracle.diffErrVarAbsoluteOneFacet( ...
                    sigma2_piX_e, sigma2_piY_e, sigmaXY_pi_e, ...
                    sigma2_iX, sigma2_iY, sigmaXY_i, nprime_iX, nprime_iY));
        end

        % ---- Two-facet difference score, p x i x o, CES only ---------------------
        %
        %  These take STRUCTS rather than a positional list. The CES D-coefficient needs 21
        %  components and 4 replicate counts; as positional arguments that is 25 slots in
        %  which a transposition produces a wrong-but-plausible number silently. With named
        %  fields, a mistake is a "reference to non-existent field" error that names the
        %  field. For an oracle, failing loudly at the call site is worth more than matching
        %  the positional style of the one-facet methods above.
        %
        %  Required fields of c (all VARIANCES, except the sigmaXY_* which are COVARIANCES):
        %    sigma2_pX     sigma2_pY     sigmaXY_p        person
        %    sigma2_poX    sigma2_poY    sigmaXY_po       person x occasion
        %    sigma2_piX    sigma2_piY    sigmaXY_pi       person x trial
        %    sigma2_poiX_e sigma2_poiY_e sigmaXY_poi_e    three-way + undifferentiated error
        %    sigma2_oX     sigma2_oY     sigmaXY_o        occasion main effect      (D only)
        %    sigma2_iX     sigma2_iY     sigmaXY_i        trial main effect         (D only)
        %    sigma2_oiX    sigma2_oiY    sigmaXY_oi       occasion x trial          (D only)
        %  Required fields of n: iX, iY, oX, oY.

        function v = diffErrVarRelativeCesTwoFacet(c, n)
            % SOURCE: Rocha 2026 p.13 Table 6, two-facet G-Coefficient denominator, the
            % terms beyond the universe score.
            % The three person-carrying blocks. Note the three-way block's divisors: each
            % measure over the PRODUCT of its own two counts, the covariance over the
            % PRODUCT of the two harmonic means.
            hmi = PsyRATIndependentOracle.covDivisor(n.iX, n.iY);
            hmo = PsyRATIndependentOracle.covDivisor(n.oX, n.oY);
            v = PsyRATIndependentOracle.diffBlock(c.sigma2_poX, c.sigma2_poY, c.sigmaXY_po, ...
                    n.oX, n.oY, hmo) + ...
                PsyRATIndependentOracle.diffBlock(c.sigma2_piX, c.sigma2_piY, c.sigmaXY_pi, ...
                    n.iX, n.iY, hmi) + ...
                PsyRATIndependentOracle.diffBlock(c.sigma2_poiX_e, c.sigma2_poiY_e, ...
                    c.sigmaXY_poi_e, n.iX .* n.oX, n.iY .* n.oY, hmi .* hmo);
        end

        function v = diffErrVarAbsoluteCesTwoFacet(c, n)
            % SOURCE: Rocha 2026 p.13 Table 6, two-facet D-Coefficient denominator, the
            % terms beyond the universe score.
            % Adds the three blocks with no person index: occasion main, trial main, and
            % occasion x trial. Each mirrors the divisor pattern of its person-carrying
            % counterpart -- o like po, i like pi, oi like poi,e.
            hmi = PsyRATIndependentOracle.covDivisor(n.iX, n.iY);
            hmo = PsyRATIndependentOracle.covDivisor(n.oX, n.oY);
            v = PsyRATIndependentOracle.diffErrVarRelativeCesTwoFacet(c, n) + ...
                PsyRATIndependentOracle.diffBlock(c.sigma2_oX, c.sigma2_oY, c.sigmaXY_o, ...
                    n.oX, n.oY, hmo) + ...
                PsyRATIndependentOracle.diffBlock(c.sigma2_iX, c.sigma2_iY, c.sigmaXY_i, ...
                    n.iX, n.iY, hmi) + ...
                PsyRATIndependentOracle.diffBlock(c.sigma2_oiX, c.sigma2_oiY, c.sigmaXY_oi, ...
                    n.iX .* n.oX, n.iY .* n.oY, hmi .* hmo);
        end

        function v = diffCesGTwoFacet(c, n)
            % SOURCE: Rocha 2026 p.13 Table 6, block "Two-facet Design (p x i x o)", row
            % "G-Coefficient". CES: both facets random, so nothing joins the universe score.
            uni = PsyRATIndependentOracle.diffUniverseScore(c.sigma2_pX, c.sigma2_pY, ...
                    c.sigmaXY_p);
            v = uni ./ (uni + PsyRATIndependentOracle.diffErrVarRelativeCesTwoFacet(c, n));
        end

        function v = diffCesDTwoFacet(c, n)
            % SOURCE: Rocha 2026 p.13 Table 6, block "Two-facet Design (p x i x o)", row
            % "D-Coefficient". CES with the full absolute error.
            uni = PsyRATIndependentOracle.diffUniverseScore(c.sigma2_pX, c.sigma2_pY, ...
                    c.sigmaXY_p);
            v = uni ./ (uni + PsyRATIndependentOracle.diffErrVarAbsoluteCesTwoFacet(c, n));
        end

        %% ============================================== CROSS-TABLE CUT-SCORE SCHEMA
        function v = cutScoreSchema(sigma2_p, absErrVar, mu, C)
            %CUTSCORESCHEMA  The single form every printed cut score in Rocha 2026
            %instantiates.
            %
            % SOURCE: Rocha 2026, Tables 2 (p.7), 3 (p.8), 4 (p.8) and 5 (p.9). A
            % table-by-table census found exactly ONE cut-score schema, instantiated four
            % times, one per univariate design, with ZERO relative variants: in all four
            % tables the relative cell is blank and the row carries no sub-label. Tables 6
            % and 7 (difference scores) print NO cut-score row at all, consistent with the
            % paper's own cross-reference on p.7 to "Tables 2-5".
            %
            % The schema is: take the design's own ABSOLUTE error variance -- term for
            % term and divisor for divisor, exactly as that table prints it in its "Error
            % Variances / Absolute" cell -- and add the squared cut-score deviation to
            % both numerator and denominator.
            %
            % Stating it once, in this form, is what makes the invariant testable for any
            % design PsyRAT supports, including designs the paper does not tabulate.
            offsq = (mu - C) .^ 2;
            v = (sigma2_p + offsq) ./ (sigma2_p + offsq + absErrVar);
        end

    end
end
