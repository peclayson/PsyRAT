classdef TestCase24AvailableOccasions < PsyRATTestBase
    % Frozen independent effect-incidence reference; no production oracle.
    methods (Test)
        function balancedControl(t)
            localCheck(t,'balanced');
        end
        function unequalPersonCounts(t)
            localCheck(t,'audit_irregular');
        end
        function missingOccasion(t)
            localCheck(t,'audit_missing_occasion');
        end
        function equalOccasionMeansNotPooledRows(t)
            localCheck(t,'within_person_unequal');
        end
        function separateStrataAndComponents(t)
            localCheck(t,'multiple_group_event_strata');
            localCheck(t,'multiple_group_event_distinct_components');
        end
        function collisionSafePersonOccasion(t)
            localCheck(t,'id_occasion_delimiter_collision');
        end
        function rejectAmbiguousStratumDisplayMapping(t)
            [rel,~] = localFixture('group_event_display_delimiter_collision');
            t.verifyError(@() psyrat_splits_summary(rel),...
                'psyrat_splits_summary:ambiguousStratum');
        end
        function literalNoneGroupControl(t)
            localCheck(t,'literal_none_is_real_group_label');
        end
        function originalCellRowOrderAndPhysicalLabels(t)
            localCheck(t,'physical_labels_shared_ordinal_model');
            localCheck(t,'physical_labels_disjoint_ordinal_model');
            localCheck(t,'label_renaming_and_noncontiguous_rows');
        end
        function numericOverridesKeepBalancedProjection(t)
            names = {'override_obs_only','override_nocc_only','override_both',...
                'override_fractional_rounding','explicit_same_count_still_projection'};
            for k=1:numel(names); localCheck(t,names{k}); end
        end
        function observedKeywordControl(t)
            [rel,~] = localFixture('audit_irregular');
            a=psyrat_splits_summary(rel);
            b=psyrat_splits_summary(rel,'nocc','observed');
            t.verifyEqual(a,b);
        end
        function averageVariancesBeforeRatiosAndPosteriorSummaries(t)
            [rel,f] = localFixture('audit_missing_occasion');
            % Two draws change only P, so the frozen U/error at P=4 fully
            % determine independent targets without a coefficient oracle.
            fields=fieldnames(rel.out);
            for k=1:numel(fields)
                if startsWith(fields{k},'sig_')
                    rel.out.(fields{k})=repmat(rel.out.(fields{k}),2,1);
                end
            end
            rel.out.sig_id=[2;3];
            s=psyrat_splits_summary(rel,'CI',.90);
            for c=1:3
                expected=f.strata.expected(c,:);
                absvar=expected(3)^2; relvar=expected(4)^2;
                u0=expected(1)*absvar/(1-expected(1));
                u=[u0;u0+5];
                d=u./(u+absvar); g=u./(u+relvar);
                t.verifyEqual(s.strata.coef(c).D,[min(d) mean(d) max(d)],'AbsTol',2e-12);
                t.verifyEqual(s.strata.coef(c).G,[min(g) mean(g) max(g)],'AbsTol',2e-12);
            end
        end
        function projectionMetadataRemainsExplicit(t)
            [rel,f]=localFixture('override_nocc_only');
            s=psyrat_splits_summary(rel,'nocc',f.nocc);
            d=s.strata.design;
            t.verifyEqual(d.mode,'balanced_projection');
            t.verifySubstring(d.fixed_facet_target,'common balanced reference support');
            t.verifyEqual(d.projection_nsplit,f.strata.nsplit);
            t.verifyEqual(d.projection_nocc,f.strata.nocc);
            t.verifyTrue(d.irregular);
        end
        function auditableWeightsAndCounts(t)
            [rel,~] = localFixture('within_person_unequal');
            s=psyrat_splits_summary(rel); d=s.strata(1).design;
            t.verifyEqual(d.mode,'observed_available_occasions');
            t.verifyEqual(d.row_weights,repmat([1/2;1/6;1/6;1/6],2,1),'AbsTol',1e-14);
            t.verifyEqual(d.row_split_index,repmat([1;1;2;3],2,1));
            t.verifyEqual(d.mean_factors,[1/2 1/2 1/3 79/288],'AbsTol',1e-14);
            t.verifyTrue(d.irregular);
            t.verifyFalse(d.missing_occasions);
            t.verifyTrue(contains(d.score_rule,'available occasion means'));
            [rel,~] = localFixture('audit_missing_occasion');
            s=psyrat_splits_summary(rel);
            t.verifyTrue(s.strata(1).design.missing_occasions);
        end
        function case24PresentationIsConditionalAndAuditable(t)
            [rel,~] = localFixture('audit_missing_occasion');
            s=psyrat_splits_summary(rel);
            t.verifyFalse(contains(s.gbias_note,'downward-biased'));
            t.verifyTrue(contains(s.gbias_note,'available occasion'));
            [rows,cols]=psyrat_splits_coef_rows(s,psyrat_unitlabels(3));
            t.verifyTrue(any(contains(cols,'reference')));
            t.verifyFalse(any(contains(cols,'lower bound')));
            t.verifyTrue(any(contains(string(rows(1,:)),'available occasion')));
            rowDesign=jsondecode(rows{1,strcmp(cols,'Design metadata')});
            t.verifyTrue(rowDesign.missing_occasions);
            t.verifyEqual(rowDesign.row_weights,s.strata.design.row_weights);
            vc=psyrat_splits_varcomp_concat(s);
            t.verifyTrue(ismember('design_metadata',vc.Properties.VariableNames));
            d=jsondecode(vc.design_metadata(1));
            t.verifyTrue(d.missing_occasions);
        end
        function case24PreflightNamesAvailableOccasionTarget(t)
            [rel,~]=localFixture('within_person_unequal');
            report=psyrat_preflight_validate('datatable',rel.data,'splits',3);
            ids={report.warnings.id};
            t.assertTrue(any(strcmp(ids,'preflight:splitsUnequalCounts')));
            w=report.warnings(strcmp(ids,'preflight:splitsUnequalCounts'));
            text=jsonencode(w);
            t.verifyTrue(contains(text,'available occasion means'));
            t.verifyFalse(contains(text,'harmonic mean of the counts'));
        end
        function uniqueDelimiterContainingStratumRemainsValid(t)
            [rel,~]=localFixture('balanced');
            rel.data.group=repmat({'A_;_B'},height(rel.data),1);
            rel.data.event=repmat({'C'},height(rel.data),1);
            rel.groups={'A_;_B'}; rel.events={'C'};
            rel.out.labels={'A_;_B_;_C'};
            s=psyrat_splits_summary(rel);
            t.verifyEqual(s.strata.design.group,'A_;_B');
            t.verifyEqual(s.strata.design.event,'C');
            t.verifyEqual(s.strata.coef(3).D(2),25600/48477,'AbsTol',2e-12);
        end
    end
end

function localCheck(t,name)
[rel,f]=localFixture(name);
s=psyrat_splits_summary(rel,'obs',f.obs,'nocc',f.nocc);
t.verifyNumElements(s.strata,numel(f.strata));
for k=1:numel(f.strata)
    st=s.strata(k); ref=f.strata(k);
    for c=1:3
        actual=[st.coef(c).D;st.coef(c).G;st.coef(c).SEM_abs;st.coef(c).SEM_rel];
        t.verifyEqual(actual,repmat(ref.expected(c,:)',1,3),'AbsTol',2e-12,...
            sprintf('%s stratum %d coefficient %d',name,k,c));
    end
end
end

function [rel,f]=localFixture(name)
root=fileparts(mfilename('fullpath'));
all=jsondecode(fileread(fullfile(root,'fixtures','case24_incidence_reference.json')));
f=all.fixtures(strcmp({all.fixtures.name},name));
rel.analysis='trt_splits'; rel.splits=3;
rel.data=struct2table(f.rows); rel.data.meas=(1:height(rel.data))';
rel.groups=unique({f.rows.group}); rel.events=unique({f.rows.event});
if isequal(rel.groups,{'none'}); rel.groups='none'; end
if isequal(rel.events,{'none'}); rel.events='none'; end
fld={'sig_id','sig_occ','sig_trl','sig_trlxid','sig_occxid',...
    'sig_trlxocc','sig_posxid','sig_err'};
rel.out.labels=cell(1,numel(f.strata));
for k=1:numel(f.strata)
    st=f.strata(k);
    if ischar(rel.groups) && ischar(rel.events); label='none';
    elseif ischar(rel.groups); label=st.event;
    elseif ischar(rel.events); label=st.group;
    else; label=[st.group '_;_' st.event]; end
    rel.out.labels{k}=label;
    for j=1:numel(fld); rel.out.(fld{j})(1,k)=st.sd(j); end
end
end
