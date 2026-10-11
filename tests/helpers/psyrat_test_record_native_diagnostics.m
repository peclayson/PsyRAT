function psyrat_test_record_native_diagnostics(testCase,rel)
%Live integration check and optional evidence capture, not an adequacy gate.
%Existing recovery targets/tolerances stay in their original test classes.
tables=rel.out.conv.data;
if ~iscell(tables{1}); tables={tables}; end
for k=1:numel(tables)
    c=tables{k}; m=psyrat_diagnostic_metadata(c);
    testCase.verifyEqual(m.format,'native_rank_fold');
    testCase.verifyFalse(m.degraded);
    testCase.verifyFalse(any(startsWith(c(:,1),'required_shape.')));
    testCase.verifyFalse(any(strcmp(c(:,1),'native_output_contract_unavailable')));
    testCase.verifyTrue(all(ismember(m.expected_names,c(2:end,1))));
end
rel=psyrat_checkconv(rel);
rel.out.conv.ndivergent=psyrat_total_divergences(rel.out.conv.data);
folder=getenv('PSYRAT_NATIVE_VALIDATION_DIR');
if ~isempty(folder)
    if ~isfolder(folder); mkdir(folder); end
    [~,revision]=system('git rev-parse HEAD');
    evidence_scope='Existing Gaussian point-recovery/output fixture; no MCSE or coverage certificate';
    if isa(testCase,'TestDodNativeHmcRecovery')
        evidence_scope='Pipeline regression only: existing generator has person-nested trial effects; see NATIVE_RECOVERY_LANE_REVIEW.md';
    end
    save(fullfile(folder,[class(testCase) '.mat']),'rel','revision','evidence_scope','-v7');
end
end
