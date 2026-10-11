function cleanup = psyrat_test_fitlme_injection(warnAt,failAt,valueAt,value)
%Test-only real-fit forwarding with deterministic warning/error/estimate faults.
%Temporary shims never reside on a production path after cleanup.
if nargin<3; valueAt=0; value=0; end
d=tempname; mkdir(d);
localWrite(fullfile(d,'fitlme.m'),{
    'function out=fitlme(varargin)'
    's=getappdata(0,''psyrat_test_fitlme_injection''); s.calls=s.calls+1; s.injectCov=false;'
    'setappdata(0,''psyrat_test_fitlme_injection'',s);'
    'if any(s.calls==s.failAt); error(''PsyRATTest:InjectedRefitFailure'',''Injected refit failure.''); end'
    'd=fileparts(mfilename(''fullpath'')); rmpath(d); c=onCleanup(@()addpath(d,''-begin''));'
    'out=fitlme(varargin{:});'
    'if any(s.calls==s.valueAt); out=PsyRATTestLmeProxy(out,s.value); end'
    'if any(s.calls==s.warnAt); warning(''PsyRATTest:OptimizerNonConvergence'',''Injected optimizer warning.''); end'
    'end'});
setappdata(0,'psyrat_test_fitlme_injection',struct('calls',0,'warnAt',warnAt, ...
    'failAt',failAt,'valueAt',valueAt,'value',value,'injectCov',false));
addpath(d,'-begin');
cleanup=onCleanup(@()localCleanup(d));
end
function localWrite(p,lines)
f=fopen(p,'w'); c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',lines{:});
end
function localCleanup(d)
rmpath(d); rmappdata(0,'psyrat_test_fitlme_injection'); rmdir(d,'s');
end
