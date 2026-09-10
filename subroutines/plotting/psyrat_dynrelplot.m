function fh = psyrat_dynrelplot(summary,varargin)
%Plot the dynamic/conditional reliability surface (Rast & Clayson, analysis 11).
%
% fh = psyrat_dynrelplot(summary,'gcoeff',2,'cutoff',.70)
%
%For one dimension the reliability G(z)/D(z) is drawn as a curve over the
%standardized predictor with a shaded credible band, one subplot per
%group/event stratum. For two dimensions the point-estimate surface is drawn as
%a filled contour over the two standardized predictors, one subplot per stratum.
%
%For the subject-level variants (any stratum carrying a per-participant
%ssrel_table: the dynamic difference designs 13/19/21, their two-facet siblings,
%the subject-level non-difference designs 26/27, and the DoD designs 28/29 in
%either family) the curve/surface is the typical-person population reference,
%and each participant's own phi at their own standardized dimension value is
%overlaid as a point.
%
%Inputs
% summary - struct from psyrat_dynrel_summary.
%
%Optional Inputs
% gcoeff - 1 = dependability D(z) (absolute), 2 = generalizability G(z)
%          (relative; default, matching Rast's dynamic G_t).
% cutoff - optional reference reliability value drawn as a horizontal line
%          (one dimension) or a highlighted contour level (two dimensions).
%          Default [] (none).
%
%Output
% fh - figure handle (Tag 'psyrat_output').

% Copyright (C) 2016-2026 Peter E. Clayson
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

gcoeff = psyrat_opt(varargin,'gcoeff',2);
cutoff = psyrat_opt(varargin,'cutoff',[]);
if ~any(gcoeff == [1 2])
    error('varargin:gcoeff','WARNING: gcoeff is invalid. Valid values are 1 or 2\n');
end

if gcoeff == 1
    coefffield = 'D'; ylbl = 'Dependability'; pref = 'Dynamic Dependability';
else
    coefffield = 'G'; ylbl = 'Generalizability'; pref = 'Dynamic Generalizability';
end

ndim = summary.ndim;
nstrata = numel(summary.strata);
fsize = 16;

%subplot grid
xplots = ceil(sqrt(nstrata));
yplots = ceil(nstrata / xplots);

fh = figure;
fh.Tag = 'psyrat_output';
set(gcf,'NumberTitle','Off');
fh.Position = [125 500 950 520];
set(gcf,'Name',sprintf('%s as a function of the dimension(s)',pref));

dim1lab = sprintf('%s (standardized)',char(string(summary.dim_names{1})));

for s = 1:nstrata
    st = summary.strata(s);
    coeff = st.(coefffield);
    subplot(yplots,xplots,s);

    if ndim == 1
        z = st.z1(:)';
        pt = coeff.pt(:)'; ll = coeff.ll(:)'; ul = coeff.ul(:)';
        %shaded credible band
        patch([z fliplr(z)],[ll fliplr(ul)],[0.3 0.5 0.9],...
            'FaceAlpha',0.20,'EdgeColor','none'); hold on;
        plot(z,pt,'Color',[0.10 0.25 0.65],'LineWidth',2);
        ylim([0 1]); xlim([min(z) max(z)]);
        set(gca,'fontsize',fsize);
        %the dimension name is user text (a column header), so it is printed
        %literally like the event and group names (Interpreter none)
        xlabel(dim1lab,'FontSize',fsize,'Interpreter','none');
        ylabel(ylbl,'FontSize',fsize);
        if ~isempty(cutoff)
            psyrat_addhline(cutoff,'Color','b','LineStyle',':');
        end
        %subject-level overlay: each participant's phi at their own z
        if isfield(st,'ssrel_table') && ~isempty(st.ssrel_table)
            tab = st.ssrel_table;
            if gcoeff == 1; sp = tab.dep_pt; else; sp = tab.gen_pt; end
            scatter(tab.z1,sp,36,[0.85 0.33 0.10],'filled',...
                'MarkerEdgeColor','k','MarkerFaceAlpha',0.70);
        end
    else
        z1 = st.z1(:)'; z2 = st.z2(:)';
        %point-estimate surface; rows index z1, cols index z2 -> transpose so
        %x = z1, y = z2 in the image.
        contourf(z1,z2,coeff.pt',12,'LineColor','none'); hold on;
        clim([0 1]); colormap(parula);
        cb = colorbar; cb.Label.String = ylbl; cb.Label.FontSize = fsize-2;
        set(gca,'fontsize',fsize);
        xlabel(dim1lab,'FontSize',fsize,'Interpreter','none');
        ylabel(sprintf('%s (standardized)',char(string(summary.dim_names{2}))),...
            'FontSize',fsize,'Interpreter','none');
        if ~isempty(cutoff)
            contour(z1,z2,coeff.pt',[cutoff cutoff],'Color','w',...
                'LineWidth',2,'LineStyle','--');
        end
        %subject-level overlay: each participant at their own (z1,z2), colored
        %by phi on the same scale as the reference surface.
        if isfield(st,'ssrel_table') && ~isempty(st.ssrel_table)
            tab = st.ssrel_table;
            if gcoeff == 1; sp = tab.dep_pt; else; sp = tab.gen_pt; end
            scatter(tab.z1,tab.z2,40,sp,'filled','MarkerEdgeColor','k');
        end
    end

    %n' annotation; the trial + occasion two-facet variant also carries an
    %occasion n' (st.nocc), shown alongside the trial n'.
    %The gamma difference variant carries PER-EVENT n' (st.obs_ev = [n1 n2]) and
    %computes the surface at those exact counts, so show both when they differ -
    %the scalar st.obs is their mean and would name a count neither event has.
    if isfield(st,'obs_ev') && numel(st.obs_ev) == 4
        %DoD variants: all four per-cell counts (their mean rarely names a
        %count any cell has)
        nstr = sprintf('n''=%d/%d/%d/%d',st.obs_ev(1),st.obs_ev(2),...
            st.obs_ev(3),st.obs_ev(4));
        if isfield(st,'nocc') && ~isempty(st.nocc)
            nstr = sprintf('%s, occ=%d',nstr,st.nocc);
        end
    elseif isfield(st,'obs_ev') && numel(st.obs_ev) == 2 && st.obs_ev(1) ~= st.obs_ev(2)
        nstr = sprintf('n''=%d/%d',st.obs_ev(1),st.obs_ev(2));
    elseif isfield(st,'nocc') && ~isempty(st.nocc)
        nstr = sprintf('n''=%d, occ=%d',st.obs,st.nocc);
    else
        nstr = sprintf('n''=%d',st.obs);
    end
    ttl = st.label;
    if ~strcmpi(ttl,'none')
        ttl = strrep(ttl,'_;_',' / ');
        title(sprintf('%s (%s)',ttl,nstr),'FontSize',fsize+2,...
            'Interpreter','none');
    else
        title(nstr,'FontSize',fsize+2);
    end
    hold off;
end

%offer a one-click figure export (PNG/PDF) for parity with the standard viewers.
%The button is hidden during export so it is not captured in the saved image.
uicontrol(fh,'Style','push','String','Save figure (PNG/PDF)',...
    'Units','pixels','Position',[10 10 160 28],...
    'Callback',@(src,~) local_savefig(fh,src));

end

function local_savefig(fh,btn)
%Export the reliability-surface figure to PNG or PDF. exportgraphics captures the
%full figure (every stratum subplot); the export button is hidden during the
%write so it is not included in the saved image, and restored afterward (even on
%error) via onCleanup.
[fn,pth] = uiputfile({'*.png','PNG image (*.png)';'*.pdf','PDF (*.pdf)'},...
    'Save reliability surface figure','dynrel_reliability_surface.png');
if isequal(fn,0); return; end
outfile = fullfile(pth,fn);
vis = get(btn,'Visible');
set(btn,'Visible','off');
drawnow;
restorebtn = onCleanup(@() set(btn,'Visible',vis));
exportgraphics(fh,outfile,'Resolution',150);
fprintf('Saved reliability surface figure to %s\n',outfile);
end
