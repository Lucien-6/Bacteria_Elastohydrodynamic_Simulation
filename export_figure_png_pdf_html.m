function export_figure_png_pdf_html(fig, pngPath, pdfPath, htmlPath)
%EXPORT_FIGURE_PNG_PDF_HTML  Save figure as PNG (600 dpi), PDF (vector), optional HTML.
%
%   export_figure_png_pdf_html(fig, pngPath, pdfPath, htmlPath)
%
%   Interactive HTML (pan / zoom / rotate in a web canvas) is written only
%   when htmlPath is nonempty and MATLAB R2026a or later; see exportgraphics.
%   Pass [] or '' for htmlPath to skip HTML entirely.
%
%  Author: Lucien (lucien-6@qq.com)
%  Created: 2026-04-22

exportgraphics(fig, pngPath, 'Resolution', 600);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

if nargin < 4 || isempty(htmlPath)
    return
end

if ~local_matlab_release_at_least_r2026a()
    return
end

try
    exportgraphics(fig, htmlPath);
catch ME
    warning('export_figure_png_pdf_html:htmlFailed', ...
        'Interactive HTML export failed: %s', ME.message);
end
end


function tf = local_matlab_release_at_least_r2026a()
tf = false;
try
    tf = ~isMATLABReleaseOlderThan('R2026a');
catch
end
end
