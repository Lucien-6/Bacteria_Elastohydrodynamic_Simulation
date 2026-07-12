function export_figure_png_pdf_html(fig, pngPath, pdfPath, htmlPath)
%EXPORT_FIGURE_PNG_PDF_HTML  Save figure as PNG (600 dpi), PDF (vector), optional HTML.
%
%  Description:
%    Exports a MATLAB figure to three formats simultaneously:
%    PNG at 600 dpi resolution, PDF with vector graphics, and optionally
%    an interactive HTML file (pan/zoom/rotate in a web canvas) when
%    MATLAB R2026a or later is available.  The function has no return
%    value and writes directly to the supplied file paths.  PNG and PDF
%    export are always attempted.  HTML export is skipped when htmlPath
%    is empty or the MATLAB release predates R2026a; if a supported HTML
%    export is attempted but fails, a MATLAB warning is emitted.
%
%  Syntax:
%    export_figure_png_pdf_html(fig, pngPath, pdfPath, htmlPath)
%
%  Inputs:
%    fig      - MATLAB figure handle to export
%    pngPath  - Output file path for PNG image
%    pdfPath  - Output file path for PDF document
%    htmlPath - Output file path for interactive HTML, or [] / '' to skip
%
%  Outputs:
%    None.
%
%  Side Effects:
%    Writes PNG and PDF files, may write an HTML file, and may emit
%    warning export_figure_png_pdf_html:htmlFailed on HTML export failure.
%
%  Author: Lucien (lucien-6@qq.com)
%  Created: 2026-04-22
%  Updated: 2026-07-12

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
