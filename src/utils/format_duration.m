function str = format_duration(seconds_in, style)
%FORMAT_DURATION  Convert duration in seconds to a human-readable string
%
%  Description:
%    Unified duration formatter shared across ODE progress bars, simulation
%    logs, visualisation loops, and parallel batch drivers.  Rounds before
%    splitting into days/hours/minutes/seconds so rollover edge cases such
%    as '1 min 60 s' or '0:60' can never occur.
%
%  Syntax:
%    str = format_duration(seconds_in)
%    str = format_duration(seconds_in, style)
%
%  Inputs:
%    seconds_in - Non-negative duration in seconds (scalar numeric).
%                 Negative, NaN, or Inf inputs return '--'.
%    style      - (Optional) Output style:
%                   'long'    (default) '3.7 s', '2 min 5.3 s',
%                             '1 h 2 min 5.1 s', '1 d 3 h 15 min'
%                   'compact' Progress-bar friendly: '3.7 s', '2 min 05 s',
%                             '1 h 02 min', '1 d 03 h'
%                   'hms'     Clock style: '0:04', '2:05', '1:02:05',
%                             '1d 03:15:00'
%
%  Outputs:
%    str - Formatted char row vector
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-11-24

if nargin < 2 || isempty(style)
    style = 'long';
end

if isempty(seconds_in) || ~isscalar(seconds_in) || ~isfinite(seconds_in) || seconds_in < 0
    str = '--';
    return;
end
if seconds_in >= 1000 * 86400
    str = '>999 d';
    return;
end

switch lower(char(style))
    case 'long'
        % Round to nearest 0.1 s first so 59.96 s rolls over cleanly to 1 min 0.0 s.
        tenths = round(seconds_in * 10);
        if tenths < 600
            str = sprintf('%.1f s', tenths / 10);
        elseif tenths < 36000
            m = floor(tenths / 600);
            s = (tenths - m * 600) / 10;
            str = sprintf('%d min %.1f s', m, s);
        elseif tenths < 864000
            h = floor(tenths / 36000);
            rem_t = tenths - h * 36000;
            m = floor(rem_t / 600);
            s = (rem_t - m * 600) / 10;
            str = sprintf('%d h %d min %.1f s', h, m, s);
        else
            total_min = round(seconds_in / 60);
            d = floor(total_min / 1440);
            rem_m = total_min - d * 1440;
            h = floor(rem_m / 60);
            m = rem_m - h * 60;
            str = sprintf('%d d %d h %d min', d, h, m);
        end

    case 'compact'
        if seconds_in < 59.95
            str = sprintf('%.1f s', seconds_in);
            return;
        end
        total_s = round(seconds_in);
        if total_s < 3600
            m = floor(total_s / 60);
            s = total_s - m * 60;
            str = sprintf('%d min %02d s', m, s);
        elseif total_s < 86400
            total_min = round(total_s / 60);
            if total_min >= 1440
                d = floor(total_min / 1440);
                h = floor((total_min - d * 1440) / 60);
                str = sprintf('%d d %02d h', d, h);
            else
                h = floor(total_min / 60);
                m = total_min - h * 60;
                str = sprintf('%d h %02d min', h, m);
            end
        else
            total_h = round(seconds_in / 3600);
            d = floor(total_h / 24);
            h = total_h - d * 24;
            str = sprintf('%d d %02d h', d, h);
        end

    case 'hms'
        total_s = round(seconds_in);
        if total_s < 3600
            m = floor(total_s / 60);
            s = total_s - m * 60;
            str = sprintf('%d:%02d', m, s);
        elseif total_s < 86400
            h = floor(total_s / 3600);
            rem_s = total_s - h * 3600;
            m = floor(rem_s / 60);
            s = rem_s - m * 60;
            str = sprintf('%d:%02d:%02d', h, m, s);
        else
            d = floor(total_s / 86400);
            rem_s = total_s - d * 86400;
            h = floor(rem_s / 3600);
            rem_s = rem_s - h * 3600;
            m = floor(rem_s / 60);
            s = rem_s - m * 60;
            str = sprintf('%dd %02d:%02d:%02d', d, h, m, s);
        end

    otherwise
        error('format_duration:UnknownStyle', ...
            'Unknown style ''%s''. Expected ''long'', ''compact'', or ''hms''.', char(style));
end

end
