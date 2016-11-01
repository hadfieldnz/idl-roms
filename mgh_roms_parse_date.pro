; svn $Id$
;+
; NAME:
;   MGH_ROMS_PARSE_DATE
;
; PURPOSE:
;   Parse a date-time string as written to a ROMS log file and
;   return the result as a Julian date
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   result = MGH_ROMS_PARSE_DATE(dt)
;
; POSITIONAL PARAMETERS:
;   dt (input, string scalar or array)
;     One or more strings of the form "Wednesday - June 1, 2005 - 11:58:54 AM"
;
; RETURN VALUE:
;   The function returns a double-precision floating-point result with
;   the same shape as the input
;
; PROCEDURE:
;   Break up the string, reassemble it and interpret the result using
;   IDL's calendar formats. It's quite clever, really, and fairly
;   robust.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2006-03:
;     Written.
;-

function mgh_roms_parse_date, dt

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(dt) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'dt'

   result = mgh_reproduce(!values.d_nan, dt)

   for i=0,n_elements(dt)-1 do begin

      tt = !values.d_nan

      s = strsplit(dt[i], '[ -,:]', /EXTRACT)

      if n_elements(s) ge 8 then begin
         fmt = '(I4.4,1X,A3,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,A2)'
         ; Work-around for bug in handling of times with AM/PM strings:
         if fix(s[4]) eq 12 then s[4] = '00'
         ss = string(s[3], s[1], s[2], s[4], s[5], s[6], s[7], FORMAT=temporary(fmt))
         fmt = '(C(CYI4.4,1X,CMoA3,1X,CDI2.2,1X,CHI2.2,1X,CMI2.2,1X,CSI2.2,1X,CAPA))'
         reads, ss, tt, FORMAT=temporary(fmt)
      end

      result[i] = tt

   endfor

   return, result

end

