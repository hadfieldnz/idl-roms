;+
; NAME:
;   MGH_ROMS_SUMMARIZE_BENCHMARK
;
; PURPOSE:
;   Read one or more ROMS log files and print out info on elapsed time.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_SUMMARIZE_BENCHMARK, file_pattern
;
; POSITIONAL PARAMETERS:
;   file_pattern (input, scalar string)
;     Pattern for input file names.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2007-05:
;     Written.
;-

pro mgh_roms_summarize_benchmark, file_pattern

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(file_pattern) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file_pattern'

   file = file_search(file_pattern, COUNT=n_file)

   if n_file eq 0 then $
        message, 'No matching files found'

   for f=0,n_file-1 do begin

      info = mgh_roms_read_log(file[f])

      print, file_basename(file[f]), round((info.date1-info.date0)*24.0D*3600.0D), FORMAT='(%"%s:  %0.0d")'

   endfor

end

