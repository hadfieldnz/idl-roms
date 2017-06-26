;+
; NAME:
;   MGH_ROMS_PROCESS_LOG
;
; PURPOSE:
;   Read info from a string vector representing the contents of a ROMS log
;   (standard output) file.
;
; CALLING SEQUENCE:
;   result = mgh_roms_process_log(text)
;
; POSITIONAL PARAMETERS:
;   text (input, string vector)
;     The log file contents
;
; RETURN VALUE:
;   The function returns a structure with a whole load of different tags.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2015-04:
;     Written.
;   Mark Hadfield, 2017-05:
;     Added support for the new string-based date output (ticket 724).
;   Mark Hadfield, 2017-06:
;     Fixed up code for reading float source locations (process_src option).
;     Lines are now processed in a WHILE loop instead of a FOR loop, allowing
;     two lines to be process together.
;-
function mgh_roms_process_log, text

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(text) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'text'

   if size(text, /TYPE) ne 7 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'text'

   prm_info = {dstart: 0D, time_ref: 0D}
   src_info = {g: 0, c: 0, n: 0, ft0: 0D, fx0: 0D, fy0: 0D, fz0: 0D, fdt: 0D, fdx: 0D, fdy: 0D, fdz: 0D}
   step_info = {grid: 0, step: 0, time: 0D, ke: 0D, pe: 0D, te: 0D, nv: 0D}
   cfl_info = {cu: 0D, cv: 0D, cw: 0D, smax: 0D}

   src_list = list()
   step_list = list()
   cfl_list = list()

   process_prm = !false
   process_src = !false
   process_step = !false

   n_line = n_elements(text)

   result = {n_line: n_line, version: '', date0: !values.d_nan, date1: !values.d_nan, dt: !values.d_nan, n_flt: 0, prm: prm_info}

   i_line = 0

   while i_line lt n_line do begin

      line = text[i_line]

      ;; Search for various lines with summary info

      if strmatch(line, string(replicate(32B, 20))+'*-*-* ?M') then begin
         result.date0 = mgh_roms_parse_date(strmid(line, 20))
         i_line++
         continue
      endif

      if strmatch(line, '*ROMS/TOMS version*') then begin
         result.version = strtrim(strmid(line, 44), 2)
         i_line++
         continue
      endif

      if strmatch(line, '*ROMS/TOMS: DONE*') then begin
         result.date1 = mgh_roms_parse_date(strmid(line, 20))
         i_line++
         continue
      endif

      ;; Search for lines indicating what sort of info we can expect

      if strmatch(line, '*Physical Parameters, Grid:*', /FOLD_CASE) then begin
         process_prm = !true
         process_src = !false
         process_step = !false
         i_line++
         continue
      endif

      if strmatch(line, '*Floats Initial Locations*', /FOLD_CASE) then begin
         process_prm = !false
         process_src = !true
         process_step = !false
         i_line++
         continue
      endif

      if process_src && strmatch(line, '*Nfloats*Number of float trajectories to compute*') then begin
         s = strsplit(line, /EXTRACT)
         result.n_flt = long(s[0])
         process_src = !false
         i_line++
         continue
      endif

      if strmatch(line, '*STEP*KINETIC*POTEN*TOTAL*NET*', /FOLD_CASE) then begin
         process_prm = !false
         process_src = !false
         process_step = !true
         i_line++
         continue
      endif

      ;; Process physical parameters info

      if process_prm then begin

         ;; Split the line into substrings at spans of white space.
         s = strsplit(line, /EXTRACT)

         ;; If there are at least 2 elements and the first is numeric, then
         ;; look for numeric data
         if n_elements(s) ge 2 && mgh_str_isnumber(s[0]) then begin
            if strmatch(s[1], 'dstart', /FOLD_CASE) then $
               result.prm.dstart = double(s[0])
            if strmatch(s[1], 'time_ref', /FOLD_CASE) then begin
               tt = long(s[0])
               result.prm.time_ref = mgh_dt_julday(YEAR=tt/10000, MONTH=(tt/100) mod 100, DAY=tt mod 100)
            endif
            i_line++
            continue
         endif

      endif

      ;; Process float-source info

      if process_src then begin

         ;; Split the line into substrings at spans of white space.
         s = strsplit(line, /EXTRACT)

         ;; If there are 11 elements, all numeric, then we have
         ;; single-line float-source data
         if n_elements(s) eq 11 && min(mgh_str_isnumber(s)) gt 0 then begin
            reads, line, src_info
            src_list->Add, src_info
            i_line++
            continue
         endif

         ;; If there are 7 elements, all numeric, then we have
         ;; double-line float-source data
         if n_elements(s) eq 7 && min(mgh_str_isnumber(s)) gt 0 then begin
            reads, line+text[i_line+1], src_info
            src_list->Add, src_info
            i_line = i_line+2
            continue
         endif

      endif

      ;; Process time-step info

      if process_step then begin

         ;; The following is based on the discussion in
         ;;   https://www.myroms.org/projects/src/ticket/724
         ;; I haven't considered time_ref eq -1 or -2 because I'm not sure what they mean.
         if result.prm.time_ref eq 0 then result.prm.time_ref = mgh_dt_julday('0001-01-01')

         ;; Split the line into substrings at spans of white space.
         s = strsplit(line, /EXTRACT)

         ;; If there are 7 elements, the second (index 1) matches
         ;; ????-??-??, the third (index 2) matches
         ;; "??:??:??.??", and the others are all numeric, then we have time
         ;; step data with date & time in "yyyy-mm-dd hh:mm:ss.ss" format.

         if n_elements(s) ge 7 && $
            strmatch(s[1], '????-??-??') && $
            strmatch(s[2], '??:??:??.??') && $
            min(mgh_str_isnumber([s[0],s[3:*]])) gt 0 then begin
            ;; Ad hoc fix for a Y2K leap year bug, see
            ;; https://www.myroms.org/forum/viewtopic.php?f=19&t=4521
            if s[1] eq '1999-03-**' then s[1] = '2000-02-28'
            if s[1] eq '1999-03-00' then s[1] = '2000-02-29'
            time = mgh_dt_julday(s[1]+' '+s[2]) - result.prm.time_ref
            line = strjoin(['0',s[0],mgh_format_float(time),s[3:*]], ' ')
            reads, line, step_info
            step_list->Add, step_info
            i_line++
            continue
         endif

         ;; If there are 7 elements, the third (index 2) matches
         ;; "??:??:??", and the others are all numeric, then we have time
         ;; step data with time in "day hh:mm:ss" format.

         if n_elements(s) ge 7 && $
            strmatch(s[2], '??:??:??') && $
            min(mgh_str_isnumber([s[0:1],s[3:*]])) gt 0 then begin
            dd = 0 & hh = 0 & mm = 0 & ss = 0
            reads, s[1], dd
            reads, s[2], hh, mm, ss, FORMAT='(I2,X,I2,X,I2)'
            line = strjoin(['0',s[0],string(dd+(hh+(mm+ss/60D)/60D)/24D),s[3:*]], ' ')
            reads, line, step_info
            step_list->Add, step_info
            i_line++
            continue
         endif

         ;; If there are 8 elements, the fourth (index 3) matches
         ;; "??:??:??", and the others are all numeric, then we have time
         ;; step data with time in "day hh:mm:ss" format, but with a
         ;; leading grid number (COAWST)

         if n_elements(s) ge 8 && $
            strmatch(s[3], '??:??:??') && $
            min(mgh_str_isnumber([s[0:2],s[4:*]])) gt 0 then begin
            dd = 0 & hh = 0 & mm = 0 & ss = 0
            reads, s[2], dd
            reads, s[3], hh, mm, ss, FORMAT='(I2,X,I2,X,I2)'
            line = strjoin([s[0:1],string(dd+(hh+(mm+ss/60D)/60D)/24D),s[4:*]], ' ')
            reads, line, step_info
            step_list->Add, step_info
            i_line++
            continue
         endif

         ;; If there are 5 elements, the first (index 0) matches
         ;; "(*,*,*)", and the others are all numeric, then we have
         ;; CFL data.

         if n_elements(s) ge 5 && $
            strmatch(s[0], "(*,*,*)") then begin
            reads, strjoin(s[1:*], ' '), cfl_info
            cfl_list->Add, cfl_info
            i_line++
            continue
         endif

         ;; If there are 6 elements, all numeric, then we have time
         ;; step data with time in numeric format

         if n_elements(s) ge 6 && min(mgh_str_isnumber(s)) gt 0 then begin
            reads, line, step_info
            step_list->Add, step_info
            i_line++
            continue
         endif

      endif

      i_line++

   endwhile

   ;; Post-process header info

   result.dt = double(result.date1-result.date0)*24*3600

   ;; Append float-source and time-step info, if available

   if src_list.Count() gt 0 then $
      result = create_struct(result, 'src', src_list->ToArray())

   if step_list.Count() gt 0 then $
      result = create_struct(result, 'step', step_list->ToArray())

   if cfl_list.Count() gt 0 then $
      result = create_struct(result, 'cfl', cfl_list->ToArray())

   ;; We're done here

   return, result

end

