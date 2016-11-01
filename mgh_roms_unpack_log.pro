;+
; NAME:
;   MGH_ROMS_UNPACK_LOG
;
; PURPOSE:
;   From a ZIP file containing ROMS output, unpack any ROMS log files,
;   and return their contents in a string array.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   result = MGH_ROMS_UNPACK_LOG(file)
;
; POSITIONAL PARAMETERS:
;   file (input, scalar string)
;     Input file name (read-only)
;
; RETURN VALUE:
;   The function returns the contents of the log file.
;
; RESTRICTIONS:
;   - Makes assumptions about the expected file names.
;   - If there are two matching files in the archive, then the
;     temporary file will contain the contents of both.
;
; DEPENDENCIES:
;   - Requires the CMUNIQUE_ID function.
;   - The Info-Zip command-line utility "unzip" must be in the path. Under
;     Windows this can be either a Windows-native or Cygwin executable. If
;     the executable is not on the path it can be called via a wrapper batch
;     file.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2010-11:
;     Written.
;   Mark Hadfield, 2015-05:
;     - Now returns log file contents rather than the name of a temporary
;       file containing them.
;-
function mgh_roms_unpack_log, file

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(file) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file'

  if n_elements(file) gt 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'file'

  if ~ isa(file, 'STRING') then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'file'

  if ~ file_test(file, /READ) then $
    message, 'The input file cannot be read: '+file

  iswin = strcmp(!version.os_family, 'Windows', /FOLD_CASE)

  fmt = '(%"unzip -p \"%s\" rom*.log ?rom*.log qrom???.o*")'
  cmd = string(FORMAT=fmt, file)
  if iswin then begin
    spawn, /HIDE, cmd, stdout, stderr
  endif else begin
    spawn, cmd, stdout, stderr
  endelse

  return, stdout

end
