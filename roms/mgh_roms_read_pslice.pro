;+
; NAME:
;   MGH_ROMS_READ_PSLICE
;
; PURPOSE:
;   For a specified grid, read grid-vertex data from a file and return
;   the corresponding MGHromsHistory::PsliceGrid structure
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2019-05:
;     Written as a generalised form of acee_inter_curvi_pslice.
;-
function mgh_roms_read_pslice, ogrd, file, NAME=name

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(ogrd) eq 0 then $
      message, 'Undefined variable: ogrd'

   if n_elements(file) eq 0 then $
      message, 'Undefined variable: file'

   if ~ file_test(file, /READ) then $
      message, 'File cannot be read: '+file

   vert = list()

   line = ''  &  data = lonarr(2)

   openr, lun, file, /GET_LUN

   while ~eof(lun) do begin
      readf, lun, line
      if mgh_str_iswhite(line) then continue
      if strmatch(line, '#*') then continue
      reads, line, data
      vert->Add, data
   endwhile

   free_lun, lun

   vert = vert->ToArray()

   return, ogrd->PsliceGrid(NAME=name, VERT_XI=vert[*,0], VERT_ETA=vert[*,1])

end
