;+
; NAME:
;   MGH_ROMS_READ_FLTIN
;
; PURPOSE:
;   Read info from a ROMS floats-input file.
;
;   Note that function MGH_ROMS_READ_LOG also reads float-source info, but
;   from the output log file. There are 2 reasons why MGH_ROMS_READ_FLTIN
;   might be preferred: it can read trailing comments that are ignored by
;   ROMS; the info in the log file might have been written at lower precision.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_READ_FLTIN(file)
;
; POSITIONAL PARAMETERS:
;   file (input, scalar string)
;     Input file name (read-only)
;
; KEYWORD PARAMETERS:
;   READ_COMMENT (input, switch)
;     This keyword specifies whether the function attempts to read
;     string data at the end of each source-definition line and save
;     it in the result's comment field.  The default is to read the
;     comment, but this will fail if there is no string data after the
;     last numeric value (FDZ). Set the keyword to zero to override
;     the default.
;
; RETURN VALUE:
;   The function returns a structure with a whole load of different tags.
;
; RESTRICTIONS:
;   Currently handles only one grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-05:
;     Written.
;   Mark Hadfield, 2009-11:
;     - Added READ_COMMENT keyword.
;     - Position and time data now double-precision.
;-
function mgh_roms_read_fltin, file, READ_COMMENT=read_comment

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(file) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file'

   if n_elements(file) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'file'

   if size(file, /TYPE) ne 7 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'file'

   if ~ file_test(file, /READ) then $
        message, 'The input file cannot be read: '+file

   if n_elements(read_comment) eq 0 then read_comment = 1B

   openr, lun, file, /GET_LUN, COMPRESS=strmatch(file, '*.gz')

   line = ''
   n_line = 0
   l_float = 0B
   n_float = 0

   if keyword_set(read_comment) then begin
      src_info = {g: 0, c: 0, t: 0, n: 0, $
                  ft0: 0D, fx0: 0D, fy0: 0D, fz0: 0D, $
                  fdt: 0D, fdx: 0D, fdy: 0D, fdz: 0D, comment: ''}
   endif else begin
      src_info = {g: 0, c: 0, t: 0, n: 0, $
                  ft0: 0D, fx0: 0D, fy0: 0D, fz0: 0D, $
                  fdt: 0D, fdx: 0D, fdy: 0D, fdz: 0D}
   endelse

   stage = 0

   osrc = obj_new('MGH_Vector')

   while ~ eof(lun) do begin

      readf, lun, line

      n_line += 1

      if (stage eq 0) && strmatch(line, '*POS*=*', /FOLD_CASE) then begin
         stage = 1
         continue
      end

      case stage of

         0: begin
            if strmatch(line, '*NFLOATS*=*', /FOLD_CASE) then begin
               s = strsplit(line, /EXTRACT)
               reads, s[2], n_float
            endif
            if strmatch(line, '*LFLOATS*=*', /FOLD_CASE) then begin
               s = strsplit(line, /EXTRACT)
               n_grid = n_elements(s)-2
               if n_grid ne 1 then $
                    message, 'Invalid number of grids'
               l_float = strmatch(s[2], 'T*', /FOLD_CASE)
            endif
         end

         1: begin

            if mgh_str_iswhite(line) then break
            if strmatch(line, '!*', /FOLD_CASE) then break
            reads, line, src_info
            if keyword_set(read_comment) then $
                 src_info.comment = strtrim(src_info.comment, 2)
            osrc->Add, src_info

         end

      endcase

   endwhile

   free_lun, lun

   n_src = osrc->Count()

   result = create_struct('n_line', n_line, 'l_float', l_float, $
                          'n_float', n_float, 'n_src', n_src)

   if n_src gt 0 then $
        result = create_struct(result, 'src', osrc->ToArray())

   obj_destroy, osrc

   return, result

end

