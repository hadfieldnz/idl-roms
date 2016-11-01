; svn $Id$
;+
; NAME:
;   MGH_ROMS_FILE_TYPE
;
; PURPOSE:
;   Open a ROMS netCDF file (or sequence thereof), read and return the value
;   of the global attribute "type"
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_FILE_TYPE(file)
;
; POSITIONAL PARAMETERS:
;   file (input, scalar string)
;     Input file name (read-only)
;
; RETURN VALUE:
;   The function returns a scalar string containing the type attribute. Values
;   currently supported by ROMS include:
;
;     "ROMS hystory file" (sic)
;     "ROMS restart file"
;     "ROMS station file"
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-04:
;     Written.
;   Mark Hadfield, 2010-09:
;     The value returned when the type attribute does not exist is
;     now an empty string.
;-

function mgh_roms_file_type, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if size(file, /TNAME) ne 'STRING' then $
        message, 'Invalid or missing file argument'

   onc = obj_new('mghncsequence', file)

   result = ''
   if onc->HasAtt('type', /GLOBAL) then $
        result = onc->AttGet('type', /GLOBAL)

   obj_destroy, onc

   return, result

end

