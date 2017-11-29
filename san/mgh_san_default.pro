;+
; NAME:
;   MGH_SAN_DEFAULT
;
; PURPOSE:
;   The MGH_SAN system presupposes the existence of an !MGH_SAN system variable
;   pointing to an IDL hash or orderedhash object, which defines one or more
;   directory trees, or volumes, each one being accessed via a hash key (the
;   volume label, though it doesn't actually have to be a string).
;
;   This function returns the default volume name. This can be specified with
;   another system variable, !MGH_SAN_DEFAULT. Otherwise the function returns
;   the first key in the hash object. Note that for an orderedhash object
;   the order of the keys is preserved, but for a hash object it is not.
;
; CALLING SEQUENCE:
;   result = mgh_san_default()
;
; DEPENDENCIES:
;   Requires a system variable, !MGH_SAN, to be defined, as described
;   in the documentation.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2017-11:
;     Written.
;-
function mgh_san_default

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   defsysv, '!MGH_SAN', EXISTS=exists

   if ~ exists then $
      message, 'The !MGH_SAN system variable needs to be set.'

   defsysv, '!MGH_SAN_DEFAULT', EXISTS=exists

   if exists then begin

      return, !MGH_SAN_DEFAULT

   endif else begin

      return, (!MGH_SAN.keys())[0]

   endelse

end
