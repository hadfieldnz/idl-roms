;+
; NAME:
;   MGH_ROMS_INTEGRAL_AREA
;
; PURPOSE:
;   Calculate area integral of a ROMS scalar field
;
; POSITIONAL ARGUMENTS:
;   a (input, numeric 2-D array)
;     Scalar field, defined at rho points
;
;   pm (input, numeric 2-D array)
;     Grid metric in xi direction, defined at rho points.
;
;   pn (input, numeric 2-D array)
;     Grid metric in eta direction, defined at rho points.
;
; KEYWORD ARGUMENTS:
;   MASK (input, numeric 2-D array)
;     Mask defined at rho points. The can be the usual land-sea mask
;     (land=0, sea=1) but it can also have fractional values to allow
;     the partial contribution of the cell values to the integral.
;
;   AREA (output, double-precision scalar)
;     Volume occupied by the scalar field in m^2.
;
; RETURN VALUE:
;   The function returns a double-precision scalar containing the integral.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-07:
;     Written as MGH_ROMS_AREA_INTEGRAL.
;   Mark Hadfield, 2014-06:
;     - Renamed MGH_ROMS_INTEGRAL_AREA.
;     - Ensure correct results with fractional mask values.
;   Mark Hadfield, 2014-08:
;     - Added NAN keyword for TOTAL function.
;-
function mgh_roms_integral_area, a, pm, pn, $
     AREA=area, MASK=mask

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if size(a, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'a'

  if size(pm, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pm'

  if size(pn, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pn'

  if size(a, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'a'

  if size(pm, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

  if size(pn, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

  dim = size(a, /DIMENSIONS)

  if ~ array_equal(size(pm, /DIMENSIONS), dim) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pm'

  if ~ array_equal(size(pn, /DIMENSIONS), dim) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pn'

  has_mask = n_elements(mask) gt 0

  if has_mask then begin
    if ~ array_equal(size(mask, /DIMENSIONS), dim) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'mask'
  endif

  ;; Calculate integral and (optionally) area

  fac = has_mask ? mask : 1

  result = total(a*fac/(pm*pn), /DOUBLE, /NaN)

  if arg_present(area) then area = total(fac/(pm*pn), /DOUBLE)

  return, result

end
