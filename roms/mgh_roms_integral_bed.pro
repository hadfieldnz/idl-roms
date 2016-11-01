;+
; NAME:
;   MGH_ROMS_INTEGRAL_BED
;
; PURPOSE:
;   Calculate the bed integral of a ROMS sediment variable. This must
;   be a variable like sandmass_??/mudmass_??/bed_thickness which
;   is integrated by summing over the bed layers.
;
; POSITIONAL ARGUMENTS:
;   a (input, numeric 2-D or 3-D array)
;     Scalar field, defined at rho points. It will normally have 3 dimensions
;     (xi,eta,bed) but the code also has to allow for a unit trailing dimension
;     to be omitted when the number of bed layers is 1.
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
; RETURN VALUE:
;   The function returns a double-precision scalar containing the integral.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-07:
;     Written as MGH_ROMS_BED_INTEGRAL.
;   Mark Hadfield, 2014-06:
;     - Renamed MGH_ROMS_INTEGRAL_BED.
;     - Ensure correct results with fractional mask values.
;   Mark Hadfield, 2016-07:
;     - Reformatted.
;-
function mgh_roms_integral_bed, a, pm, pn, MASK=mask

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

   ;; The number of dimensions of a
   if size(a, /N_DIMENSIONS) lt 2 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'a'

   if size(pm, /N_DIMENSIONS) ne 2 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

   if size(pn, /N_DIMENSIONS) ne 2 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

   dim = size(a, /DIMENSIONS)

   if ~ array_equal(size(pm, /DIMENSIONS), dim[0:1]) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pm'

   if ~ array_equal(size(pn, /DIMENSIONS), dim[0:1]) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pn'

   has_mask = n_elements(mask) gt 0

   if has_mask then begin
      if ~ array_equal(size(mask, /DIMENSIONS), dim[0:1]) then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'mask'
   endif

   ;; Accumulate integral.

   result = 0D

   for j=1,dim[1]-2 do begin
      for i=1,dim[0]-2 do begin
         fac = has_mask ? mask[i,j] : 1
         if fac gt 0 then begin
            result = result + total(fac*a[i,j,*]/(pm[i,j]*pn[i,j]), /DOUBLE)
         endif
      endfor
   endfor

   return, result

end
