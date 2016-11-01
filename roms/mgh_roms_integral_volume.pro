;+
; NAME:
;   MGH_ROMS_INTEGRAL_VOLUME
;
; PURPOSE:
;   Calculate volume integral of a ROMS scalar field
;
; POSITIONAL ARGUMENTS:
;   a (input, numeric 3-D array)
;     Scalar field, defined at rho points
;
;   pm (input, numeric 2-D array)
;     Grid metric in xi direction, defined at rho points.
;
;   pn (input, numeric 2-D array)
;     Grid metric in eta direction, defined at rho points.
;
;   h (input, numeric 2-D array)
;     Bathymetry, defined at rho points.
;
;   scoord (input, structure)
;     S-coordinate parameters
;
; KEYWORD ARGUMENTS:
;   AREA (output, double-precision scalar)
;     Area occupied by the vertical projection of the scalar field in m^2.
;
;   MASK (input, numeric 2-D array)
;     Mask defined at rho points. The can be the usual land-sea mask
;     (land=0, sea=1) but it can also have fractional values to allow
;     the partial contribution of the cell values to the integral.
;
;   VOLUME (output, double-precision scalar)
;     Volume occupied by the scalar field in m^3.
;
;   ZETA (input, numeric 2-D array)
;     Sea surface height.
;
; RETURN VALUE:
;   The function returns a double-precision scalar containing the integral.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2011-08:
;     Written as MGH_ROMS_VOLUME_INTEGRAL.
;   Mark Hadfield, 2014-06:
;     - Renamed MGH_ROMS_INTEGRAL_VOLUME.
;     - Ensure correct results with fractional mask values.
;-
function mgh_roms_integral_volume, a, pm, pn, h, scoord, $
     AREA=area, MASK=mask, VOLUME=volume, ZETA=zeta

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

  if size(h, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'h'

  if size(a, /N_DIMENSIONS) ne 3 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'a'

  if size(pm, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

  if size(pn, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

  if size(h, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'h'

  dim = size(a, /DIMENSIONS)

  if ~ array_equal(size(pm, /DIMENSIONS), dim[0:1]) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pm'

  if ~ array_equal(size(pn, /DIMENSIONS), dim[0:1]) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'pn'

  if ~ array_equal(size(h, /DIMENSIONS), dim[0:1]) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'h'

  has_mask = n_elements(mask) gt 0

  if has_mask then begin
    if ~ array_equal(size(mask, /DIMENSIONS), dim[0:1]) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'mask'
  endif

  has_zeta = n_elements(zeta) gt 0

  if has_zeta then begin
    if ~ array_equal(size(zeta, /DIMENSIONS), dim[0:1]) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'a', 'mask'
  endif

  ;; Optionally calculate area
  if arg_present(area) then begin
    b = 1
    if has_mask then b *= mask
    area = total(b/(pm*pn), /DOUBLE)
  endif

  ;; Optionally calculate volume

  if arg_present(volume) then begin
    b = h
    if has_zeta then b += zeta
    if has_mask then b *= mask
    volume = total(b/(pm*pn), /DOUBLE, /NAN)
  endif

  ;; Process s-coordinate info then calculate s, cs, and later z, at
  ;; w levels.

  theta_s = scoord.theta_s
  theta_b = scoord.theta_b
  hc = scoord.hc
  vstretch = mgh_struct_has_tag(scoord, 'vstretch') ? scoord.vstretch : 1
  vtransform = mgh_struct_has_tag(scoord, 'vtransform') ? scoord.vtransform : 1

  s = mgh_range(-1, 0, N_ELEMENTS=dim[2]+1)

  cs = mgh_roms_s_to_cs(s, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

  ;; Accumulate integral.

  result = 0D

  for j=1,dim[1]-2 do begin
    for i=1,dim[0]-2 do begin
      fac = has_mask ? mask[i,j] : 1
      if fac gt 0 then begin
        if has_zeta then begin
          zz = mgh_roms_s_to_z(s, h[i,j], ZETA=zeta[i,j], $
                               CS=cs, HC=hc, VTRANSFORM=vtransform)
        endif else begin
          zz = mgh_roms_s_to_z(s, h[i,j], $
                               CS=cs, HC=hc, VTRANSFORM=vtransform)
        endelse
        result = result + total(fac*a[i,j,*]*mgh_diff(zz)/(pm[i,j]*pn[i,j]), /DOUBLE)
      endif
    endfor
  endfor

  return, result

end
