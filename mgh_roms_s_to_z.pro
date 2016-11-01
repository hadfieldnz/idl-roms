;+
; NAME:
;   MGH_ROMS_S_TO_Z
;
; PURPOSE:
;   This function calculates the depth for a given set of s-coordinate
;   values.
;
; CATEGORY:
;   Ocean models.
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_S_TO_Z(S, H)
;
; POSITIONAL PARAMETERS:
;   s (input, scalar or array)
;     S-coordinate values, 0 <= S <= 1.
;
;   h (input, scalar)
;     Water column depth relative to Z = 0
;
; KEYWORD PARAMETERS:
;   ZETA (input, scalar)
;     Sea surface height above Z=0. Default is 0
;
;   HCRIT (input, scalar)
;     Critical depth. Default is H, which reverts to a sigma
;     coordinate.
;
;   CS (input/output, scalar/array matching s)
;     A set of stretching coefficients matching the S values. If these
;     are supplied they are used, otherwise they are calulated from
;     theta_s and theta_b and are available on output.
;
;   THETA_S (input, scalar)
;     Surface control parameter, used only if CS not supplied. Default
;     is 0.
;
;   THETA_B (input, scalar)
;     Bottom control parameter, used only if CS not supplied. Default
;     is 0.
;
; EXPLANATION:
;   The s coordinate is defined in Song & Haidvogel, 1994, J Comp Phys 115, 228-244
;   and in the ROMS Users Manual.
;
; SEE ALSO:
;   MGH_ROMS_S_TO_CS
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1998-09:
;     Written.
;   Mark Hadfield, 2001-07:
;     Updated for IDL 5.4.
;   Mark Hadfield, 2002-08:
;     Fixed error in expression for z: first term, zeta*(1+s), was
;     written zeta*(1+theta_s).
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options.
;   Mark Hadfield, 2014-07:
;     Reformatted.
;-
function mgh_roms_s_to_z, s, h, $
     ZETA=zeta, CS=cs, HC=hc, THETA_S=theta_s, THETA_B=theta_b, $
     VSTRETCH=vstretch, VTRANSFORM=vtransform

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(s) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 's'

  if n_elements(h) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'h'

  if n_elements(h) ne 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'h'

  if n_elements(hc) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'hc'

  if n_elements(vtransform) eq 0 then vtransform = 1

  ;; If CS data are not supplied, calculate them from
  ;; the vertical-stretching parameters.

  if n_elements(cs) eq 0 then begin
    cs = mgh_roms_s_to_cs(s, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)
  endif else begin
    if n_elements(theta_s) gt 0 then $
      message, 'THETA_S not required when CS is supplied'
    if n_elements(theta_b) gt 0 then $
      message, 'THETA_B not required when CS is supplied'
    if n_elements(vstretch) gt 0 then $
      message, 'VSTRETCH not required when CS is supplied'
  endelse

  ;; Calculate & return z.

  my_zeta = n_elements(zeta) gt 0 ? zeta : 0

  case vtransform of
    1: begin
      z0 = hc*(s-cs) + h*cs;
      z = z0 + my_zeta*(1+z0/h)
    end
    2: begin
      z0 = (hc*s+h*cs)/(h+hc)
      z = my_zeta + (my_zeta+h)*z0
    end
  endcase

  return, z

end
