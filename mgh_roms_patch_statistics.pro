;+
; FUNCTION NAME:
;   MGH_ROMS_PATCH_STATISTICS
;
; PURPOSE:
;   This function calculates and returns the area and volume of a polygonal
;   subset of the domain.
;
;   The polygonal subset is defined with the MGHromsHistory::PatchGrid method.
;   This defaults to encompassing the entire domain.
;
; CALLING SEQUENCE
;   result = mgh_roms_patch_statistics(ohis)
;
; RETURN VALUE:
;   The function returns the result in a structure.
;
; POSITIONAL PARAMETERS:
;   ohis (input, object)
;     A reference to an MGHromsHistory object.
;
; KEYWORD PARAMETERS:
;   LONLAT (input, switch)
;    This keyword specicifies whether the vertext locations (VERTX and VERTY)
;    are defined in (x,y) or (lon,lat). The default is 1 if the history
;    file has (lon,lat) data and 0 otherwise.
;
;   VERTX (input, numeric vector)
;   VERTY (input, numeric vector)
;     Vertex locations defining the polygonal subset.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2014-07:
;     Written.
;-
function mgh_roms_patch_statistics, ohis, $
     LONLAT=lonlat, VERTX=vertx, VERTY=verty

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Process file argument

  if n_elements(ohis) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'

  if n_elements(ohis) gt 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'ohis'

  if ~ obj_valid(ohis) gt 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ohis'

  ;; Process the info defining the polygonal Patch. The default will enclose
  ;;  all interior points

  patch = ohis->PatchGrid(LONLAT=lonlat, VERTX=vertx, VERTY=verty)

  xr0 = patch.xi_range[0]
  xr1 = patch.xi_range[1]
  xrn = xr1 - xr0 + 1

  er0 = patch.eta_range[0]
  er1 = patch.eta_range[1]
  ern = er1 - er0 + 1

  ;; Get horizontal grid data

  mask = patch.frac

  h = ohis->VarGet('h', OFFSET=[xr0,er0], COUNT=[xrn,ern])

  pm = ohis->VarGet('pm', OFFSET=[xr0,er0], COUNT=[xrn,ern])
  pn = ohis->VarGet('pn', OFFSET=[xr0,er0], COUNT=[xrn,ern])

  return, {area: total(mask/(pm*pn), /DOUBLE, /NAN), $
           vol0: total(mask*h/(pm*pn), /DOUBLE, /NAN)}

end
