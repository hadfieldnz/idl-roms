;+
; NAME:
;   MGH_ROMS_BRY_NEW
;
; PURPOSE:
;   Create a ROMS boundary file, given a ROMS grid file or similar
;   for grid information.
;
; CALLING SEQUENCE:
;   MGH_ROMS_BRY_NEW, file_bry, file_grd
;
; POSITIONAL PARAMETERS:
;   file_bry (input, scalar string)
;     The name of a ROMS boundary file to be written.
;
;   file_grd (input, scalar string)
;     The name of a ROMS climatology file to be read.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-03:
;     Written.
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options, but introduced a bug: VSTRETCH and VTRANSFORM not read
;     from file!
;   Mark Hadfield, 2009-10:
;     Fixed bug introduced in 2009-04
;   Mark Hadfield, 2013-10:
;     - Used this for the first time (apparently) in 4 years. Fixed typo; re-indented
;       and cleaned up.
;     - Labelled this routine obsolete, as new boundary files can be created with
;       MGH_ROMS_CLM_NEW.
;-
pro mgh_roms_bry_new, file_bry, file_grd, $
     SCOORD_PARAMETERS=scoord_parameters, VERTICAL=vertical

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  compile_opt OBSOLETE

  ;; Check file name etc

  if n_elements(file_bry) eq 0 then $
    message, 'Name for boundary file not supplied'

  if n_elements(file_grd) eq 0 then $
    message, 'Name for template file not supplied'

  if ~ file_test(file_grd, /READ) then $
    message, 'Template file cannot be read'

  ;; Open the files

  message, /INFORM, string(FORMAT='(%"Opening grid file %s")', file_grd)

  ogrd = obj_new('MGHncReadFile', file_grd)

  message, /INFORM, string(FORMAT='(%"Creating boundary file %s")', file_bry)

  obry = obj_new('MGHncFile', file_bry, /CREATE, /CLOBBER)

  obry->AttAdd, /GLOBAL, 'type', 'ROMS boundary file'

  ;; Specify defaults--could get s-coordinate defaults from template,
  ;; if it has them.

  if n_elements(vertical) eq 0 then vertical = 1B

  if keyword_set(vertical) then begin

    if n_elements(scoord_parameters) eq 0 then $
      scoord_parameters = {n: 20, hc: 20.0D, theta_s: 4.0D, theta_b: 0.2D, $
      vstretch: 2, vtransform: 2}

    sc = scoord_parameters

    if sc.vtransform eq 1 && (sc.hc gt min(ogrd->VarGet('h'))) then $
      message, 'S-coordinate critical depth exceeds minimum depth'

  endif

  ;; Specify dimensions

  message, /INFORM, 'Specifying boundary file dimensions'

  obry->DimCopy, ogrd, $
    ['xi_rho','eta_rho','xi_u','eta_u','xi_v','eta_v']

  if keyword_set(vertical) then begin
    obry->DimAdd, 's_rho', sc.n
    obry->DimAdd, 's_w', sc.n + 1
  endif

  ;; Specify grid variables

  if keyword_set(vertical) then begin

    message, /INFORM, 'Specifying boundary file vertical grid variables'

    obry->VarAdd, 'theta_s'
    obry->VarAdd, 'theta_b'
    obry->VarAdd, 'hc'
    obry->VarAdd, 'Vstretching'
    obry->VarAdd, 'Vtransform'
    obry->VarAdd, 'sc_r', ['s_rho']
    obry->VarAdd, 'Cs_r', ['s_rho']
    obry->VarAdd, 'sc_w', ['s_w']
    obry->VarAdd, 'Cs_w', ['s_w']

  end

  ;; Write vertical grid data

  if keyword_set(vertical) then begin

    message, /INFORM, 'Writing vertical grid data'

    obry->VarPut, 'theta_b', sc.theta_b
    obry->VarPut, 'theta_s', sc.theta_s
    obry->VarPut, 'hc', sc.hc
    obry->VarPut, 'Vstretching', sc.vstretch
    obry->VarPut, 'Vtransform', sc.vtransform

    sc_w = mgh_range(-1, 0, N_ELEMENTS=sc.n+1)
    sc_r = mgh_stagger(sc_w, DELTA=-1)

    cs_r = mgh_roms_s_to_cs(sc_r, THETA_B=sc.theta_b, THETA_S=sc.theta_s, VSTRETCH=sc.vstretch)

    obry->VarPut, 'sc_r', sc_r
    obry->VarPut, 'Cs_r', cs_r

  endif

  ;; Clean up

  obj_destroy, [obry,ogrd]

end

