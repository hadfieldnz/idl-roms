;+
; NAME:
;   MGH_ROMS_CLM_NEW
;
; PURPOSE:
;   Create an empty climatology file given a ROMS grid file (or
;   similar) for grid info.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_NEW, file_clm, file_grd
;
; POSITIONAL ARGUMENTSS:
;   file_clm (input, scalar string)
;     The name of a ROMS climate file to be created (write-only)
;
;   history (input, string scalar/vector OR object reference)
;     The ROMS history-file sequence from which grid data are to be read, either
;     as a list of file names or an MGHromsHistory object.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-01:
;     Written.
;   Mark Hadfield, 2004-09:
;     Added check to ensure S-coordinate critical depth does not exceed
;     minimum depth in bathymetry.
;   Mark Hadfield, 2007-01:
;     S-coordinate variables renamed from sc_r & sc_w to s_rho & s_w.
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options.
;   Mark Hadfield, 2011-07:
;     Added N_BED keyword
;   Mark Hadfield, 2013-08:
;     - Added TYPE keyword
;     - Extra grid variables copied to output file: lon_u, lat_u, lon_v, lat_v.
;   Mark Hadfield, 2013-08:
;     - The history file may now be passed to this routine as an
;       MGHromsHistory object.
;   Mark Hadfield, 2014-03:
;     - Files now created with NETCDF_FORMAT keyword.
;-
pro mgh_roms_clm_new, file_clm, history, $
     N_BED=n_bed, SCOORD_PARAMETERS=scoord_parameters, $
     TYPE=type, VERTICAL=vertical

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Check file name etc

  if n_elements(file_clm) eq 0 then $
    message, 'Name for climate file not supplied'

  if n_elements(history) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'

  ;; Open the files

  message, /INFORM, string(FORMAT='(%"Creating climate file %s")', file_clm)

  oclm = obj_new('MGHncFile', file_clm, /CREATE, /CLOBBER, /NETCDF4_FORMAT)

  case size(history, /TNAME) of
    'STRING': begin
      fmt = '(%"Opening history file sequence beginning with %s")'
      message, /INFORM, string(FORMAT=fmt, history[0])
      ohis = obj_new('MGHromsHistory', history)
      history_destroy = 1B
    end
    'OBJREF': begin
      ohis = history
    end
    else: $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
  endcase

  ;; Specify defaults

  if n_elements(n_bed) eq 0 then n_bed = 0

  if n_elements(vertical) eq 0 then vertical = 1B

  if n_elements(type) eq 0 then type = 'ROMS climatology file'

  if keyword_set(vertical) then begin

    if n_elements(scoord_parameters) eq 0 then $
      scoord_parameters = {n: 20, hc: 20.0D, theta_s: 4.0D, theta_b: 0.2D, $
      vstretch: 2, vtransform: 2}

    sp = scoord_parameters

    if sp.vtransform eq 1 && (sp.hc gt min(ohis->VarGet('h'))) then $
      message, 'S-coordinate critical depth exceeds minimum depth for transform 1'

  endif

  ;; Specify dimensions

  message, /INFORM, 'Specifying climate file dimensions'

  oclm->DimCopy, ohis, $
    ['xi_rho','eta_rho','xi_u','eta_u','xi_v','eta_v']

  if keyword_set(vertical) then begin
    oclm->DimAdd, 's_rho', sp.n
    oclm->DimAdd, 's_w', sp.n + 1
  endif

  if n_bed gt 0 then begin
    oclm->DimAdd, 'Nbed', n_bed
  endif

  ;; Specify variables

  vars_grd = ['xl','el','h','f','pm','pn','lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v','angle']

  if ohis->HasVar('mask_rho') && ohis->HasVar('mask_v') && ohis->HasVar('mask_v') then $
    vars_grd = [vars_grd,'mask_rho','mask_u','mask_v']

  oclm->VarCopy, ohis, /DEFINITION, vars_grd

  if keyword_set(vertical) then begin
    oclm->VarAdd, 'theta_s', /DOUBLE
    oclm->VarAdd, 'theta_b', /DOUBLE
    oclm->VarAdd, 'hc', /DOUBLE
    oclm->VarAdd, 's_rho', ['s_rho'], /DOUBLE
    oclm->VarAdd, 'Cs_r', ['s_rho'], /DOUBLE
    oclm->VarAdd, 's_w', ['s_w'], /DOUBLE
    oclm->VarAdd, 'Cs_w', ['s_w'], /DOUBLE
    oclm->VarAdd, 'Vstretching', /SHORT
    oclm->VarAdd, 'Vtransform', /SHORT
  end

  oclm->VarCopy, ohis, /ATTRIBUTES, vars_grd

  oclm->AttAdd, /GLOBAL, 'type', type

  fmt = '(%"Created by procedure mgh_roms_clm_new at %s")'
  oclm->AttAdd, /GLOBAL, 'history', $
    string(FORMAT=fmt, mgh_dt_string(mgh_dt_now()))

  ;; Copy horizontal grid data

  message, /INFORM, 'Copying grid data'

  oclm->VarCopy, ohis, /DATA, vars_grd

  ;; Write vertical grid data

  if keyword_set(vertical) then begin

    message, /INFORM, 'Writing vertical grid data'

    oclm->VarPut, 'hc', sp.hc
    oclm->VarPut, 'theta_b', sp.theta_b
    oclm->VarPut, 'theta_s', sp.theta_s
    oclm->VarPut, 'Vtransform', sp.vtransform
    oclm->VarPut, 'Vstretching', sp.vstretch

    s_w = mgh_range(-1D, 0D, N_ELEMENTS=sp.n+1)
    s_r = mgh_stagger(s_w, DELTA=-1)

    cs_r = mgh_roms_s_to_cs(s_r, $
      THETA_B=sp.theta_b, THETA_S=sp.theta_s, $
      VSTRETCH=sp.vstretch)
    cs_w = mgh_roms_s_to_cs(s_w, $
      THETA_B=sp.theta_b, THETA_S=sp.theta_s, $
      VSTRETCH=sp.vstretch)

    oclm->VarPut, 's_rho', s_r
    oclm->VarPut, 's_w', s_w
    oclm->VarPut, 'Cs_r', cs_r
    oclm->VarPut, 'Cs_w', cs_w

  endif

  ;; Clean up

  if keyword_set(history_destroy) then obj_destroy, ohis
  obj_destroy, oclm

end

