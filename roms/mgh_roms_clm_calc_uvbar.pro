;+
; NAME:
;   MGH_ROMS_CLM_CALC_UVBAR
;
; PURPOSE:
;   Given a ROMS climatology file with 3D velocity data, calculate &
;   write vertical averages
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_CALC_UVBAR, file_clm
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file to be processed.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-05:
;     Written.
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options, but introduced a bug: VSTRETCH and VTRANSFORM not read
;     from file!
;   Mark Hadfield, 2009-10:
;     - Removed calls to widget_event(/NOWAIT).
;     - Fixed bug introduced in 2009-04
;   Mark Hadfield, 2009-11:
;     - Fixed bug: 2D zeta array passed to MGH_ROMS_S_TO_Z instead of scalar.
;   Mark Hadfield, 2010-02:
;     Added PACK_DATA functionality.
;   Mark Hadfield, 2013-07:
;     Default TIME_NAME now 'ocean_time'
;   Mark Hadfield, 2013-10:
;     Cleaned up and reformatted while writing MGH_ROMS_BRY_CALC_UVBAR.
;   Mark Hadfield, 2015-11:
;     Corrected a latent bug in setting dimensions of ubar and vbar: the bug
;     (using dim instead of dim[0:1]) is latent because the array dimensions
;     are forced to 2D by the subraction of a 2D expression.
;   Mark Hadfield, 2017-07:
;     Added a missing_value attribute for packed data, for the benefit of applications
;     that do not respect valid_range.
;-
pro mgh_roms_clm_calc_uvbar, file_clm, $
     PACK_DATA=pack_data, TIME_NAME=time_name

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Check file name etc

  if n_elements(file_clm) eq 0 then $
    message, 'Name for climatology file not supplied'

  if ~ file_test(file_clm, /READ) then $
    message, 'Climatology file cannot be read'

  if ~ file_test(file_clm, /WRITE) then $
    message, 'Climatology file cannot be written to'

  ;; Process keyword parameters

  if n_elements(time_name) eq 0 then time_name = 'ocean_time'

  ;; Open the file

  fmt = '(%"Modifying climatology file %s")'
  message, /INFORM, string(FORMAT=fmt, file_clm)

  oclm = obj_new('MGHncFile', file_clm, /MODIFY)

  ;; Check that the required time dimension & variable are available

  if ~ oclm->HasDim(time_name) then message, 'Time dimension missing'
  if ~ oclm->HasVar(time_name) then message, 'Time variable missing'

  n_time = oclm->DimInfo(time_name, /DIMSIZE)

  ;; Get grid data

  dim = [oclm->DimInfo('xi_rho', /DIMSIZE), $
         oclm->DimInfo('eta_rho', /DIMSIZE), $
         oclm->DimInfo('s_rho', /DIMSIZE)]

  h = oclm->VarGet('h')

  s_w = oclm->VarGet('s_w')

  theta_s =oclm->VarGet('theta_s')
  theta_b = oclm->VarGet('theta_b')

  hc = oclm->VarGet('hc')

  vstretch = oclm->HasVar('Vstretching') ? oclm->VarGet('Vstretching') : 1
  vtransform = oclm->HasVar('Vtransform') ? oclm->VarGet('Vtransform') : 1

  h_u = mgh_roms_stagger(h, TO_GRID='u')
  h_v = mgh_roms_stagger(h, TO_GRID='v')

  cs = mgh_roms_s_to_cs(s_w, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

  ;; Add vertical-average velocity variables

  uv_range = [-10,10]
  if keyword_set(pack_data) then begin
    pack_range = [-32766S,32767S]
    scale = mgh_norm_coord(pack_range, uv_range)
  endif

  message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'ubar')
  if keyword_set(pack_data) then begin
    oclm->VarAdd, 'ubar', ['xi_u','eta_u',time_name], /SHORT
    oclm->attadd, 'ubar', 'add_offset', scale[0]
    oclm->attadd, 'ubar', 'scale_factor', scale[1]
    oclm->attadd, 'ubar', 'valid_range', pack_range
    oclm->attadd, 'ubar', 'missing_value', pack_range[0]-1S
  endif else begin
    oclm->VarAdd, 'ubar', ['xi_u','eta_u',time_name], /FLOAT
    oclm->attadd, 'ubar', 'valid_range', uv_range
  endelse
  oclm->AttAdd, 'ubar', 'time', time_name

  message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'vbar')
  if keyword_set(pack_data) then begin
    oclm->VarAdd, 'vbar', ['xi_v','eta_v',time_name], /SHORT
    oclm->attadd, 'vbar', 'add_offset', scale[0]
    oclm->attadd, 'vbar', 'scale_factor', scale[1]
    oclm->attadd, 'vbar', 'valid_range', pack_range
    oclm->attadd, 'vbar', 'missing_value', pack_range[0]-1S
  endif else begin
    oclm->VarAdd, 'vbar', ['xi_v','eta_v',time_name], /FLOAT
    oclm->attadd, 'vbar', 'valid_range', uv_range
  endelse
  oclm->AttAdd, 'vbar', 'time', time_name

  ;; For each time, retrieve velocity and zeta if available)

  for r=0,n_time-1 do begin

    message, /INFORM, string(FORMaT='(%"Getting u, v & zeta for record %d")', r)

    uu = oclm->VarGet('u', OFFSET=[0,0,0,r], COUNT=[0,0,0,1], /AUTOSCALE)
    vv = oclm->VarGet('v', OFFSET=[0,0,0,r], COUNT=[0,0,0,1], /AUTOSCALE)

    if oclm->HasVar('zeta') then begin
      zeta = oclm->VarGet('zeta', OFFSET=[0,0,r], COUNT=[0,0,1], /AUTOSCALE)
      zeta_u = mgh_roms_stagger(zeta, TO_GRID='u')
      zeta_v = mgh_roms_stagger(zeta, TO_GRID='v')
      mgh_undefine, zeta
    endif else begin
      zeta_u = fltarr(dim[0:1]-[1,0])
      zeta_v = fltarr(dim[0:1]-[0,1])
    endelse

    message, /INFORM, string(FORMaT='(%"Calculating ubar & vbar for record %d")', r)

    ubar = fltarr(dim[0:1]-[1,0])
    vbar = fltarr(dim[0:1]-[0,1])

    for i=0,dim[0]-2 do begin
      for j=0,dim[1]-1 do begin
        zz = mgh_roms_s_to_z(s_w, h_u[i,j], ZETA=zeta_u[i,j], CS=cs, HC=hc, VTRANSFORM=vtransform)
        dz = mgh_diff(zz)
        ubar[i,j] = total(dz*reform(uu[i,j,*]))/total(dz)
        mgh_undefine, zz, dz
      endfor
    endfor

    for i=0,dim[0]-1 do begin
      for j=0,dim[1]-2 do begin
        zz = mgh_roms_s_to_z(s_w, h_v[i,j], ZETA=zeta_v[i,j], CS=cs, HC=hc, VTRANSFORM=vtransform)
        dz = mgh_diff(zz)
        vbar[i,j] = total(dz*reform(vv[i,j,*]))/total(dz)
        mgh_undefine, zz, dz
      endfor
    endfor

    if keyword_set(pack_data) then begin
      ubar = round((temporary(ubar)-scale[0])/scale[1])
      vbar = round((temporary(vbar)-scale[0])/scale[1])
    endif

    message, /INFORM, string(FORMaT='(%"Writing ubar & vbar for record %d")', r)

    oclm->VarPut, 'ubar', temporary(ubar), OFFSET=[0,0,r]
    oclm->VarPut, 'vbar', temporary(vbar), OFFSET=[0,0,r]

  endfor

  ;; Close the output file

  obj_destroy, oclm

end

