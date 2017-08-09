;+
; NAME:
;   MGH_ROMS_bry_CALC_UVBAR
;
; PURPOSE:
;   Given a ROMS boundary file with 3D velocity data, calculate &
;   write vertical averages
;
; CALLING SEQUENCE:
;   MGH_ROMS_bry_CALC_UVBAR, file_bry
;
; POSITIONAL PARAMETERS:
;   file_bry (input, scalar string)
;     The name of a ROMS boundary file to be processed.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-10:
;     Written.
;   Mark Hadfield, 2017-07:
;     Added a missing_value attribute for packed data, for the benefit of applications
;     that do not respect valid_range.
;-
pro mgh_roms_bry_calc_uvbar, file_bry, $
     PACK_DATA=pack_data, TIME_NAME=time_name

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Check file name etc

  if n_elements(file_bry) eq 0 then $
    message, 'Name for boundary file not supplied'

  if ~ file_test(file_bry, /READ) then $
    message, 'Boundary file cannot be read'

  if ~ file_test(file_bry, /WRITE) then $
    message, 'Boundary file cannot be written to'

  ;; Process keyword parameters

  if n_elements(time_name) eq 0 then time_name = 'ocean_time'

  ;; Open the file

  fmt = '(%"Modifying boundary file %s")'
  message, /INFORM, string(FORMAT=fmt, file_bry)

  obry = obj_new('MGHncFile', file_bry, /MODIFY)

  ;; Check that the required time dimension & variable are available

  if ~ obry->HasDim(time_name) then message, 'Time dimension missing'
  if ~ obry->HasVar(time_name) then message, 'Time variable missing'

  n_time = obry->DimInfo(time_name, /DIMSIZE)

  ;; Get grid data

  dim = [obry->DimInfo('xi_rho', /DIMSIZE), $
         obry->DimInfo('eta_rho', /DIMSIZE), $
         obry->DimInfo('s_rho', /DIMSIZE)]

  h = obry->VarGet('h')

  h_u = mgh_roms_stagger(h, TO_GRID='u')
  h_v = mgh_roms_stagger(h, TO_GRID='v')

  s_w = obry->VarGet('s_w')

  theta_s =obry->VarGet('theta_s')
  theta_b = obry->VarGet('theta_b')

  hc = obry->VarGet('hc')

  vstretch = obry->HasVar('Vstretching') ? obry->VarGet('Vstretching') : 1
  vtransform = obry->HasVar('Vtransform') ? obry->VarGet('Vtransform') : 1

  cs = mgh_roms_s_to_cs(s_w, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

  ;; Check for the existence of velocity and zeta data on each boundary

  bname = ['south','east','north','west']

  has_uv = bytarr(4)
  has_zeta = bytarr(4)

  for b=0,3 do begin

    has_uv[b] = obry->HasVar('u_'+bname[b]) && obry->HasVar('v_'+bname[b])
    has_zeta[b] = obry->HasVar('zeta_'+bname[b])

  endfor

  ;; For each boundary, add the appropriate vertical velocity variables,
  ;; provided the respective velocity variables are available

  uv_range = [-10,10]
  if keyword_set(pack_data) then begin
    pack_range = [-32766S,32767S]
    uv_scale = mgh_norm_coord(pack_range, uv_range)
  endif

  for b=0,3 do begin

    if has_uv[b] then begin

      var_ubar = 'ubar_'+bname[b]
      var_vbar = 'vbar_'+bname[b]

      case b mod 2 of
        0: begin
          dim_u = 'xi_u'
          dim_v = 'xi_v'
        end
        1: begin
          dim_u = 'eta_u'
          dim_v = 'eta_v'
        end
      endcase

      message, /INFORM, 'Adding variable '+var_ubar
      if keyword_set(pack_data) then begin
        obry->VarAdd, var_ubar, [dim_u,time_name], /SHORT
        obry->AttAdd, var_ubar, 'add_offset', uv_scale[0]
        obry->AttAdd, var_ubar, 'scale_factor', uv_scale[1]
        obry->AttAdd, var_ubar, 'valid_range', pack_range
        obry->AttAdd, var_ubar, 'missing_value', pack_range[0]-1S
      endif else begin
        obry->VarAdd, var_ubar, [dim_u,time_name], /FLOAT
        obry->AttAdd, var_ubar, 'valid_range', uv_range
      endelse
      obry->AttAdd, var_ubar, 'time', time_name

      message, /INFORM, 'Adding variable '+var_vbar
      if keyword_set(pack_data) then begin
        obry->VarAdd, var_vbar, [dim_v,time_name], /SHORT
        obry->AttAdd, var_vbar, 'add_offset', uv_scale[0]
        obry->AttAdd, var_vbar, 'scale_factor', uv_scale[1]
        obry->AttAdd, var_vbar, 'valid_range', pack_range
        obry->AttAdd, var_vbar, 'missing_value', pack_range[0]-1S
      endif else begin
        obry->VarAdd, var_vbar, [dim_v,time_name], /FLOAT
        obry->AttAdd, var_vbar, 'valid_range', uv_range
      endelse
      obry->AttAdd, var_vbar, 'time', time_name

    endif

  endfor

  ;; For each time, retrieve velocity and zeta if available)

  for r=0,n_time-1 do begin

    message, /INFORM, string(FORMaT='(%"Reading velocity data for record %d")', r)

    for b=0,3 do begin

      if has_uv[b] then begin

        var_u = 'u_'+bname[b]
        var_v = 'v_'+bname[b]

        var_ubar = 'ubar_'+bname[b]
        var_vbar = 'vbar_'+bname[b]

        ub = obry->VarGet(var_u, OFFSET=[0,0,r], COUNT=[0,0,1], /AUTOSCALE)
        vb = obry->VarGet(var_v, OFFSET=[0,0,r], COUNT=[0,0,1], /AUTOSCALE)

        dim_ub = size(ub, /DIMENSIONS)
        dim_vb = size(vb, /DIMENSIONS)

        if has_zeta[b] then begin
          var_zeta = 'zeta_'+bname[b]
          zeta = obry->VarGet(var_zeta, OFFSET=[0,r], COUNT=[0,1], /AUTOSCALE)
          ;; Note that assigning boundary zeta values to normal velocity points
          ;; implies a half-grid-cell shift in space.
          case b mod 2 of
            0: begin
              zeta_u = mgh_stagger(zeta, DELTA=[-1])
              zeta_v = zeta
            end
            1: begin
              zeta_u = zeta
              zeta_v = mgh_stagger(zeta, DELTA=[-1])
            end
          endcase
          mgh_undefine, zeta
        endif else begin
          zeta_u = fltarr(dim_ub[0])
          zeta_v = fltarr(dim_vb[0])
        endelse

        message, /INFORM, string(FORMaT='(%"Calculating ubar & vbar for record %d")', r)

        ubar = fltarr(dim_ub[0])
        vbar = fltarr(dim_vb[0])

        ;; We need depths at the boundary velocity points.

        case b of
          0: begin
            hb_u = h_u[*,0]
            hb_v = h_v[*,0]
          end
          1: begin
            hb_u = reform(h_u[-1,*])
            hb_v = reform(h_v[-1,*])
          end
          2: begin
            hb_u = h_u[*,-1]
            hb_v = h_v[*,-1]
          end
          3: begin
            hb_u = reform(h_u[0,*])
            hb_v = reform(h_v[0,*])
          end
        endcase

        for i=0,dim_ub[0]-1 do begin
          zz = mgh_roms_s_to_z(s_w, hb_u[i], ZETA=zeta_u[i], CS=cs, HC=hc, VTRANSFORM=vtransform)
          dz = mgh_diff(zz)
          ubar[i] = total(dz*reform(ub[i,*]))/total(dz)
          mgh_undefine, zz, dz
        endfor

        for i=0,dim_vb[0]-1 do begin
          zz = mgh_roms_s_to_z(s_w, hb_v[i], ZETA=zeta_v[i], CS=cs, HC=hc, VTRANSFORM=vtransform)
          dz = mgh_diff(zz)
          vbar[i] = total(dz*reform(vb[i,*]))/total(dz)
          mgh_undefine, zz, dz
        endfor

        if keyword_set(pack_data) then begin
          ubar = round((temporary(ubar)-uv_scale[0])/uv_scale[1])
          vbar = round((temporary(vbar)-uv_scale[0])/uv_scale[1])
        endif else begin
          ubar = temporary(ubar) > uv_range[0] < uv_range[1]
          vbar = temporary(vbar) > uv_range[0] < uv_range[1]
        endelse

        message, /INFORM, string(FORMaT='(%"Writing ubar & vbar for record %d")', r)

        obry->VarPut, var_ubar, temporary(ubar), OFFSET=[0,r]
        obry->VarPut, var_vbar, temporary(vbar), OFFSET=[0,r]

      endif

    endfor

  endfor

  ;; Close the output file

  obj_destroy, obry

end

