;+
; NAME:
;   MGH_ROMS_BRY_LOAD_ROMS_REGRID
;
; PURPOSE:
;   Into a ROMS boundary file, load regridded data from a ROMS
;   history-file sequence or similar.
;
; CALLING SEQUENCE:
;   mgh_roms_bry_load_roms_regrid, file_bry, history
;
; POSITIONAL PARAMETERS:
;   file_bry (input, string scalar)
;     The name of a ROMS boundary file into which the data are to
;     be written. This file must contain grid information: we'll call
;     this grid the "inner grid".
;
;   history (input, string scalar/vector OR object reference)
;     The ROMS history-file sequence from which data are to be read, either
;     as a list of file names or an MGHromsHistory object. This must contain grid
;     information: we'll call this grid the "outer grid".
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2006-11:
;     Written.
;   Mark Hadfield, 2013-10:
;     - Ported functionality from MGH_ROMS_CLM_LOAD_ROMS_REGRID.
;     - Now reads the necessary grid data from the boundary files.
;   Mark Hadfield, 2015-10:
;     - Reformatted.
;-
pro mgh_roms_bry_load_roms_regrid_var2d_fill, values, TRIANGLES=triangles

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Fill 2D scalar data on the outer grid. The triangles keyowrd allows
   ;; triangle data to be stored between invocations.

   values = mgh_fill2d(values, METHOD=1, TRIANGLES=triangles)

end

function mgh_roms_bry_load_roms_regrid_var2d_boundary, values, rloc

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Interpolate 2D scalar data from the outer grid to the boundary slice

   return, interpolate(values, reform(rloc[0,*]), reform(rloc[1,*]))

end

function mgh_roms_bry_load_roms_regrid_var3d_intermediate, values, zloc_on

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE
   compile_opt HIDDEN

   ;; Interpolate 3D scalar data from the outer grid to the intermediate grid

   dim = size(zloc_on, /DIMENSIONS)

   ii = rebin(reform(findgen(dim[0]), [dim[0],1,1]), dim)
   jj = rebin(reform(findgen(dim[1]), [1,dim[1],1]), dim)

   return, interpolate(values, temporary(ii), temporary(jj), zloc_on, MISSING=!values.f_nan)

end

pro mgh_roms_bry_load_roms_regrid_var3d_fill, values

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE
   compile_opt HIDDEN

   ;; Fill 3D scalar data (on the intermediate grid). We use
   ;; the default method here (fill along latitude circles first, etc) because
   ;; it is *much* faster than the alternatives. It will give biased results where
   ;; the outer grid is not aligned east-west.

   dim = size(values, /DIMENSIONS)

   for k=0,dim[2]-1 do values[*,*,k] = mgh_fill2d(values[*,*,k], METHOD=0)

end

function mgh_roms_bry_load_roms_regrid_var3d_boundary, values, rloc, zloc_ni

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE
   compile_opt HIDDEN

   ;; Interpolate 3D scalar data (on the intermediate grid) to the boundary locations

   dim = size(zloc_ni, /DIMENSIONS)

   ii = rebin(reform(rloc[0,*], [dim[0],1]), dim)
   jj = rebin(reform(rloc[1,*], [dim[0],1]), dim)

   return, interpolate(values, temporary(ii), temporary(jj), zloc_ni)

end

pro mgh_roms_bry_load_roms_regrid_show_indices, rloc, r2loc

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   for b=0,3 do begin
      loc = *rloc[b]
      if b eq 0 then begin
         xrange = [min(loc[0,*]),max(loc[0,*])]
         yrange = [min(loc[1,*]),max(loc[1,*])]
      endif else begin
         xrange = [min(loc[0,*]) < xrange[0],max(loc[0,*]) > xrange[1]]
         yrange = [min(loc[1,*]) < yrange[0],max(loc[1,*]) > yrange[1]]
      endelse
      mgh_undefine, loc
   endfor

   !null = plot([0], /NODATA, XRANGE=xrange+[-1,1], YRANGE=yrange+[-1,1])
   for b=0,3 do begin
      loc = *r2loc[b]
      case (b mod 2) of
         0: begin
            !null = plot(loc[0,*,0], loc[1,*,0], /OVER, COLOR=mgh_color('blue'))
            !null = plot(loc[0,*,1], loc[1,*,1], /OVER, COLOR=mgh_color('red'))
         end
         1: begin
            !null = plot(loc[0,0,*], loc[1,0,*], /OVER, COLOR=mgh_color('blue'))
            !null = plot(loc[0,1,*], loc[1,1,*], /OVER, COLOR=mgh_color('red'))
         end
      endcase
   endfor

end

pro mgh_roms_bry_load_roms_regrid, file_bry, history, $
     AUTOSCALE_TIME=autoscale_time, PACK_DATA=pack_data, $
     N_DYE=n_dye, N_MUD=n_mud, N_SAND=n_sand, $
     RECORDS=records, VERTICAL=vertical, $
     TIME_NAME=time_name, TIME_SHIFT=time_shift

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   cj = complex(0, 1)

   ;; Check arguments

   if n_elements(file_bry) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file_bry'

   if ~ file_test(file_bry, /READ) then $
      message, 'Boundary file cannot be read'

   if ~ file_test(file_bry, /WRITE) then $
      message, 'Boundary file cannot be written to'

   if n_elements(history) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'

   if n_elements(time_name) eq 0 then time_name = 'ocean_time'

   ;; Open the files

   fmt = '(%"Opening boundary file %s for modification")'
   message, /INFORM, string(FORMAT=fmt, file_bry)

   obry = obj_new('MGHncFile', file_bry, /MODIFY)

   case 1B of
      isa(history, 'STRING'): begin
         fmt = '(%"Opening history file sequence beginning with %s")'
         message, /INFORM, string(FORMAT=fmt, history[0])
         ohis = obj_new('MGHromsHistory', history)
         history_destroy = 1B
      end
      isa(history, 'OBJREF'): begin
         ohis = history
      end
      else: $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
   endcase

   ;; Establish records to be read

   if n_elements(records) eq 0 then begin
      ohis->GetProperty, N_RECORDS=n_rec
      records = lindgen(n_rec)
   endif

   if n_elements(time_shift) eq 0 then $
      time_shift = 0

   ;; If number of dye, mud or sand variables is unspecified, determine
   ;; them from the history file

   if n_elements(n_dye) eq 0 then begin
      n_dye = 0
      while 1B do begin
         var = string(FOrMaT='(%"dye_%2.2d")', n_dye+1)
         if ~ ohis->HasVar(var) then break
         n_dye += 1
      endwhile
      mgh_undefine, var
   endif

   if n_elements(n_mud) eq 0 then begin
      n_mud = 0
      while 1B do begin
         var = string(FOrMaT='(%"mud_%2.2d")', n_mud+1)
         if ~ ohis->HasVar(var) then break
         n_mud += 1
      endwhile
      mgh_undefine, var
   endif

   if n_elements(n_sand) eq 0 then begin
      n_sand = 0
      while 1B do begin
         var = string(FOrMaT='(%"sand_%2.2d")', n_sand+1)
         if ~ ohis->HasVar(var) then break
         n_sand += 1
      endwhile
      mgh_undefine, var
   endif

   ;; Check if boundary file has a vertical dimension

   if n_elements(vertical) eq 0 then vertical = obry->HasDim('s_rho')

   ;; Specify valid ranges for variables (and scaling parameters if applicable)
   ;; for the netCDF variables.

   zeta_range = [-10.0,10.0]
   uv_range = [-10.0,10.0]
   if keyword_set(vertical) then begin
      temp_range = [-5.0,45.0]
      salt_range = [0.0,40.0]
      dye_range = [-0.5,1.5]
      mud_range = [0.0,10.0]
      sand_range = [0.0,10.0]
   endif

   if keyword_set(pack_data) then begin
      pack_range = [-32766S,32767S]
      zeta_scale = mgh_norm_coord(pack_range, zeta_range)
      uv_scale = mgh_norm_coord(pack_range, uv_range)
      if keyword_set(vertical) then begin
         temp_scale = mgh_norm_coord(pack_range, temp_range)
         salt_scale = mgh_norm_coord(pack_range, salt_range)
         dye_scale = mgh_norm_coord(pack_range, dye_range)
         mud_scale = mgh_norm_coord(pack_range, mud_range)
         sand_scale = mgh_norm_coord(pack_range, sand_range)
      endif
   endif

   ;; Get grid data

   message, /INFORM, 'Getting inner and outer grid data'

   grdi_dim = [obry->DimInfo('xi_rho', /DIMSIZE), $
      obry->DimInfo('eta_rho', /DIMSIZE), $
      obry->DimInfo('s_rho', /DIMSIZE)]
   grdi_lon = obry->VarGet('lon_rho')
   grdi_lat = obry->VarGet('lat_rho')
   grdi_mask = obry->VarGet('mask_rho')
   grdi_angle = obry->VarGet('angle')

   if vertical then begin

      grdi_h = obry->VarGet('h')
      grdi_s_rho = obry->VarGet('s_rho')
      grdi_theta_s = obry->VarGet('theta_s')
      grdi_theta_b = obry->VarGet('theta_b')
      grdi_hc = obry->VarGet('hc')
      grdi_vstretch = obry->HasVar('Vstretching') ? obry->VarGet('Vstretching') : 1
      grdi_vtransform = obry->HasVar('Vtransform') ? obry->VarGet('Vtransform') : 1
      grdi_cs = mgh_roms_s_to_cs(grdi_s_rho, THETA_S=grdi_theta_s, THETA_B=grdi_theta_b, VSTRETCH=grdi_vstretch)

   endif

   grdo_dim = [ohis->DimInfo('xi_rho', /DIMSIZE), $
      ohis->DimInfo('eta_rho', /DIMSIZE), $
      ohis->DimInfo('s_rho', /DIMSIZE)]
   grdo_lon = ohis->VarGet('lon_rho')
   grdo_lat = ohis->VarGet('lat_rho')
   grdo_mask = ohis->VarGet('mask_rho')
   grdo_angle = ohis->VarGet('angle')

   if vertical then begin

      grdo_h = ohis->VarGet('h')
      grdo_s_rho = ohis->VarGet('s_rho')
      grdo_theta_s = ohis->VarGet('theta_s')
      grdo_theta_b = ohis->VarGet('theta_b')
      grdo_hc = ohis->VarGet('hc')
      grdo_vstretch = ohis->HasVar('Vstretching') ? ohis->VarGet('Vstretching') : 1
      grdo_vtransform = ohis->HasVar('Vtransform') ? ohis->VarGet('Vtransform') : 1
      grdo_cs = mgh_roms_s_to_cs(grdo_s_rho, THETA_S=grdo_theta_s, THETA_B=grdo_theta_b, VSTRETCH=grdo_vstretch)

   endif

   ;; Boundary names

   bname = ['south','east','north','west']

   ;; Add boundary variables

   ;; ...zeta

   for b=0,3 do begin
      var = 'zeta_'+bname[b]
      case b mod 2 of
         0: hdim = 'xi_rho'
         1: hdim = 'eta_rho'
      endcase
      message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
      if keyword_set(pack_data) then begin
         obry->VarAdd, var, [hdim,time_name], /SHORT
         obry->AttAdd, var, 'add_offset', zeta_scale[0]
         obry->AttAdd, var, 'scale_factor', zeta_scale[1]
         obry->AttAdd, var, 'valid_range', pack_range
      endif else begin
         obry->VarAdd, var, [hdim,time_name], /FLOAT
         obry->AttAdd, var, 'valid_range', zeta_range
      endelse
      obry->AttAdd, var, 'time', time_name
   endfor


   if keyword_set(vertical) then begin

      ;; ...temperature

      for b=0,3 do begin
         var = 'temp_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_rho'
            1: hdim = 'eta_rho'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /SHORT
            obry->AttAdd, var, 'add_offset', temp_scale[0]
            obry->AttAdd, var, 'scale_factor', temp_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', temp_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

      ;; ...salinity

      for b=0,3 do begin
         var = 'salt_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_rho'
            1: hdim = 'eta_rho'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /SHORT
            obry->AttAdd, var, 'add_offset', salt_scale[0]
            obry->AttAdd, var, 'scale_factor', salt_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', salt_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

      ;; ...u

      for b=0,3 do begin
         var = 'u_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_u'
            1: hdim = 'eta_u'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /SHORT
            obry->AttAdd, var, 'add_offset', uv_scale[0]
            obry->AttAdd, var, 'scale_factor', uv_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', uv_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

      ;; ...v

      for b=0,3 do begin
         var = 'v_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_v'
            1: hdim = 'eta_v'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /SHORT
            obry->AttAdd, var, 'add_offset', uv_scale[0]
            obry->AttAdd, var, 'scale_factor', uv_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,'s_rho',time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', uv_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

   endif else begin

      ;; ...ubar

      for b=0,3 do begin
         var = 'ubar_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_u'
            1: hdim = 'eta_u'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,time_name], /SHORT
            obry->AttAdd, var, 'add_offset', uv_scale[0]
            obry->AttAdd, var, 'scale_factor', uv_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', uv_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

      ;; ...vbar

      for b=0,3 do begin
         var = 'vbar_'+bname[b]
         case b mod 2 of
            0: hdim = 'xi_v'
            1: hdim = 'eta_v'
         endcase
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', var)
         if keyword_set(pack_data) then begin
            obry->VarAdd, var, [hdim,time_name], /SHORT
            obry->AttAdd, var, 'add_offset', uv_scale[0]
            obry->AttAdd, var, 'scale_factor', uv_scale[1]
            obry->AttAdd, var, 'valid_range', pack_range
         endif else begin
            obry->VarAdd, var, [hdim,time_name], /FLOAT
            obry->AttAdd, var, 'valid_range', uv_range
         endelse
         obry->AttAdd, var, 'time', time_name
      endfor

   endelse

   ;; To interpolate scalars, we need the location of the boundary rho
   ;; points in the outer rho grid. The location arrays (below) are
   ;; dimensioned [2,n_xi_rho] for south & north and [2,n_eta_rho] for east and
   ;; west.

   rloc = ptrarr(4)

   for b=0,3 do begin
      case b of
         0: begin
            loni = grdi_lon[*,0]
            lati = grdi_lat[*,0]
         end
         1: begin
            loni = reform(grdi_lon[-1,*])
            lati = reform(grdi_lat[-1,*])
         end
         2: begin
            loni = grdi_lon[*,-1]
            lati = grdi_lat[*,-1]
         end
         3: begin
            loni = reform(grdi_lon[0,*])
            lati = reform(grdi_lat[0,*])
         end
      endcase
      rloc[b] = ptr_new(mgh_locate2(grdo_lon, grdo_lat, XOUT=loni, YOUT=lati, /DOUBLE))
      mgh_undefine, loni, lati
   endfor

   ;; To handle velocities we need the location of the
   ;; two outermost rows/columns of rho points, i.e. the boundary
   ;; points plus the ones inside that, in the outer rho grid.  The
   ;; location arrays (below) are dimensioned [2,n_xi_rho,2] for south &
   ;; north and [2,2,n_eta_rho] for east and west. The extra dimension is
   ;; for the extra row/column and it is reversed on the east and
   ;; north boundaries so index 0 stands for the outer row/column and
   ;; index 1 stands for the next one in.

   r2loc = ptrarr(4)

   for b=0,3 do begin
      case b of
         0: begin
            loni = grdi_lon[*,0:1]
            lati = grdi_lat[*,0:1]
         end
         1: begin
            loni = reverse(grdi_lon[-2:-1,*], 1)
            lati = reverse(grdi_lat[-2:-1,*], 1)
         end
         2: begin
            loni = reverse(grdi_lon[*,-2:-1], 2)
            lati = reverse(grdi_lat[*,-2:-1], 2)
         end
         3: begin
            loni = grdi_lon[0:1,*]
            lati = grdi_lat[0:1,*]
         end
      endcase
      r2loc[b] = ptr_new(mgh_locate2(grdo_lon, grdo_lat, XOUT=loni, YOUT=lati, /DOUBLE))
      mgh_undefine, loni, lati
   endfor

   ;; For 3D *scalar* interpolation, we will use an intermediate
   ;; z-level grid, horizontally the same as the outer grid
   ;; and spanning the depth range of the inner grid. This is done
   ;; to accommodate large differences between the inner and outer grid
   ;; topography: filling of the outer grid data is done on the horizontal
   ;; surfaces of the intermediate grid. The relationship between these grids
   ;; is specified by variables zloc_on (the vertical position of intermediate grid
   ;; points on the outer grid) and zloc_ni (the vertical position of
   ;; inner grid points on the intermediate grid).

   ;; For 3D *velocity* interpolation, we will go straight from the outer
   ;; to the inner grid, with variable zloc (the vertical position of inner grid
   ;; points on the outer grid).

   if keyword_set(vertical) then begin

      message, /INFORM, 'Vertical grid location'

      ;; Outer grid heights.

      grdo_z = dblarr(grdo_dim)
      for i=0,grdo_dim[0]-1 do begin
         for j=0,grdo_dim[1]-1 do begin
            grdo_z[i,j,*] = mgh_roms_s_to_z(grdo_s_rho, grdo_h[i,j], HC=grdo_hc, CS=grdo_cs, VTRANSFORM=grdo_vtransform)
         endfor
      endfor

      ;; For 3D velocity interpolation, locate boundary grid points
      ;; (as specified by r2loc) on the outer grid. Clip zloc to the
      ;; vertical extent of the outer grid.

      z2loc = ptrarr(4)

      for b=0,3 do begin
         case b of
            0: h = grdi_h[*,0:1]
            1: h = reverse(grdi_h[-2:-1,*], 1)
            2: h = reverse(grdi_h[*,-2:-1], 2)
            3: h = grdi_h[0:1,*]
         endcase
         loc = *r2loc[b]
         hh = interpolate(grdo_h, reform(loc[0,*,*]), reform(loc[1,*,*]))
         dim = size(h, /DIMENSIONS)
         ll = fltarr([dim,grdi_dim[2]])
         for j=0,dim[1]-1 do begin
            for i=0,dim[0]-1 do begin
               zz = mgh_roms_s_to_z(grdo_s_rho, hh[i,j], HC=grdo_hc, CS=grdo_cs, VTRANSFORM=grdo_vtransform)
               z = mgh_roms_s_to_z(grdi_s_rho, h[i,j], HC=grdi_hc, CS=grdi_cs, VTRANSFORM=grdi_vtransform)
               ll[i,j,*] = mgh_locate(zz, XOUT=z, /EXTRAPOLATE)
            endfor
         endfor
         ll = ll > 0 < (grdo_dim[2]-1)
         z2loc[b] = ptr_new(ll)
      endfor

      ;; For 3D scalar interpolation, locate intermediate grid points
      ;; on the outer grid.

      ss = mgh_range(-1D, 0D, N_ELEMENTS=2*grdi_dim[2]+1)

      zz = mgh_roms_s_to_z(ss, max(grdi_h), HC=grdi_hc, THETA_S=grdi_theta_s, THETA_B=grdi_theta_b, $
         VSTRETCH=grdi_vstretch, VTRANSFORM=grdi_vtransform)

      zloc_on = dblarr([grdo_dim[0:1],n_elements(zz)])
      for i=0,grdo_dim[0]-1 do begin
         for j=0,grdo_dim[1]-1 do begin
            ;        z = mgh_roms_s_to_z(grdo_s_rho, grdo_h[i,j], HC=grdo_hc, CS=grdo_cs, VTRANSFORM=grdo_vtransform)
            zloc_on[i,j,*] = mgh_locate(reform(grdo_z[i,j,*]), XOUT=zz, /EXTRAPOLATE)
         endfor
      endfor

      ;; Intermediate grid points that have been located above the outer
      ;; grid are "pushed down". Intermediate grid points that have been
      ;; located below the outer grid are set to missing.

      zloc_on =  zloc_on < (grdo_dim[2]-1)
      l_neg = where(zloc_on lt 0, n_neg)
      if n_neg gt 0 then zloc_on[l_neg] = !values.f_nan

      ;; Now locate boundary grid points (as specified by rloc) on the
      ;; intermediate grid.

      zloc_ni = ptrarr(4)

      for b=0,3 do begin
         case b of
            0: h = grdi_h[*,0]
            1: h = reform(grdi_h[-1,*])
            2: h = grdi_h[*,-1]
            3: h = reform(grdi_h[0,*])
         endcase
         zl1 = fltarr(n_elements(h), grdi_dim[2])
         for i=0,n_elements(h)-1 do begin
            z = mgh_roms_s_to_z(grdi_s_rho, h[i], HC=grdi_hc, CS=grdi_cs, VTRANSFORM=grdi_vtransform)
            zl1[i,*] = mgh_locate(zz, XOUT=z, /EXTRAPOLATE)
         endfor
         zloc_ni[b] = ptr_new(zl1, /NO_COPY)
      endfor

   endif

   ;; Show plots of the indices

   if keyword_set(show_indices) then $
      mgh_roms_bry_load_roms_regrid_show_indices, rloc, r2loc

   ;; Calculate difference in angle between inner and outer grids
   ;; on boundary points.

   r2ang = ptrarr(4)

   for b=0,3 do begin

      ;; Inner grid angle
      case b of
         0: ai = grdi_angle[*,0:1]
         1: ai = reverse(grdi_angle[-2:-1,*], 1)
         2: ai = reverse(grdi_angle[*,-2:-1], 2)
         3: ai = grdi_angle[0:1,*]
      endcase

      ;; Outer grid angle
      loc = *r2loc[b]
      ao = interpolate(grdo_angle, reform(loc[0,*,*]), reform(loc[1,*,*]))

      ;; Difference
      r2ang[b] = ptr_new(temporary(ai) - temporary(ao))

   endfor

   ;; Work through records, reading, interpolating &
   ;; writing.

   for r=0,n_elements(records)-1 do begin

      fmt = '(%"Loading boundary data from record %d into record %d")'
      message, /INFORM, string(FORMAT=temporary(fmt), records[r], r)

      ;; Optionally convert time according to the attributes of the
      ;; input and output files

      time_his = ohis->VarGet(time_name, AUTOSCALE=0, OFFSET=records[r], COUNT=[1])

      if keyword_set(autoscale_time) then begin

         tu_his = {scale: 1D, offset: 0D}
         if ohis->HasAtt(time_name, 'units') then $
            tu_his = mgh_dt_units(ohis->AttGet(time_name, 'units'))

         tu_bry = {scale: 1D, offset: 0D}
         if obry->HasAtt(time_name, 'units') then $
            tu_bry = mgh_dt_units(obry->AttGet(time_name, 'units'))

         ;; Take the time offset in the history file into account if
         ;; and only if the time offset in the climatology file is
         ;; non-zero.. This is done because ROMS currently adds a
         ;; gratuitous reference time of "0001-01-01" when TIME_REF is
         ;; set to 0, leading to a tu_his.offset value of 1721423.5.

         if tu_bry.offset ne 0 then begin
            time_bry = (time_his*tu_his.scale+tu_his.offset-tu_bry.offset)/tu_bry.scale
         endif else begin
            time_bry = time_his*tu_his.scale/tu_bry.scale
         endelse

      endif else begin

         time_bry = time_his

      endelse

      case n_elements(time_shift) of
         0:
         1: time_bry += time_shift[0]
         n_elements(records): time_bry += time_shift[r]
      endcase

      fmt = '(%"Input time is %f; output time is %f")'
      message, /INFORM, string(FORMAT=temporary(fmt), time_his, time_bry)

      obry->VarPut, time_name, time_bry, OFFSET=[r]

      ;; Variables...

      ;; ...zeta

      zeta = ohis->VarGet('zeta', /AUTOSCALE, OFFSET=[0,0,records[r]], COUNT=[0,0,1])

      mgh_roms_bry_load_roms_regrid_var2d_fill, zeta, TRIANGLES=seta_triangles

      for b=0,3 do begin
         zb = mgh_roms_bry_load_roms_regrid_var2d_boundary(zeta, *rloc[b])
         if keyword_set(pack_data) then begin
            zb = round((temporary(zb)-zeta_scale[0])/zeta_scale[1])
            zb = temporary(zb) > pack_range[0] < pack_range[1]
         endif else begin
            zb = temporary(zb) > zeta_range[0] < zeta_range[1]
         endelse
         obry->VarPut, 'zeta_'+bname[b], OFFSET=[0,r], zb
         mgh_undefine, zb
      endfor

      mgh_undefine, zeta

      if keyword_set(vertical) then begin

         ;; ...temperature

         temp = ohis->VarGet('temp', /AUTOSCALE, OFFSET=[0,0,0,records[r]], COUNT=[0,0,0,1])

         temp = mgh_roms_bry_load_roms_regrid_var3d_intermediate(temp, zloc_on)

         mgh_roms_bry_load_roms_regrid_var3d_fill, temp

         for b=0,3 do begin
            tb = mgh_roms_bry_load_roms_regrid_var3d_boundary(temp, *rloc[b], *zloc_ni[b])
            if keyword_set(pack_data) then begin
               tb = round((temporary(tb)-temp_scale[0])/temp_scale[1])
               tb = temporary(tb) > pack_range[0] < pack_range[1]
            endif else begin
               tb = temporary(tb) > temp_range[0] < temp_range[1]
            endelse
            obry->VarPut, 'temp_'+bname[b], OFFSET=[0,0,r], tb
            mgh_undefine, tb
         endfor

         mgh_undefine, temp

         ;; ...salinity

         salt = ohis->VarGet('salt', /AUTOSCALE, OFFSET=[0,0,0,records[r]], COUNT=[0,0,0,1])

         salt = mgh_roms_bry_load_roms_regrid_var3d_intermediate(salt, zloc_on)

         mgh_roms_bry_load_roms_regrid_var3d_fill, salt

         for b=0,3 do begin
            sb = mgh_roms_bry_load_roms_regrid_var3d_boundary(salt, *rloc[b], *zloc_ni[b])
            if keyword_set(pack_data) then begin
               sb = round((temporary(sb)-salt_scale[0])/salt_scale[1])
               sb = temporary(sb) > pack_range[0] < pack_range[1]
            endif else begin
               sb = temporary(sb) > salt_range[0] < salt_range[1]
            endelse
            obry->VarPut, 'salt_'+bname[b], OFFSET=[0,0,r], sb
            mgh_undefine, sb
         endfor

         mgh_undefine, salt

         if n_dye+n_mud+n_sand gt 0 then message, 'Not yet!'

         ;; ...baroclinic velocities

         ;; ......get on outer u & v grids
         u = ohis->VarGet('u', OFFSET=[0,0,0,records[r]], COUNT=[0,0,0,1])
         v = ohis->VarGet('v', OFFSET=[0,0,0,records[r]], COUNT=[0,0,0,1])

         ;; ......fill
         l_miss = where(~ finite(u), n_miss)
         if n_miss gt 0 then u[l_miss] = 0
         l_miss = where(~ finite(v), n_miss)
         if n_miss gt 0 then v[l_miss] = 0

         ;; ......onto outer rho grid
         u = mgh_stagger(temporary(u), DELTA=[1,0,0])
         v = mgh_stagger(temporary(v), DELTA=[0,1,0])

         for b=0,3 do begin

            loc = *r2loc[b]
            zloc = *z2loc[b]
            dim = size(zloc, /DIMENSIONS)

            ;; ......interpolate onto boundary rho points
            xloc3d = rebin(reform(loc[0,*,*], [dim[0:1],1]), dim)
            yloc3d = rebin(reform(loc[1,*,*], [dim[0:1],1]), dim)
            ub = interpolate(u, xloc3d, yloc3d, zloc)
            vb = interpolate(v, xloc3d, yloc3d, zloc)
            mgh_undefine, xloc3d, yloc3d

            ;; ......rotate thru difference in angle between outer and inner grids
            ang3d = rebin(reform(*r2ang[b], [dim[0:1],1]), dim)
            uv = complex(temporary(ub), temporary(vb))*exp(-cj*ang3d)
            ub = real_part(uv)
            vb = imaginary(uv)
            mgh_undefine, uv, ang3d

            ;; ......onto boundary u & v points
            case (b mod 2) of
               0: begin
                  ub = mgh_stagger(reform(ub[*,0,*]), DELTA=[-1,0])
                  vb = mgh_avg(vb, 2)
               end
               1: begin
                  ub = mgh_avg(ub, 1)
                  vb = mgh_stagger(reform(vb[0,*,*]), DELTA=[-1,0])
               end
            endcase

            ;; ......pack data
            if keyword_set(pack_data) then begin
               ub = round((temporary(ub)-uv_scale[0])/uv_scale[1])
               vb = round((temporary(vb)-uv_scale[0])/uv_scale[1])
               ub = temporary(ub) > pack_range[0] < pack_range[1]
               vb = temporary(vb) > pack_range[0] < pack_range[1]
            endif else begin
               ub = temporary(ub) > uv_range[0] < uv_range[1]
               vb = temporary(vb) > uv_range[0] < uv_range[1]
            endelse

            obry->VarPut, 'u_'+bname[b], OFFSET=[0,0,r], ub
            obry->VarPut, 'v_'+bname[b], OFFSET=[0,0,r], vb

         endfor

      endif else begin

         ;; ...barotropic velocities

         ;; ......get on outer u & v grids
         ubar = ohis->VarGet('ubar', OFFSET=[0,0,records[r]], COUNT=[0,0,1])
         vbar = ohis->VarGet('vbar', OFFSET=[0,0,records[r]], COUNT=[0,0,1])

         ;; ......fill
         l_miss = where(~ finite(ubar), n_miss)
         if n_miss gt 0 then ubar[l_miss] = 0
         l_miss = where(~ finite(vbar), n_miss)
         if n_miss gt 0 then vbar[l_miss] = 0

         ;; ......onto outer rho grid
         ubar = mgh_stagger(temporary(ubar), DELTA=[1,0])
         vbar = mgh_stagger(temporary(vbar), DELTA=[0,1])

         for b=0,3 do begin

            loc = *r2loc[b]

            ;; ......interpolate onto boundary rho points
            ub = interpolate(ubar, reform(loc[0,*,*]), reform(loc[1,*,*]))
            vb = interpolate(vbar, reform(loc[0,*,*]), reform(loc[1,*,*]))

            ;; ......rotate thru difference in angle between outer and inner grids
            ang = *r2ang[b]
            uvb = complex(temporary(ub), temporary(vb))*exp(-cj*ang)
            ub = real_part(uvb)
            vb = imaginary(uvb)
            mgh_undefine, uvb, ang

            ;; ......onto boundary u & v points
            case (b mod 2) of
               0: begin
                  ub = mgh_stagger(ub[*,0], DELTA=-1)
                  vb = mgh_avg(vb, 2)
               end
               1: begin
                  ub = mgh_avg(ub, 1)
                  vb = mgh_stagger(reform(vb[0,*]), DELTA=-1)
               end
            endcase

            ;; ......pack data
            if keyword_set(pack_data) then begin
               ub = round((temporary(ub)-uv_scale[0])/uv_scale[1])
               vb = round((temporary(vb)-uv_scale[0])/uv_scale[1])
               ub = temporary(ub) > pack_range[0] < pack_range[1]
               vb = temporary(vb) > pack_range[0] < pack_range[1]
            endif else begin
               ub = temporary(ub) > uv_range[0] < uv_range[1]
               vb = temporary(vb) > uv_range[0] < uv_range[1]
            endelse

            obry->VarPut, 'ubar_'+bname[b], OFFSET=[0,r], ub
            obry->VarPut, 'vbar_'+bname[b], OFFSET=[0,r], vb

         endfor

      endelse

   endfor

   ;; Clean up

   ptr_free, rloc, r2loc, r2ang
   if vertical then ptr_free, zloc_ni

   message, /INFORM, 'Closing files'

   if keyword_set(history_destroy) then obj_destroy, ohis
   obj_destroy, obry

end

