;+
; NAME:
;   MGH_ROMS_SCALAR_HARMONICS_TO_NCDF
;
; PURPOSE:
;   Get scalar data from a ROMS file, perform an annual-harmonic analysis
;   and save the results to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   mgh_roms_scalar_harmonics_to_ncdf, history, file_out
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string array
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   file_out (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
; KEYWORD PARAMETERS:
;   DEPTH
;     Set this keyword to a numeric vector specifying the depth(s) of z
;     surface(s) on which data are to be analysed. This keyword should
;     be specified only for variables having a vertical dimension and
;     it cannot be used together with LEVEL or SIGMA.
;
;   LEVEL
;     Set this keyword to an integer vector to specify the
;     s-coordinate level(s) on which data are to be analysed.  This
;     keyword should be specified only for variables having a vertical
;     dimension and it cannot be used together with DEPTH or SIGMA.
;
;   SIGMA
;     Set this keyword to a numeric vector specifying the sigma
;     surface(s) on which data are to be analysed. This keyword should
;     be specified only for variables having a vertical dimension and
;     it cannot be used together with DEPTH or LEVEL.
;
;   TIME_RANGE (input, numeric 2-element vector)
;     Time interval (in days) over which to perform the tidal analysis.
;     Default is [0,tmax], where tmin and tmax are minimum and
;     maximum times.
;
;   VARIABLE (input, string scalar)
;     The name of a scalar variable to be analysed. Default is 'zeta'.
;
; TO DO:
;   The square tiling stuff that was introduced for efficiency reasons is
;   very intricate (but now tested). Perhaps it would be better to
;   work in rows?
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-04:
;     Written.
;-
pro mgh_roms_scalar_harmonics_to_ncdf, history, file_out, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     CYCLE_LENGTH=cycle_length, N_HARMONICS=n_harmonics, $
     TIME_RANGE=time_range, VARIABLE=variable

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         history_file = ohis
      end
      'OBJREF': begin
         ohis = history
         history_file = history
         history_destroy = 0
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
        message, 'Name for out file not supplied'

   ;; Process keywords

   if n_elements(variable) eq 0 then $
        variable = 'zeta'

   if n_elements(n_harmonics) eq 0 then $
        n_harmonics = 4

   if n_elements(cycle_length) eq 0 then $
        cycle_length = 365.2422D0

   ;; Check out the dimensionality of the variables

   dim = ohis->VarDims(variable)

   if min([strlen(dim.horizontal),strlen(dim.time)]) eq 0 then $
        message, 'The variable must have horizontal & time dimensions'
   if ~ array_equal(dim.horizontal, ['xi_rho','eta_rho']) then $
        message, 'I can handle only variables on the horizontal RHO grid right now.'

   ;; Process keywords relating to the vertical grid

   dim_vertical = dim.vertical
   has_vertical = strlen(dim_vertical) gt 0

   if has_vertical then begin
      if ohis->HasVar(dim_vertical) then begin
         s = ohis->VarGet(dim_vertical)
      endif else begin
         n_s_rho = ohis->DimInfo('s_rho', /DIMSIZE)
         s = mgh_range(-1, 0, STRIDE=1.0D/n_s_rho)
         if dim_vertical eq 's_rho' then $
	      s = mgh_stagger(s, DELTA=-1)
      endelse
      n_key = (n_elements(depth) gt 0) + (n_elements(levels) gt 0) + $
              (n_elements(sigmas) gt 0)
      if n_key gt 1 then $
           message, 'The DEPTH, LEVEL & SIGMA keywords cannot be used together'
      if n_key eq 0 then $
           level = lindgen(n_elements(s))
      mgh_undefine, n_key
      ;; By now we have exactly one of the depth-related keywords defined.
      case 1B of
         n_elements(level) gt 0: n_vertical = n_elements(level)
         n_elements(depth) gt 0: n_vertical = n_elements(depth)
         n_elements(sigma) gt 0: n_vertical = n_elements(sigma)
      endcase
   endif else begin
      fmt = '(%"The %s keyword is not required or allowed when the ' + $
            'variable has no vertical ension")'
      if n_elements(level) gt 0 then message, string(FORMAT=fmt, 'LEVEL')
      if n_elements(depth) gt 0 then message, string(FORMAT=fmt, 'DEPTH')
      if n_elements(sigma) gt 0 then message, string(FORMAT=fmt, 'SIGMA')
   endelse

   ;; Get time data from the history file and process the TIME_RANGE parameter.

   dim_time = dim.time

   case 1B of
      ohis->HasVar(dim_time): $
           time_var = dim_time
      ohis->HasVar('ocean_time'): $
           time_var = 'ocean_time'
      ohis->HasVar('scrum_time'): $
           time_var = 'scrum_time'
      else: $
           message, 'Time variable not found'
   endcase

   time = ohis->VarGet(time_var, AUTOSCALE=0)

   time_units = {scale: 1D/(24D*3600D), offset: 0D}
   if ohis->MGHncSequence::HasAtt(time_var, 'units') then $
        time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))

   time = time*time_units.scale

   time_ref = time_units.offset

   if n_elements(time_range) eq 0 then time_range = mgh_minmax(time) + [3,0]

   rrange = mgh_subset(time, time_range+[-1,1]*1.0D-6)

   rr0 = rrange[0]
   rr1 = rrange[1]
   rrn = rr1 - rr0 + 1

   time = time[rr0:rr1]

   ;; Does the history file have a mask?

   has_mask = ohis->HasVar('mask_rho')

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   onc_out = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   ;; A bit of header info

   fmt = '(%"Scalar harmonics from ROMS dataset %s")'
   onc_out->AttAdd, /GLOBAL, 'long_name', $
        string(mgh_get_property(ohis, /NAME), FORMAT=fmt)

   onc_out->AttAdd, /GLOBAL, 'time_range', time_range
   onc_out->AttAdd, /GLOBAL, 'cycle_length', cycle_length

   ;; Create netCDF dimensions. The horizontal grid is the ROMS
   ;; interior rho grid

   hdim = [ohis->DimInfo('xi_rho', /DIMSIZE)-2, $
           ohis->DimInfo('eta_rho', /DIMSIZE)-2]

   onc_out->DimAdd, 'xi', hdim[0]
   onc_out->DimAdd, 'eta', hdim[1]

   if has_vertical then begin
      case 1B of
         n_elements(level) gt 0: dim_vertical = 's'
         n_elements(depth) gt 0: dim_vertical = 'depth'
         n_elements(sigma) gt 0: dim_vertical = 'sigma'
      endcase
      onc_out->DimAdd, dim_vertical, n_vertical
   endif

   onc_out->DimAdd, 'harmonic', 2*n_harmonics+1

   ;; Create variables

   onc_out->VarAdd, 'lon', ['xi','eta'], /FLOAT
   onc_out->AttAdd, 'lon', 'long_name', 'longitude'
   onc_out->AttAdd, 'lon', 'units', 'degree_east'

   onc_out->VarAdd, 'lat', ['xi','eta'], /FLOAT
   onc_out->AttAdd, 'lat', 'long_name', 'latitude'
   onc_out->AttAdd, 'lat', 'units', 'degree_north'

   if has_mask then begin
      onc_out->VarAdd, 'mask', ['xi','eta'], /BYTE
      onc_out->AttAdd, 'mask', 'long_name', 'mask'
      onc_out->AttAdd, 'mask', 'option_0', 'land'
      onc_out->AttAdd, 'mask', 'option_1', 'water'
   endif

   onc_out->VarAdd, 'h', ['xi','eta'], /FLOAT
   onc_out->AttAdd, 'h', 'long_name', 'bathymetry'
   onc_out->AttAdd, 'h', 'units', 'meter'

   if has_vertical then begin
      onc_out->VarAdd, dim_vertical, dim_vertical
   endif

   ;; Create netCDF variable for results

   coeff_name = mgh_str_vanilla(variable)+'_coefficients'
   stdev_name = mgh_str_vanilla(variable)+'_residual_sd'

   fill_real = 1.E10

   dim = ['xi','eta']
   if has_vertical then dim = [dim,dim_vertical]
   coeff_dim = [dim,'harmonic']
   stdev_dim = dim

   onc_out->VarAdd, coeff_name, coeff_dim, /FLOAT
   onc_out->AttAdd, coeff_name, '_FillValue', fill_real

   onc_out->VarAdd, stdev_name, stdev_dim, /FLOAT
   onc_out->AttAdd, stdev_name, '_FillValue', fill_real

   ;; Add parameter & grid data

   lon = ohis->VarGet('lon_rho', OFFSET=[1,1], COUNT=hdim)
   onc_out->VarPut, 'lon', lon

   lat = ohis->VarGet('lat_rho', OFFSET=[1,1], COUNT=hdim)
   onc_out->VarPut, 'lat', lat

   h = ohis->VarGet('h', OFFSET=[1,1], COUNT=hdim)
   onc_out->VarPut, 'h', h

   if has_mask then begin
      mask = round(ohis->VarGet('mask_rho', OFFSET=[1,1], COUNT=hdim))
      onc_out->VarPut, 'mask', mask
   endif

   if has_vertical then begin
      case 1B of
         n_elements(level) gt 0: $
              onc_out->VarPut, dim_vertical, s[level]
         n_elements(depth) gt 0: $
              onc_out->VarPut, dim_vertical, depth
         n_elements(sigma) gt 0: $
              onc_out->VarPut, dim_vertical, sigma
      endcase
   endif

   ;; Set up result array.

   dim = [hdim]
   if has_vertical then dim = [dim,n_vertical]
   coeff = make_array([dim,2*n_harmonics+1], VALUE=fill_real)
   stdev = make_array(dim, VALUE=fill_real)

   ;; Retrieve data and calculate coefficients. A tiling scheme is
   ;; used for retrieving data to speed up access while avoiding large
   ;; memory requirements. Tile offsets are expressed relative to the
   ;; SW interior rho point.

   n_tile = round(hdim/128.D0) > 2

   tile_size = floor(hdim/double(n_tile))

   tile_offset_x = tile_size[0]*lindgen(n_tile[0])
   tile_offset_y = tile_size[1]*lindgen(n_tile[1])

   mgh_undefine, tile_size

   tile_count_x = [tile_offset_x[1:*],hdim[0]] - tile_offset_x
   tile_count_y = [tile_offset_y[1:*],hdim[1]] - tile_offset_y

   for l=0,n_tile[1]-1 do begin
      for k=0,n_tile[0]-1 do begin

         fmt = '(%"Processing tile %d, %d, offset %d, %d...")'
         print, FORMAT=fmt, k, l, tile_offset_x[k], tile_offset_y[l]

         xi_range = 1 + tile_offset_x[k] + [0,tile_count_x[k]-1]
         eta_range = 1 + tile_offset_y[l] + [0,tile_count_y[l]-1]

         grid = ohis->HsliceGrid(variable, $
                                 XI_RANGE=xi_range, ETA_RANGE=eta_range)

         ;; Set up an array to hold data retrieved for this tile

         dim = [tile_count_x[k],tile_count_y[l]]
         if has_vertical then dim = [dim,n_vertical]
         dim = [dim,rrn]

         data = make_array(dim, VALUE=0.)

         ;; Accumulate tile data

         print, 'Accumulating tile data'

         for r=0,rrn-1 do begin

            ss = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=0, $
                                  RECORD=rr0+r, DEPTH=depth, LEVEL=level, $
                                  SIGMA=sigma)

            if has_vertical then begin
               data[*,*,*,r] = temporary(ss)
            endif else begin
               data[*,*,r] = temporary(ss)
            endelse

         endfor

         ;; Process tile data, working on the time series at each
         ;; point

         print, 'Processing tile data'

         for j=0,tile_count_y[l]-1 do begin
            for i=0,tile_count_x[k]-1 do begin

               ii = tile_offset_x[k] + i
               jj = tile_offset_y[l] + j

               if (~ has_mask) || mask[ii,jj] then begin

                  w = replicate(1, rrn)

                  if has_vertical then begin
                     for m=0,n_vertical-1 do begin
                        s = reform(data[i,j,m,*])
			                  l_good = where(finite(s), n_good)
			                  if n_good ge rrn/2 then begin
                          ann = fltarr(2*n_harmonics+1)
                          sfit = mgh_curvefit(time, s, w, ann, $
                                              FUNCTION_NAME='mgh_fit_annual', $
                                              N_HARMONICS=n_harmonics, $
                                              CYCLE_LENGTH=cycle_length)
                          coeff[ii,jj,m,*] = ann
                          stdev[ii,jj,m] = mgh_stdev(s-sfit)
			endif
                     endfor
                  endif else begin
                     s = reform(data[i,j,*])
                     l_good = where(finite(s), n_good)
                     if n_good ge rrn/2 then begin
                       ann = fltarr(2*n_harm+1)
                       sfit = mgh_curvefit(time, s, w, ann, $
                                           FUNCTION_NAME='mgh_fit_annual', $
                                           N_HARMONICS=n_harmonics, $
                                           CYCLE_LENGTH=cycle_length)
                     endif
                     coeff[ii,jj,*] = temporary(ann)
                     stdev[ii,jj,m] = mgh_stdev(s-sfit)
                  endelse

               endif

            endfor
         endfor

      endfor
   endfor

   ;; Write data & clean up

   onc_out->VarPut, coeff_name, coeff
   onc_out->VarPut, stdev_name, stdev

   obj_destroy, onc_out

end
