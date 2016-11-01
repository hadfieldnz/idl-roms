;+
; NAME:
;   MGH_ROMS_TIDE_VECTOR_TO_NCDF
;
; PURPOSE:
;   Get vector data from a ROMS file, carry out a
;   time-stamped tidal analysis and save the results to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   mgh_roms_naive_vector_to_ncdf, history, file_out
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
;   CONSTITUENTS (input, string vector)
;     A list of tidal constituents. EAch one must appear in the NAME
;     field of the structure returned by MGH_T_GETCONSTS. Default is
;     ['M2'].
;
;   DEPTHS
;     Set this keyword to a numeric vector specifying the depth(s) of z
;     surface(s) on which data are to be analysed. This keyword should
;     be specified only for variables having a vertical dimension and
;     it cannot be used together with LEVELS or SIGMAS.
;
;   LEVELS
;     Set this keyword to an integer vector to specify the
;     s-coordinate level(s) on which data are to be analysed.  This
;     keyword should be specified only for variables having a vertical
;     dimension and it cannot be used together with DEPTHS or SIGMAS.
;
;   SIGMAS
;     Set this keyword to a numeric vector specifying the sigma
;     surface(s) on which data are to be analysed. This keyword should
;     be specified only for variables having a vertical dimension and
;     it cannot be used together with DEPTHS or SIGMAS.
;
;   TIME_RANGE (input, numeric 2-element vector)
;     Time interval (in days) over which to perform the tidal analysis.
;     Default is [tmin+3,tmax], where tmin and tmax are minimum and
;     maximum times.
;
;   VARIABLE (input, string 2-element vector)
;     The names of a pair of u & v type variables. Default is
;     ['ubar','vbar'].
;
; TO DO:
;   The square tiling stuff that was introduced for efficiency reasons is
;   very intricate. Perhaps it would be better to work in rows?
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-11:
;     Written.
;   Mark Hadfield, 2009-06:
;     - Added support for multiple, user-specifiable constituents and
;       user-specifiable variables.
;     - Data is now extracted using the MGHromsHistory class's HSliceGrid and
;       HsliceData methods, which allows the routine to accept DEPTHS, LEVELS or
;       SIGMAS keywords for 3D variables.
;   Mark Hadfield, 2009-08:
;     - Found and fixed a very sneaky little bug in tiling: eta_range calculation
;       had index 1 (one) instead of l (ell) in tile_count_y!
;   Mark Hadfield, 2009-09:
;     - Reduced target tile size from 128 to 32 to avoid out-of-memory errors when
;       processing 3D velocity data.
;-
pro mgh_roms_tide_vector_to_ncdf, history, file_out, $
     CONSTITUENTS=constituents, DEPTHS=depths, LEVELS=levels, SIGMAS=sigmas, $
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
         history_destroy = 1
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

   ;; Process VARIABLE keyword

   if n_elements(variable) eq 0 then $
        variable = ['ubar','vbar']

   ;; Process constituent list.

   if n_elements(constituents) eq 0 then constituents = ['M2']

   n_con = n_elements(constituents)

   ;; Check out the dimensionality of the variables

   dim_u = ohis->VarDims(variable[0])
   dim_v = ohis->VarDims(variable[1])

   if min([strlen(dim_u.horizontal),strlen(dim_u.time)]) eq 0 then $
        message, 'The U variable must have horizontal & time dimensions'
   if min([strlen(dim_v.horizontal),strlen(dim_v.time)]) eq 0 then $
        message, 'The V variable must have horizontal & time dimensions'
   if dim_u.time ne dim_v.time then $
        message, 'The U & V variables have different time dimensions'
   if dim_u.vertical ne dim_v.vertical then $
        message, 'The U & V variables have different vertical dimensions'
   if ~ array_equal(dim_u.horizontal, ['xi_u','eta_u']) then $
        message, 'The U variable has incorrect horizontal dimensions'
   if ~ array_equal(dim_v.horizontal, ['xi_v','eta_v']) then $
        message, 'The V variable has incorrect horizontal dimensions'

   ;; Process keywords relating to the vertical grid

   dim_vertical = dim_u.vertical
   has_vertical = strlen(dim_vertical) gt 0

   case has_vertical of

      0B: begin
         fmt = '(%"The %s keyword is not required or allowed when the ' + $
               'variable has no vertical dimension")'
         if n_elements(levels) gt 0 then message, string(FORMAT=fmt, 'LEVELS')
         if n_elements(depths) gt 0 then message, string(FORMAT=fmt, 'DEPTHS')
         if n_elements(sigmas) gt 0 then message, string(FORMAT=fmt, 'SIGMAS')
      end

      1B: begin

         theta_s = ohis->VarGet('theta_s')
         theta_b = ohis->VarGet('theta_b')
         hc = ohis->VarGet('hc')
         s = ohis->VarGet(dim_vertical)

         n_key = (n_elements(depths) gt 0) + (n_elements(levels) gt 0) + $
                 (n_elements(sigmas) gt 0)
         if n_key gt 1 then $
              message, 'The DEPTHS, LEVELS & SIGMAS keywords cannot be used together'
         if n_key eq 0 then $
              levels = lindgen(n_elements(s))
         mgh_undefine, n_key

         ;; By now we have exactly one of the depth-related keywords defined.
         case 1B of
            n_elements(levels) gt 0: n_vertical = n_elements(levels)
            n_elements(depths) gt 0: n_vertical = n_elements(depths)
            n_elements(sigmas) gt 0: n_vertical = n_elements(sigmas)
         endcase

      end

   endcase

   ;; Get time data from the history file and process the TIME_RANGE parameter.

   dim_time = dim_u.time

   case 1B of
      ohis->HasVar(dim_time): $
           time_var = dim_time
      ohis->HasVar('ocean_time'): $
           time_var = 'ocean_time'
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

   fmt = '(%"Vector tidal analysis (time-stamped) from ROMS dataset %s")'
   onc_out->AttAdd, /GLOBAL, 'long_name', $
        string(mgh_get_property(ohis, /NAME), FORMAT=fmt)

   onc_out->AttAdd, /GLOBAL, 'time_range', time_range

   ;; Create netCDF dimensions. The horizontal grid is the ROMS
   ;; interior rho grid

   hdim = [ohis->DimInfo('xi_rho', /DIMSIZE)-2, $
           ohis->DimInfo('eta_rho', /DIMSIZE)-2]

   onc_out->DimAdd, 'xi', hdim[0]
   onc_out->DimAdd, 'eta', hdim[1]

   if has_vertical then begin
      case 1B of
         n_elements(levels) gt 0: dim_vertical = 's'
         n_elements(depths) gt 0: dim_vertical = 'depth'
         n_elements(sigmas) gt 0: dim_vertical = 'sigma'
      endcase
      onc_out->DimAdd, dim_vertical, n_vertical
   endif

   onc_out->DimAdd, 'constituent', n_con
   onc_out->DimAdd, 'cindex', 1+2*n_con
   onc_out->DimAdd, 'name', 4

   ;; Create variables

   onc_out->VarAdd, 'cname', ['name','constituent'], /CHAR
   onc_out->AttAdd, 'cname', 'long_name', 'Constituent name'

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

   onc_out->VarAdd, 'angle', ['xi','eta'], /FLOAT
   onc_out->AttAdd, 'angle', 'long_name', 'angle between XI-axis and EAST'
   onc_out->AttAdd, 'angle', 'units', 'radian'

   if has_vertical then begin
      onc_out->VarAdd, dim_vertical, dim_vertical
   endif

   ;; Create result variables

   coeff_name = 'coeff_'+variable

   fill_real = 1.E10

   var_dim = ['cindex','xi','eta']
   if has_vertical then var_dim = [var_dim,dim_vertical]

   for v=0,1 do begin
      onc_out->VarAdd, coeff_name[v], var_dim, /FLOAT
      if ohis->HasAtt(variable[v], 'long_name') then begin
         att = ohis->AttGet(variable[v], 'long_name')+' coefficient'
         onc_out->AttAdd, coeff_name[v], 'long_name', temporary(att)
      endif
      if ohis->HasAtt(variable[v], 'units') then begin
         att = ohis->AttGet(variable[v], 'units')
         onc_out->AttAdd, coeff_name[v], 'units', temporary(att)
      endif
      onc_out->AttAdd, coeff_name[v], '_FillValue', fill_real
   endfor

   ;; Add constituent & grid data

   onc_out->VarPut, 'cname', constituents

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

   angle = ohis->VarGet('angle', OFFSET=[1,1], COUNT=hdim)
   onc_out->VarPut, 'angle', angle

   if has_vertical then begin
      case 1B of
         n_elements(levels) gt 0: $
              onc_out->VarPut, dim_vertical, s[levels]
         n_elements(depths) gt 0: $
              onc_out->VarPut, dim_vertical, depths
         n_elements(sigmas) gt 0: $
              onc_out->VarPut, dim_vertical, sigmas
      endcase
   endif

   ;; Set up coefficient array. If there is a trailing vertical
   ;; dimension of size 1, strip it off now to avoid surprises later.

   cdim = [1+2*n_con,hdim]
   if has_vertical then cdim = [cdim,n_vertical]

   coeff = reform(make_array(cdim, VALUE=complex(fill_real, fill_real)))

   ;; Create an MGHtideParameters object to use in the analysis

   oprm = obj_new('MGHtideParameters', /DEFAULT_PARAMETERS)

   ;; Retrieve data and calculate coefficients. A tiling scheme is
   ;; used for retrieving data to speed up access while avoiding large
   ;; memory requirements. Tile offsets are expressed relative to the
   ;; SW interior rho point.

   n_tile = round(hdim/32.D0) > 2

   tile_size = floor(hdim/double(n_tile))

   tile_offset_x = tile_size[0]*lindgen(n_tile[0])
   tile_offset_y = tile_size[1]*lindgen(n_tile[1])

   mgh_undefine, tile_size

   tile_count_x = [tile_offset_x[1:*],hdim[0]] - tile_offset_x
   tile_count_y = [tile_offset_y[1:*],hdim[1]] - tile_offset_y

   for l=0,n_tile[1]-1 do begin
      for k=0,n_tile[0]-1 do begin

         fmt = '(%"Processing tile %d, %d, offset %d, %d")'
         print, FORMAT=fmt, k, l, tile_offset_x[k], tile_offset_y[l]

         xi_range = 1 + tile_offset_x[k] + [0,tile_count_x[k]-1]
         eta_range = 1 + tile_offset_y[l] + [0,tile_count_y[l]-1]

         grid_u = ohis->HsliceGrid(variable[0], $
                                   XI_RANGE=xi_range+[-1,0], ETA_RANGE=eta_range)
         grid_v = ohis->HsliceGrid(variable[1], $
                                   XI_RANGE=xi_range, ETA_RANGE=eta_range+[-1,0])

         for r=0,rrn-1 do begin

            ;; Get data. Note that if the number of DEPTHS, LEVELS, or
            ;; SIGMAS is 1, then the trailing unit dimension will be
            ;; stripped off.
            uu = ohis->HsliceData(variable[0], GRID=grid_u, MASK_VALUE=0, $
                                  RECORD=rr0+r, DEPTHS=depths, LEVELS=levels, $
                                  SIGMAS=sigmas)
            vv = ohis->HsliceData(variable[1], GRID=grid_v, MASK_VALUE=0, $
                                  RECORD=rr0+r, DEPTHS=depths, LEVELS=levels, $
                                  SIGMAS=sigmas)

            uu = mgh_stagger(uu, DELTA=[-1,0,0])
            vv = mgh_stagger(vv, DELTA=[0,-1,0])

            if r eq 0 then begin
               u = fltarr([size(uu, /DIMENSIONS),rrn])
               v = fltarr([size(vv, /DIMENSIONS),rrn])
               n_vdim = size(u, /N_DIMENSIONS)
            endif

            case n_vdim of
               3: begin
                  u[*,*,r] = temporary(uu)
                  v[*,*,r] = temporary(vv)
               end
               4: begin
                  u[*,*,*,r] = temporary(uu)
                  v[*,*,*,r] = temporary(vv)
               end
            endcase

         endfor

         for j=0,tile_count_y[l]-1 do begin
            for i=0,tile_count_x[k]-1 do begin

               ii = tile_offset_x[k] + i
               jj = tile_offset_y[l] + j

               if (~ has_mask) || mask[ii,jj] then begin
                  case n_vdim of
                     3: begin
                        mgh_tide_vector_analysis, $
                             oprm, time+time_ref, reform(u[i,j,*]), reform(v[i,j,*]), cc, $
                             CONSTITUENTS=constituents
                        coeff[*,ii,jj] = cc
                     end
                     4: begin
                        for m=0,n_vertical-1 do begin
                           mgh_tide_vector_analysis, $
                                oprm, time+time_ref, reform(u[i,j,m,*]), reform(v[i,j,m,*]), cc, $
                                CONSTITUENTS=constituents
                           coeff[*,ii,jj,m] = cc
                        endfor
                     end
                  endcase
               endif

            endfor
         endfor

      endfor
   endfor

   ;; Write data & clean up

   onc_out->VarPut, coeff_name[0], real_part(coeff)
   onc_out->VarPut, coeff_name[1], imaginary(coeff)

   obj_destroy, oprm
   obj_destroy, onc_out
   if history_destroy then obj_destroy, ohis

end
