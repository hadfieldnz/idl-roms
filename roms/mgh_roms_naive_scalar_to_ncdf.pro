;+
; NAME:
;   MGH_ROMS_NAIVE_SCALAR_TO_NCDF
;
; PURPOSE:
;   Get scalar data from a ROMS file, carry out a naive
;   (non-time-stamped) tidal analysis and save the results to a netCDF
;   file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   mgh_roms_naive_scalar_to_ncdf, history, file_out
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
;   VARIABLE (input, string scalar)
;     The name of a scalar variable to be analysed. Default is 'zeta'.
;
; TO DO:
;   The square tiling stuff that was introduced for efficiency reasons is
;   very intricate. Perhaps it would be better to work in rows?
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-11:
;     Written.
;   Mark Hadfield, 2009-02:
;     - Added support for multiple, user-specifiable constituents and
;       user-specifiable variables.
;     - Data is now extracted using the MGHromsHistory class's HSliceGrid and
;       HsliceData methods, which allows the routine to accept DEPTHS, LEVELS or
;       SIGMAS keywords for 3D variables.
;   Mark Hadfield, 2009-08:
;     - Found and fixed a very sneaky little bug in tiling: eta_range calculation
;       had index 1 (one) instead of l (ell) in tile_count_y!
;   Mark Hadfield, 2010-07:
;     Found and fixed a couple of bugs:
;     - A vector analysis was carried out rather than a scalar analysis for
;       vertically varying data. (I don't think this part of the code was ever
;       exercised.)
;     - The coefficient array was complex-valued. (This had no effect as the
;       imaginary part was never calculated and ignored when written to the
;       netCDF file.)
;     Added code to write out the rms magnitude of the fit's residual.
;-
pro mgh_roms_naive_scalar_to_ncdf, history, file_out, $
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
        variable = 'zeta'

   ;; Process constituent list. The analysis routine fits the data
   ;; with one or more sinusoids and requires a frequency for each
   ;; one. Here we assume the sinusoids correspond to tidal
   ;; constituents and we get frequencies from tidal constant
   ;; tables. Phase information associated with the constituents is
   ;; ignored.

   if n_elements(constituents) eq 0 then constituents = ['M2']

   n_con = n_elements(constituents)

   const = mgh_t_getconsts()

   l_con = lonarr(n_con)
   for c=0,n_con-1 do begin
      l_con[c] = where(strmatch(const.name, constituents[c], /FOLD_CASE))
      if l_con[c] lt 0 then $
           message, 'Could not find constituent: '+constituents[c]
   endfor

   ;; Convert to d^-1
   frequency = 24.0D*const.freq[l_con]

   mgh_undefine, const

   ;; Check out the dimensionality of the variables

   dim = ohis->VarDims(variable)

   if min([strlen(dim.horizontal),strlen(dim.time)]) eq 0 then $
        message, 'The variable must have horizontal & time dimensions'
   if ~ array_equal(dim.horizontal, ['xi_rho','eta_rho']) then $
        message, 'I can handle only variables on the horizontal RHO grid right now.'

   ;; Process keywords relating to the vertical grid

   dim_vertical = dim.vertical
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

   dim_time = dim.time

   time_var = ohis->TimeVarName(dim_time)
   if isa(time_var, /NULL) then message, 'Time variable not found'

   time_units = {scale: 1D/(24D*3600D), offset: 0D}
   if ohis->MGHncSequence::HasAtt(time_var, 'units') then $
      time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))

   time = ohis->VarGet(time_var, AUTOSCALE=0)*time_units.scale

   if n_elements(time_range) eq 0 then time_range = mgh_minmax(time) + [3,0]

   rrange = mgh_subset(time, time_range+[-1,1]*1.0D-6, EMPTY=empty)

   if empty then message, "Record range is empty"

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

   fmt = '(%"Scalar tidal analysis (naive) from ROMS dataset %s")'
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

   if has_vertical then begin
      onc_out->VarAdd, dim_vertical, dim_vertical
   endif

   ;; Create coefficient variable

   coeff_name = 'coeff_'+variable

   fill_real = 1.E10

   var_dim = ['cindex','xi','eta']
   if has_vertical then var_dim = [var_dim,dim_vertical]

   onc_out->VarAdd, coeff_name, var_dim, /FLOAT
   if ohis->HasAtt(variable, 'long_name') then begin
      att = ohis->AttGet(variable, 'long_name')+' coefficient'
      onc_out->AttAdd, coeff_name, 'long_name', temporary(att)
   endif
   if ohis->HasAtt(variable, 'units') then begin
      att = ohis->AttGet(variable, 'units')
      onc_out->AttAdd, coeff_name, 'units', temporary(att)
   endif
   onc_out->AttAdd, coeff_name, '_FillValue', fill_real

   ;; Create residual magnitude variable

   res_name = 'res_'+variable

   fill_real = 1.E10

   var_dim = ['xi','eta']
   if has_vertical then var_dim = [var_dim,dim_vertical]

   onc_out->VarAdd, res_name, var_dim, /FLOAT
   if ohis->HasAtt(variable, 'long_name') then begin
      att = ohis->AttGet(variable, 'long_name')+' residual rms'
      onc_out->AttAdd, res_name, 'long_name', temporary(att)
   endif
   if ohis->HasAtt(variable, 'units') then begin
      att = ohis->AttGet(variable, 'units')
      onc_out->AttAdd, res_name, 'units', temporary(att)
   endif
   onc_out->AttAdd, res_name, '_FillValue', fill_real

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

   ;; Set up coefficient array.

   cdim = [1+2*n_con,hdim]
   if has_vertical then cdim = [cdim,n_vertical]

   coeff = reform(make_array(cdim, VALUE=fill_real))

   ;; Set up residual magnitude array.

   rdim = [hdim]
   if has_vertical then rdim = [rdim,n_vertical]

   res = reform(make_array(rdim, VALUE=fill_real))

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

         fmt = '(%"Processing tile %d, %d, offset %d, %d")'
         print, FORMAT=fmt, k, l, tile_offset_x[k], tile_offset_y[l]

         xi_range = 1 + tile_offset_x[k] + [0,tile_count_x[k]-1]
         eta_range = 1 + tile_offset_y[l] + [0,tile_count_y[l]-1]

         grid = ohis->HsliceGrid(variable, $
                                 XI_RANGE=xi_range, ETA_RANGE=eta_range)

         for r=0,rrn-1 do begin

            ;; Get data. Note that if the number of DEPTHS, LEVELS, or
            ;; SIGMAS is 1, then the trailing unit dimension will be
            ;; stripped off.
            ss = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=0, $
                                  RECORD=rr0+r, DEPTHS=depths, LEVELS=levels, $
                                  SIGMAS=sigmas)

            if r eq 0 then begin
               s = fltarr([size(ss, /DIMENSIONS),rrn])
               n_vdim = size(s, /N_DIMENSIONS)
            endif

            case n_vdim of
               3: begin
                  s[*,*,r] = temporary(ss)
               end
               4: begin
                  s[*,*,*,r] = temporary(ss)
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
                        mgh_naive_scalar_analysis, $
                             time, reform(s[i,j,*]), cc, $
                             FREQUENCY=frequency, RESIDUAL=residual
                        coeff[*,ii,jj] = cc
                        res[ii,jj] = sqrt(mean(residual^2))
                     end
                     4: begin
                        for m=0,n_vertical-1 do begin
                           mgh_naive_scalar_analysis, $
                                time, reform(s[i,j,m,*]), cc, $
                                FREQUENCY=frequency, RESIDUAL=residual
                           coeff[*,ii,jj,m] = cc
                           res[ii,jj,mm] = sqrt(mean(residual^2))
                        endfor
                     end
                  endcase
               endif

            endfor
         endfor

      endfor
   endfor

   ;; Write data & clean up

   onc_out->VarPut, coeff_name, coeff

   onc_out->VarPut, res_name, res

   obj_destroy, onc_out
   if history_destroy then obj_destroy, ohis

end
