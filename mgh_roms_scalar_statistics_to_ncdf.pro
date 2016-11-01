;+
; NAME:
;   MGH_ROMS_SCALAR_STATISTICS_TO_NCDF
;
; PURPOSE:
;   Get scalar data from a ROMS file, calculate one or more statistical
;   quantities and save the results to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   mgh_roms_scalar_statistics_to_ncdf, history, file_out
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
;   PARAMETERS (input, string vector)
;     A list of statistical parameters. Default is ['mean'].
;
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
;     it cannot be used together with DEPTH or SIGMA.
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
;   Mark Hadfield, 2010-07:
;     Written, based on the file-handling and data-handling framework
;     of MGH_ROMS_TIDE_SCALAR_TO_NCDF and the statistics code of the
;     Mgh_Roms_Plot_Hstats class.
;   Mark Hadfield, 2011-12:
;     - Got rid of the history_destroy functionality.
;     - Variable names now converted with MGH_STR_VANILLA in forming
;       netcDF variable names.
;   Mark Hadfield, 2013-04:
;     - Replaced DEPTHS, LEVELS & SIGMAS keywords with DEPTH,
;       LEVEL & SIGMA.
;-
pro mgh_roms_scalar_statistics_to_ncdf, history, file_out, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     PARAMETERS=parameters, TIME_RANGE=time_range, VARIABLE=variable

    compile_opt DEFINT32
    compile_opt STRICTARR
    compile_opt STRICTARRSUBS
    compile_opt LOGICAL_PREDICATE

    ;; Process history argument.

    case size(history, /TNAME) of
       'STRING': begin
          ohis = obj_new('MGHromsHistory', history)
       end
       'OBJREF': begin
          ohis = history
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

    ;; Process parameter list.

    if n_elements(parameters) eq 0 then parameters = ['mean']

    n_par = n_elements(parameters)

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
       s = ohis->VarGet(dim_vertical)
       n_key = (n_elements(depth) gt 0) + (n_elements(level) gt 0) + $
          (n_elements(sigma) gt 0)
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
          'variable has no vertical dimension")'
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

    fmt = '(%"Scalar statistics from ROMS dataset %s")'
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
          n_elements(level) gt 0: dim_vertical = 's'
          n_elements(depth) gt 0: dim_vertical = 'depth'
          n_elements(sigma) gt 0: dim_vertical = 'sigma'
       endcase
       onc_out->DimAdd, dim_vertical, n_vertical
    endif

    onc_out->DimAdd, 'parameter', n_par
    onc_out->DimAdd, 'name', max(strlen(parameters))

    ;; Create variables

    onc_out->VarAdd, 'parameter', ['name','parameter'], /CHAR
    onc_out->AttAdd, 'parameter', 'long_name', 'Parameter name'

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

    ;; Create netCDF variable for result

    result_name = mgh_str_vanilla(variable)+'_statistics'

    fill_real = 1.E10

    result_dim = ['xi','eta']
    if has_vertical then result_dim = [result_dim,dim_vertical]
    result_dim = [result_dim,'parameter']

    onc_out->VarAdd, result_name, result_dim, /FLOAT
    onc_out->AttAdd, result_name, '_FillValue', fill_real

    ;; Add parameter & grid data

    onc_out->VarPut, 'parameter', parameters

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
    dim = [dim,n_par]

    result = make_array(dim, VALUE=fill_real)

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

                   for p=0,n_par-1 do begin

                      par = parameters[p]

                      case 1B of
                         strmatch(par, 'mean', /FOLD_CASE): begin
                            if has_vertical then begin
                               result[ii,jj,*,p] = mgh_avg(reform(data[i,j,*,*]), 1)
                            endif else begin
                               result[ii,jj,p] = mgh_avg(reform(data[i,j,*]))
                            endelse
                         end
                         strmatch(par, 'median', /FOLD_CASE): begin
                            if has_vertical then begin
                               result[ii,jj,*,p] = median(reform(data[i,j,*,*]), DIM=1)
                            endif else begin
                               result[ii,jj,p] = median(reform(data[i,j,*]))
                            endelse
                         end
                         strmatch(par, 'standard deviation', /FOLD_CASE): begin
                            if has_vertical then begin
                               result[ii,jj,*,p] = mgh_stdev(reform(data[i,j,*,*]), 1)
                            endif else begin
                               result[ii,jj,p] = mgh_stdev(reform(data[i,j,*]))
                            endelse
                         end
                         strmatch(par, 'max', /FOLD_CASE): begin
                            if has_vertical then begin
                               result[ii,jj,*,p] = max(reform(data[i,j,*,*]), DIM=1)
                            endif else begin
                               result[ii,jj,p] = max(reform(data[i,j,*]))
                            endelse
                         end
                         strmatch(par, 'min', /FOLD_CASE): begin
                            if has_vertical then begin
                               result[ii,jj,p] = min(reform(data[i,j,*]))
                            endif else begin
                               result[ii,jj,*,p] = min(reform(data[i,j,*,*]), DIM=1)
                            endelse
                         end
                         strmatch(par, 'percentile*', /FOLD_CASE): begin
                            ;; It's a bit inefficient to do this for
                            ;; every grid point!
                            pp = strsplit(par, /EXTRACT)
                            threshold = n_elements(pp) gt 1 ? float(pp[1]) : 50
                            if has_vertical then begin
                               for m=0,n_vertical-1 do begin
                                  result[ii,jj,m,p] = $
                                     mgh_percentile(reform(data[i,j,m,*]), THRESHOLD=threshold, /NAN)
                               endfor
                            endif else begin
                               result[ii,jj,p] = $
                                  mgh_percentile(reform(data[i,j,*]), THRESHOLD=threshold, /NAN)
                            endelse
                         end
                         strmatch(par, 'fraction > ?*', /FOLD_CASE): begin
                            ;; It's a bit inefficient to do this for
                            ;; every grid point!
                            pp = strsplit(par, /EXTRACT)
                            threshold = float(pp[2])
                            if has_vertical then begin
                               for m=0,n_vertical-1 do begin
                                  dd = reform(data[i,j,m,*])
                                  l_good = where(finite(dd), n_good)
                                  if n_good gt 0 then begin
                                     dd = dd[l_good]
                                     result[ii,jj,m,p] = total(dd gt threshold)/n_elements(dd)
                                  endif
                               endfor
                            endif else begin
                               dd = reform(data[i,j,*])
                               l_good = where(finite(dd), n_good)
                               if n_good gt 0 then begin
                                  dd = dd[l_good]
                                  result[ii,jj,p] = total(dd gt threshold)/n_elements(dd)
                               endif
                            endelse
                         end
                      endcase

                   endfor

                endif

             endfor
          endfor

       endfor
    endfor

    ;; Write data & clean up

    onc_out->VarPut, result_name, result

    obj_destroy, onc_out

end
