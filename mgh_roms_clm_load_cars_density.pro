;+
; NAME:
;   MGH_ROMS_CLM_LOAD_CARS_DENSITY
;
; PURPOSE:
;   Load CARS hydrography data into a ROMS climate file.
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_LOAD_CARS_DENSITY, file_clm
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file into which CARS data are to
;     be written.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-04:
;     Written, based on MGH_ROMS_CLM_LOAD_WOA_DENSITY.
;   Mark Hadfield, 2004-01:
;     Generalised to allow annual (monthly) fields. Names of time-related
;     dimensions & variables have been changed.
;   Mark Hadfield, 2004-09:
;     - Fixed bug: "flattened" temperature & salinity fields were actually
;       being set to zero.
;     - Vertical interpolation now quadratic.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;   Mark Hadfield, 2010-01:
;     - Added DATASET keyword.
;     - Handling of annual & semi-annual variation modified to fit
;       the scheme used in CARS 2006 and CARS 2009.
;   Mark Hadfield, 2010-10:
;     Reduced memory required to construct the time-varying CARS fields
;     by removing a call to REBIN. The memory requirement could be
;     further reduced by bringing the writing of data inside the loop
;     over time.
;   Mark Hadfield, 2014-01:
;     Updated indentation.
;-
function mgh_roms_clm_load_cars_density_flatten, data

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  dim = size(data, /DIMENSIONS)

  result = fltarr(dim)

  for k=0,dim[2]-1 do $
    result[*,*,k] = mgh_avg(data[*,*,k])

  return, result

end

function mgh_roms_clm_load_cars_density_var3d, data, cars, grid

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Interpolate CARS 3D data to the model grid.

  n_depth = n_elements(cars.depth)

  ;; ...first horizontally

  hdata = fltarr([grid.dim[0:1], n_depth])

  x = mgh_locate(cars.lon, XOUT=grid.lon[*])
  y = mgh_locate(cars.lat, XOUT=grid.lat[*])

  for k=0,n_depth-1 do begin
    hdata[*,*,k] = reform(interpolate(data[*,*,k], x, y, CUBIC=-0.5), grid.dim[0], grid.dim[1])
  endfor

  ;; ...then vertically to the s levels. Note that quadratic interpolation
  ;; may give (wildly) incorrect results if depths are not supplied
  ;; as floating point.

  gdata = fltarr(grid.dim)

  for i=0,grid.dim[0]-1 do begin
    for j=0,grid.dim[1]-1 do begin
      gdata[i,j,*] = interpol(reform(hdata[i,j,*]), float(cars.depth), -reform(grid.z[i,j,*]), /QUADRATIC)
    endfor
  endfor

  return, gdata

end

pro mgh_roms_clm_load_cars_density, file_clm, $
  ANNUAL=annual, DATASET=dataset, FLATTEN=flatten, MARGIN=margin, $
  TAPER_ANNUAL=taper_annual, TIME_NAME=time_name

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

  ;; Process keywords

  if n_elements(margin) eq 0 then margin = {lon: [-10,10], lat: [-10,10]}

  if n_elements(time_name) eq 0 then time_name = 'time'

  ;; Open restart file as MGHromsHistory object to get dimensions &
  ;; grid data. A bit messy this!

  fmt = '(%"Opening climatology file %s for read-only access")'
  message, /INFORM, string(FORMAT=fmt, file_clm)

  oclm = obj_new('MGHromsHistory', file_clm)

  grid = {dim: oclm->DimRho(), lon: oclm->VarGet('lon_rho'), $
          lat: oclm->VarGet('lat_rho'), z: oclm->GetZGrid()}

  obj_destroy, oclm

  ;; Open the file

  fmt = '(%"Opening climatology file %s for modification")'
  message, /INFORM, string(FORMAT=fmt, file_clm)

  oclm = obj_new('MGHncFile', file_clm, /MODIFY)

  ;; Check that the required time dimension & variable are available

  if ~ oclm->HasDim(time_name) then message, 'Time dimension missing'
  if ~ oclm->HasVar(time_name) then message, 'Time variable missing'

  ;; Add variables to hold CARS data

  if ~ oclm->HasVar('temp') then begin
    message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'temp')
    oclm->VarAdd, 'temp', ['xi_rho','eta_rho','s_rho',time_name]
    oclm->AttAdd, 'temp', 'time', time_name
  endif

  if ~ oclm->HasVar('salt') then begin
    message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'salt')
    oclm->VarAdd, 'salt', ['xi_rho','eta_rho','s_rho',time_name]
    oclm->AttAdd, 'salt', 'time', time_name
  endif

  ;; Specify boundaries over which to get hydrographic data

  bound = {lon: mgh_minmax(grid.lon)+margin.lon, lat: mgh_minmax(grid.lat)+margin.lat}

  ;; Get CARS temperature & salinity data

  message, /INFORM, 'Getting CARS temperature & salinity data'

  cars_temp = mgh_cars_volume('temperature', DATASET=dataset, BOUNDARIES=bound)
  cars_salt = mgh_cars_volume('salinity', DATASET=dataset, BOUNDARIES=bound)

  cars = {n_depth: cars_temp.n_depth, n_ann: cars_temp.n_ann, n_semiann: cars_temp.n_semiann, $
          n_lon: n_elements(cars_temp.lon), n_lat: n_elements(cars_temp.lat), $
          lon: cars_temp.lon, lat: cars_temp.lat, depth: cars_temp.depth}

  if keyword_set(taper_annual) then begin
    fac = 0.5*(1.-tanh((cars.depth[0:cars.n_ann-1]-200.)/25.))
    fac = rebin(reform(fac, [1,1,cars.n_ann]), [cars.n_lon,cars.n_lat,cars.n_ann])
    *cars_temp.an_cos *= fac
    *cars_temp.an_sin *= fac
    *cars_salt.an_cos *= fac
    *cars_salt.an_sin *= fac
    mgh_undefine, fac
    fac = 0.5*(1.-tanh((cars.depth[0:cars.n_semiann-1]-200.)/25.))
    fac = rebin(reform(fac, [1,1,cars.n_semiann]), [cars.n_lon,cars.n_lat,cars.n_semiann])
    *cars_temp.sa_cos *= fac
    *cars_temp.sa_sin *= fac
    *cars_salt.sa_cos *= fac
    *cars_salt.sa_sin *= fac
    mgh_undefine, fac
  endif

  case keyword_set(annual) of

    0B: begin
      n_time = 1
      ctemp = *cars_temp.mean
      csalt = *cars_salt.mean
    end

    1B: begin
      time = oclm->VarGet(time_name)
      n_time = n_elements(time)
      cycle = oclm->AttGet(time_name, 'cycle_length')
      ctemp = fltarr([cars.n_lon,cars.n_lat,cars.n_depth,n_time])
      csalt = fltarr([cars.n_lon,cars.n_lat,cars.n_depth,n_time])
      for t=0,n_time-1 do begin
        ctemp[*,*,*,t] = *cars_temp.mean
        csalt[*,*,*,t] = *cars_salt.mean
        ang = 2*!pi*time[t]/float(cycle)
        ctemp[*,*,0:cars.n_ann-1,t] += $
          *cars_temp.an_cos*cos(ang) + *cars_temp.an_sin*sin(ang)
        ctemp[*,*,0:cars.n_semiann-1,t] += $
          *cars_temp.sa_cos*cos(2*ang) + *cars_temp.sa_sin*sin(2*ang)
        csalt[*,*,0:cars.n_ann-1,t] += $
          *cars_salt.an_cos*cos(ang) + *cars_salt.an_sin*sin(ang)
        csalt[*,*,0:cars.n_semiann-1,t] += $
          *cars_salt.sa_cos*cos(2*ang) + *cars_salt.sa_sin*sin(2*ang)
      endfor
    end

  endcase

  mgh_undefine, cars_temp, cars_salt

  ;; Fill CARS temperature & salinity data

  message, /INFORM, 'Filling CARS temperature & salinity data'

  for t=0,n_time-1 do begin
    for k=0,cars.n_depth-1 do begin
      ctemp[*,*,k,t] = mgh_fill2d(ctemp[*,*,k,t])
      csalt[*,*,k,t] = mgh_fill2d(csalt[*,*,k,t])
    endfor
  endfor

  ;; Optionally flatten (homogenise) fields

  if keyword_set(flatten) then begin

    message, /INFORM, 'Flattening WOA temperature & salinity data'

    for t=0,n_time-1 do begin
      ctemp[0,0,0,t] = mgh_roms_clm_load_cars_density_flatten(ctemp[*,*,*,t])
      csalt[0,0,0,t] = mgh_roms_clm_load_cars_density_flatten(csalt[*,*,*,t])
    endfor

  endif

  ;; Calculate potential temperature

  message, /INFORM, 'Calculating CARS potential temperature'

  cptmp = ctemp
  for t=0,n_time-1 do begin
    for i=0,cars.n_lon-1 do begin
      for j=0,cars.n_lat-1 do begin
        cptmp[i,j,0,t] = $
          mgh_sw_ptmp(ctemp[i,j,*,t], csalt[i,j,*,t], $
          mgh_sw_pres(cars.depth, 45))
      endfor
    endfor
  endfor
  mgh_undefine, ctemp

  ;; Interpolate temperature & salinity to the model grid and write each variable to
  ;; the file

  message, /INFORM, 'Interpolating CARS pot. temp. to the model grid'

  gptmp = fltarr([grid.dim,n_time])
  for t=0,n_time-1 do $
    gptmp[0,0,0,t] = mgh_roms_clm_load_cars_density_var3d(cptmp[*,*,*,t], cars, grid)

  message, /INFORM, 'Writing pot. temp. to the output file'
  oclm->VarPut, 'temp', temporary(gptmp)

  message, /INFORM, 'Interpolating CARS salinity to the model grid'

  gsalt = fltarr([grid.dim,n_time])
  for t=0,n_time-1 do $
    gsalt[0,0,0,t] = mgh_roms_clm_load_cars_density_var3d(csalt[*,*,*,t], cars, grid)

  message, /INFORM, 'Writing salinity to the output file'
  oclm->VarPut, 'salt', temporary(gsalt)

  ;; Close the output file

  obj_destroy, oclm

end

