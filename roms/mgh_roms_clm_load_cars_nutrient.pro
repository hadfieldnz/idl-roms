; svn $Id$
;+
; NAME:
;   MGH_ROMS_CLM_LOAD_CARS_NUTRIENT
;
; PURPOSE:
;   Load CARS nutrient (nitrate only) data into a ROMS climate file.
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_LOAD_CARS_NUTRIENT, file_clm
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file into which CARS data are to
;     be written.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2004-09:
;     Written, based on MGH_ROMS_CLM_LOAD_CARS_DENSITY.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;-

function mgh_roms_clm_load_cars_nutrient_flatten, data

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

function mgh_roms_clm_load_cars_nutrient_var3d, data, cars, grid

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
      hdata[*,*,k] = reform(interpolate(data[*,*,k], x, y, CUBIC=-0.5), $
                            grid.dim[0], grid.dim[1])
   endfor

   ;; ...then vertically to the s levels. Note that quadratic interpolation
   ;; may give (wildly) incorrect results if depths are not supplied
   ;; as floating point.

   gdata = fltarr(grid.dim)

   for i=0,grid.dim[0]-1 do begin
      for j=0,grid.dim[1]-1 do begin
         gdata[i,j,*] = interpol(reform(hdata[i,j,*]), float(cars.depth), $
                                 -reform(grid.z[i,j,*]), /QUADRATIC)
      endfor
   endfor

   return, gdata

end

pro mgh_roms_clm_load_cars_nutrient, file_clm, $
     ANNUAL=annual, FLATTEN=flatten, MARGIN=margin, TIME_NAME=time_name

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

   grid = {dim: oclm->DimRho(), $
           lon: oclm->VarGet('lon_rho'), $
           lat: oclm->VarGet('lat_rho'), $
           z: oclm->GetZGrid() $
          }

   obj_destroy, oclm

   ;; Open the file

   fmt = '(%"Opening climatology file %s for modification")'
   message, /INFORM, string(FORMAT=fmt, file_clm)

   oclm = obj_new('MGHncFile', file_clm, /MODIFY)

   ;; Check that the required time dimension & variable are available

   if ~ oclm->HasDim(time_name) then message, 'Time dimension missing'
   if ~ oclm->HasVar(time_name) then message, 'Time variable missing'

   ;; Add variables to hold CARS data

   if ~ oclm->HasVar('NO3') then begin
      message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'NO3')
      oclm->VarAdd, 'NO3', ['xi_rho','eta_rho','s_rho',time_name]
      oclm->AttAdd, 'NO3', 'time', time_name
   endif

   ;; Specify boundaries over which to get nutrient data

   bound = {lon: mgh_minmax(grid.lon)+margin.lon, lat: mgh_minmax(grid.lat)+margin.lat}

   ;; Get CARS data

   message, /INFORM, 'Getting CARS nitrate data'

   cars_nitr = mgh_cars_volume('nitrate', DATASET='cars_2000', BOUNDARIES=bound)

   cars = {n_depth: cars_nitr.n_depth, n_dfit: cars_nitr.n_dfit, $
           n_lon: n_elements(cars_nitr.lon), n_lat: n_elements(cars_nitr.lat), $
           lon: cars_nitr.lon, lat: cars_nitr.lat, depth: cars_nitr.depth}

   case keyword_set(annual) of

      0B: begin
         n_time = 1
         cnitr = cars_nitr.mean
      end

      1B: begin
         time = oclm->VarGet(time_name)
         n_time = n_elements(time)
         cycle = oclm->AttGet(time_name, 'cycle_length')
         cnitr = rebin(*cars_nitr.mean, [cars.n_lon, cars.n_lat, cars.n_depth, n_time])
         for t=0,n_time-1 do begin
            ang = 2*!pi*time[t]/float(cycle)
            cnitr[*,*,0:cars.n_dfit-1,t] += $
                 *cars_nitr.an_cos*cos(ang) + *cars_nitr.an_sin*sin(ang)
         endfor
      end

   endcase

   mgh_undefine, cars_nitr

   ;; Fill CARS nutrient data. This can be quite
   ;; time-consuming.

   message, /INFORM, 'Filling CARS nutrient data'

   for t=0,n_time-1 do begin
      cnitr[0,0,0,t] = mgh_hydro2_fill(cnitr[*,*,*,t], cars)
   endfor

   ;; Optionally flatten (homogenise) fields

   if keyword_set(flatten) then begin

      message, /INFORM, 'Flattening CARS nutrient data'

      for t=0,n_time-1 do begin
         cnitr[0,0,0,t] = mgh_roms_clm_load_cars_nutrient_flatten(cnitr[*,*,*,t])
      endfor

   endif

   ;; Interpolate nutrient to the model grid

   message, /INFORM, 'Interpolating CARS nitrate to the model grid'

   gnitr = fltarr([grid.dim,n_time])
   for t=0,n_time-1 do begin
      gnitr[0,0,0,t] = $
           mgh_roms_clm_load_cars_nutrient_var3d(cnitr[*,*,*,t], cars, grid)
   endfor

   ;; Clip unphysical values

   gnitr = gnitr > 0.05

   ;; Write nutrient data to the output file

   message, /INFORM, 'Writing nitrate to the output file'
   oclm->VarPut, 'NO3', gnitr

   ;; Close the output file

   obj_destroy, oclm

end

