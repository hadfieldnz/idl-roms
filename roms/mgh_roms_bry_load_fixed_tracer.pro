;+
; NAME:
;   MGH_ROMS_BRY_LOAD_FIXED_TRACER
;
; PURPOSE:
;   Load zero data for passive or sediment tracers to a ROMS boundary file.
;
; CALLING SEQUENCE:
;   mgh_roms_bry_load_fixed_tracer, file_bry
;
; POSITIONAL PARAMETERS:
;   file_bry (input, scalar string)
;     The name of a ROMS boundary file into which data are to
;     be written.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2011-07:
;     Written as MGH_ROMS_BRY_LOAD_VANILLA_SEDIMENT.
;   Mark Hadfield, 2011-08:
;     Added support for dye (passive tracer) variables.
;   Mark Hadfield, 2013-08:
;     Renamed MGH_ROMS_BRY_LOAD_ZERO_TRACER.
;   Mark Hadfield, 2014-06:
;     - Now allows non-zero values to be specified for each tracer.
;     - Renamed MGH_ROMS_BRY_LOAD_FIXED_TRACER.
;   Mark Hadfield, 2015-03:
;     - Fixed bug: specified values (dye_value, mud_value, sand_value)
;       were being loaded on north and south boundaries, but zeroes
;       were being loaded on the east & west boundaries.
;-;   Mark Hadfield, 2019-05:
;     - The routine no longer assumes the time dimension has only one record,
;       but writes data to all records. This is a work-around for apparent
;       problems with ROMS' handling of calendar=none variables.
;-

pro mgh_roms_bry_load_fixed_tracer, file_bry, $
      N_DYE=n_dye, N_MUD=n_mud, N_SAND=n_sand, $
      DYE_VALUE=dye_value, MUD_VALUE=mud_value, SAND_VALUE=sand_value, $
      TIME_NAME=time_name

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

   ;; Process keywords

   if n_elements(time_name) eq 0 then time_name = 'time_fixed'

   if n_elements(n_dye) eq 0 then n_dye = 0

   if n_elements(n_mud) eq 0 then n_mud = 0

   if n_elements(n_sand) eq 0 then n_sand = 0

   if n_elements(dye_value) eq 0 then begin
    if n_dye gt 0 then dye_value = replicate(0.0, n_dye)
   endif

   if n_elements(mud_value) eq 0 then begin
     if n_mud gt 0 then mud_value = replicate(0.0, n_mud)
   endif

   if n_elements(sand_value) eq 0 then begin
     if n_sand gt 0 then sand_value = replicate(0.0, n_sand)
   endif

   ;; Open the file

   fmt = '(%"Opening boundary file %s for modification")'
   message, /INFORM, string(FORMAT=fmt, file_bry)

   obry = obj_new('MGHncFile', file_bry, /MODIFY)

   ;; Check that the required time dimension & variable are available.
   ;; Establish the number of time values

   if ~ obry->HasDim(time_name) then message, 'Time dimension missing'
   if ~ obry->HasVar(time_name) then message, 'Time variable missing'

   n_time = obry->DimInfo(time_name, /DIMSIZE)

   ;; Get grid dimensions

   grd_dim = {xi_rho: obry->DimInfo('xi_rho', /DiMSIZE), $
              eta_rho: obry->DimInfo('eta_rho', /DiMSIZE), $
              s_rho: obry->DimInfo('s_rho', /DiMSIZE)}

   ;; Add variables to hold data

   bname = ['south','north']

   for b=0,1 do begin
      for i=0,n_dye-1 do begin
         name = string(FORMAT='(%"dye_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['xi_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
      for i=0,n_mud-1 do begin
         name = string(FORMAT='(%"mud_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['xi_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
      for i=0,n_sand-1 do begin
         name = string(FORMAT='(%"sand_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['xi_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
   endfor

   bname = ['west','east']

   for b=0,1 do begin
      for i=0,n_dye-1 do begin
         name = string(FORMAT='(%"dye_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['eta_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
      for i=0,n_mud-1 do begin
         name = string(FORMAT='(%"mud_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['eta_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
      for i=0,n_sand-1 do begin
         name = string(FORMAT='(%"sand_%s_%2.2d")', bname[b], i+1)
         if ~ obry->HasVar(name) then begin
            message, /INFORM, string(FORMAT='(%"Adding variable %s")', name)
            obry->VarAdd, name, ['eta_rho','s_rho',time_name]
            obry->AttAdd, name, 'time', time_name
         endif
      endfor
   endfor

   ;; Load data

   bname = ['south','north']

   for b=0,1 do begin
      for i=0,n_dye-1 do begin
         name = string(FORMAT='(%"dye_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, OFFSET=[0,0,r], $
               replicate(dye_value[i], [grd_dim.xi_rho,grd_dim.s_rho])
         endfor
      endfor
      for i=0,n_mud-1 do begin
         name = string(FORMAT='(%"mud_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, OFFSET=[0,0,r], $
               replicate(mud_value[i], [grd_dim.xi_rho,grd_dim.s_rho])
         endfor
      endfor
      for i=0,n_sand-1 do begin
         name = string(FORMAT='(%"sand_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, OFFSET=[0,0,r], $
               replicate(sand_value[i], [grd_dim.xi_rho,grd_dim.s_rho])
           endfor
      endfor

   endfor

   bname = ['west','east']

   for b=0,1 do begin
      for i=0,n_dye-1 do begin
         name = string(FORMAT='(%"dye_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, OFFSET=[0,0,r], $
               replicate(dye_value[i], [grd_dim.eta_rho,grd_dim.s_rho])
         endfor
      endfor
      for i=0,n_mud-1 do begin
         name = string(FORMAT='(%"mud_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, rOFFSET=[0,0,r], $
               eplicate(mud_value[i], [grd_dim.eta_rho,grd_dim.s_rho])
         endfor
      endfor
      for i=0,n_sand-1 do begin
         name = string(FORMAT='(%"sand_%s_%2.2d")', bname[b], i+1)
         for r=0,n_time-1 do begin
            obry->VarPut, name, OFFSET=[0,0,r], $
               replicate(sand_value[i], [grd_dim.eta_rho,grd_dim.s_rho])
         endfor
      endfor

   endfor

   ;; Close the output file

   obj_destroy, obry

end

