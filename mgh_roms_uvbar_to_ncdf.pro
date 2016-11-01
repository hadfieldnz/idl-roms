; svn $Id$
;+
; NAME:
;   MGH_ROMS_UVBAR_TO_NCDF
;
; PURPOSE:
;   Get selected records of ubar & vbar from a ROMS history or averages
;   file, and write them to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   mgh_roms_uvbar_to_ncdf, history, file_out
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
;   RECORDS (input, integer vector)
;     Records to be copied
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-11:
;     Written.
;-

pro mgh_roms_uvbar_to_ncdf, history, file_out, RECORDS=records

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         history_file = ohis
         history_destroy = 1B
      end
      'OBJREF': begin
         ohis = history
         history_file = history
         history_destroy = 0B
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
        message, 'Name for out file not supplied'

   ;; Select records to be transferred

   n_time = ohis->DimInfo('ocean_time', /DIMSIZE)

   if n_elements(records) eq 0 then records = lindgen(n_time)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   ;; A bit of header info

   fmt = '(%"Snapshots of depth-integrated velocity from ROMS dataset %s")'
   oout->AttAdd, /GLOBAL, 'long_name', string(mgh_get_property(ohis, /NAME), FORMAT=fmt)

   ;; Create dimensions. The grid is the ROMS interior rho grid

   dim = [ohis->DimInfo('xi_rho', /DIMSIZE)-2, $
          ohis->DimInfo('eta_rho', /DIMSIZE)-2]

   oout->DimAdd, 'xi', dim[0]
   oout->DimAdd, 'eta', dim[1]

   oout->DimAdd, 'time'

   ;; Create grid variables

   oout->VarAdd, 'lon', ['xi','eta'], /FLOAT
   oout->AttAdd, 'lon', 'long_name', 'longitude'
   oout->AttAdd, 'lon', 'units', 'degree_east'

   oout->VarAdd, 'lat', ['xi','eta'], /FLOAT
   oout->AttAdd, 'lat', 'long_name', 'latitude'
   oout->AttAdd, 'lat', 'units', 'degree_north'

   oout->VarAdd, 'mask', ['xi','eta'], /BYTE
   oout->AttAdd, 'mask', 'long_name', 'mask'
   oout->AttAdd, 'mask', 'option_0', 'land'
   oout->AttAdd, 'mask', 'option_1', 'water'

   oout->VarAdd, 'h', ['xi','eta'], /FLOAT
   oout->AttAdd, 'h', 'long_name', 'bathymetry'
   oout->AttAdd, 'h', 'units', 'meter'

   oout->VarAdd, 'angle', ['xi','eta'], /FLOAT
   oout->AttAdd, 'angle', 'long_name', 'angle between XI-axis and EAST'
   oout->AttAdd, 'angle', 'units', 'radian'

   ;; Create time variable

   oout->VarAdd, 'time', ['time'], /DOUBLE
   oout->AttAdd, 'time', 'long_name', 'time since initialization'
   oout->AttAdd, 'time', 'units', 'day'

   ;; Create result variables

   fill_real = 1.E10

   oout->VarAdd, 'ubar', ['xi','eta','time'], /FLOAT
   oout->AttAdd, 'ubar', 'long_name', 'vertically integrated u-velocity component'
   oout->AttAdd, 'ubar', 'units', 'meter second-1'
   oout->AttAdd, 'ubar', '_FillValue', fill_real

   oout->VarAdd, 'vbar', ['xi','eta','time'], /FLOAT
   oout->AttAdd, 'vbar', 'long_name', 'vertically integrated v-velocity component'
   oout->AttAdd, 'vbar', 'units', 'meter second-1'
   oout->AttAdd, 'vbar', '_FillValue', fill_real

   ;; Add grid data

   lon = ohis->VarGet('lon_rho', OFFSET=[1,1], COUNT=dim)
   oout->VarPut, 'lon', lon

   lat = ohis->VarGet('lat_rho', OFFSET=[1,1], COUNT=dim)
   oout->VarPut, 'lat', lat

   h = ohis->VarGet('h', OFFSET=[1,1], COUNT=dim)
   oout->VarPut, 'h', h

   mask = round(ohis->VarGet('mask_rho', OFFSET=[1,1], COUNT=dim))
   oout->VarPut, 'mask', mask

   angle = ohis->VarGet('angle', OFFSET=[1,1], COUNT=dim)
   oout->VarPut, 'angle', angle

   ;; Work through records, transferring data

   for r=0,n_elements(records)-1 do begin

      rec = records[r]

      time = ohis->VarGet('ocean_time', OFFSET=[rec], COUNT=[1], /AUTOSCALE)

      ubar = ohis->VarGet('ubar', /AUTOSCALE, OFFSET=[0,1,rec], COUNT=[0,dim[1],1])
      vbar = ohis->VarGet('vbar', /AUTOSCALE, OFFSET=[1,0,rec], COUNT=[dim[0],0,1])

      l_miss = where(~ finite(ubar), n_miss)
      if n_miss gt 0 then ubar[l_miss] = 0
      mgh_undefine, l_miss, n_miss

      l_miss = where(~ finite(vbar), n_miss)
      if n_miss gt 0 then vbar[l_miss] = 0
      mgh_undefine, l_miss, n_miss

      oout->VarPut, 'time', temporary(time), OFFSET=[r]

      oout->VarPut, 'ubar', mgh_stagger(ubar, DELTA=[-1,0]), OFFSET=[0,0,r]
      oout->VarPut, 'vbar', mgh_stagger(vbar, DELTA=[0,-1]), OFFSET=[0,0,r]

   endfor

   ;; Clean up

   obj_destroy, oout
   if history_destroy then obj_destroy, ohis

end

