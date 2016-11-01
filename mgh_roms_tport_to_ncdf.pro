; svn $Id$
;+
; NAME:
;   MGH_ROMS_TPORT_TO_NCDF
;
; PURPOSE:
;   Get boundary transport data from a ROMS file and save it to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_TPORT, history, file_out
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   file_out (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-04:
;     Written.
;-

pro MGH_ROMS_TPORT_TO_NCDF, history, file_out, RECORD=record

   compile_opt DEFINT32
   compile_opt STRICTARR

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
        message, 'Name for climate file not supplied'

   ;; Get transport data

   tport = ohis->GetTransportBox(RECORD=record)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   oout->AttAdd, /GLOBAL, 'long_name', 'Boundary transport data from ROMS dataset '+mgh_get_property(ohis, /NAME)

   ;; Create dimensions

   oout->DimAdd, 's_psi', n_elements(tport.distance)
   oout->DimAdd, 's_rho', n_elements(tport.depth)

   ;; Create variables

   oout->VarAdd, 'distance', 's_psi'
   oout->AttAdd, 'distance', 'long_name', 'distance along boundary'
   oout->AttAdd, 'distance', 'units', 'm'

   oout->VarAdd, 'transport', 's_psi'
   oout->AttAdd, 'transport', 'long_name', 'inward transport'
   oout->AttAdd, 'transport', 'units', 'm^3 s^-1'

   oout->VarAdd, 'depth', 's_rho'
   oout->AttAdd, 'depth', 'long_name', 'depth'
   oout->AttAdd, 'depth', 'units', 'm'

   oout->VarAdd, 'lon', 's_rho'
   oout->AttAdd, 'lon', 'long_name', 'longitude'
   oout->AttAdd, 'lon', 'units', 'degrees E'

   oout->VarAdd, 'lat', 's_rho'
   oout->AttAdd, 'lat', 'long_name', 'latitude'
   oout->AttAdd, 'lat', 'units', 'degrees N'

   oout->VarAdd, 'bound', 's_rho', /SHORT
   oout->AttAdd, 'bound', 'long_name', 'boundary number'
   oout->AttAdd, 'bound', 'valid_range', [1,4]

   ;; Write data

   oout->VarPut, 'distance', tport.distance
   oout->VarPut, 'transport', tport.transport

   oout->VarPut, 'depth', tport.depth
   oout->VarPut, 'lon', tport.lon
   oout->VarPut, 'lat', tport.lat
   oout->VarPut, 'bound', tport.boundary

   ;; Clean up

   obj_destroy, oout

   if history_destroy then obj_destroy, ohis

end

