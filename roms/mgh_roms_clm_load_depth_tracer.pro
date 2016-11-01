;+
; NAME:
;   MGH_ROMS_CLM_LOAD_DEPTH_TRACER
;
; PURPOSE:
;   Load a passive tracer (dye) representing depth to a ROMS climatology file.
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_LOAD_DEPTH_TRACER, file_clm
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file into which data are to
;     be written.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-03:
;     Written.
;-
pro mgh_roms_clm_load_depth_tracer, file_clm, $
     TIME_NAME=time_name, VARIABLE=variable

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

   if n_elements(variable) eq 0 then variable = 'dye_01'

   if n_elements(time_name) eq 0 then time_name = 'time_fix'

   ;; Open climatology file as MGHromsHistory object to get dimensions &
   ;; grid data. A bit messy this!

   fmt = '(%"Opening climatology file %s for read-only access")'
   message, /INFORM, string(FORMAT=fmt, file_clm)

   oclm = obj_new('MGHromsHistory', file_clm)

   grid = {dim: oclm->DimRho(), $
           z: oclm->GetZGrid()}

   obj_destroy, oclm

   ;; Open the file

   fmt = '(%"Opening climatology file %s for modification")'
   message, /INFORM, string(FORMAT=fmt, file_clm)

   oclm = obj_new('MGHncFile', file_clm, /MODIFY)

   ;; Check that the required time dimension is available

   if ~ oclm->HasDim(time_name) then message, 'Time dimension missing'

   ;; Define variable

   message, /INFORM, string(FORMAT='(%"Adding variable %s")', variable)
   oclm->VarAdd, variable, ['xi_rho','eta_rho','s_rho',time_name]
   oclm->AttAdd, variable, 'time', time_name

   ;; Load data

   oclm->VarPut, variable, -grid.z

   ;; Close the output file

   obj_destroy, oclm

end
