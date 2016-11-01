; svn $Id$
;+
; NAME:
;   MGH_ROMS_CLM_TO_INI
;
; PURPOSE:
;   Given a ROMS climatology file, create an initialisation (restart)
;   file from it. Missing fields are filled with 0s.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_TO_INI, file_clm, file_ini
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file to act as a template (read-only)
;
;   file_ini (input, scalar string)
;     The name of a ROMS initialisation file to be created (write-only)
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-03:
;     Written.
;-

pro MGH_ROMS_CLM_TO_INI, file_clm, file_ini

   compile_opt DEFINT32
   compile_opt STRICTARR

   ;; Check file name etc

   if n_elements(file_clm) eq 0 then $
        message, 'Name for grid file not supplied'

   if not file_test(file_clm, /READ) then $
        message, 'Grid file cannot be read'

   if n_elements(file_ini) eq 0 then $
        message, 'Name for climate file not supplied'

   ;; Open the files

   message, /INFORM, string(FORMAT='(%"Opening climatology file %s")', file_clm)

   oclm = obj_new('MGHromsHistory', file_clm)

   message, /INFORM, string(FORMAT='(%"Creating initialisation file %s")', file_ini)

   oini = obj_new('MGHncFile', file_ini, /CREATE, /CLOBBER)

   ;; Get grid dimensions

   dim = oclm->DimRho()

   ;; Specify some global attributes

   oini->AttAdd, /GLOBAL, 'type', 'ROMS restart file'

   fmt = '(%"Created by procedure mgh_roms_grd_to_clm from %s at %s")'
   oini->AttAdd, /GLOBAL, 'history', $
                 string(FORMAT=fmt,file_clm,mgh_dt_string(mgh_dt_now()))

   ;; Copy dimensions

   message, /INFORM, 'Copying dimensions to initialisation file'

   oini->DimCopy, oclm

   ;; Copy variables

   vars = oclm->VarNames()

   ;; Filter variable list here??

   oini->VarCopy, oclm, /DEFINITION, vars
   oini->VarCopy, oclm, /ATTRIBUTES, vars

   if max(strmatch(vars, 'zeta')) eq 0 then begin
      oini->VarAdd, 'zeta', ['xi_rho','eta_rho','time']
   endif

   if max(strmatch(vars, 'ubar')) eq 0 then begin
      oini->VarAdd, 'ubar', ['xi_u','eta_u','time']
   endif

   if max(strmatch(vars, 'vbar')) eq 0 then begin
      oini->VarAdd, 'vbar', ['xi_v','eta_v','time']
   endif

   if max(strmatch(vars, 'u')) eq 0 then begin
      oini->VarAdd, 'u', ['xi_u','eta_u','s_rho','time']
   endif

   if max(strmatch(vars, 'v')) eq 0 then begin
      oini->VarAdd, 'v', ['xi_v','eta_v','s_rho','time']
   endif

   if max(strmatch(vars, 'temp')) eq 0 then begin
      oini->VarAdd, 'temp', ['xi_rho','eta_rho','s_rho','time']
   endif

   if max(strmatch(vars, 'salt')) eq 0 then begin
      oini->VarAdd, 'salt', ['xi_rho','eta_rho','s_rho','time']
   endif

   ;; Copy data, filling fields with 0s where necessary

   oini->VarCopy, oclm, /DATA, vars

   if max(strmatch(vars, 'zeta')) eq 0 then begin
      oini->VarPut, 'zeta', fltarr(dim[0:1])
   endif

   if max(strmatch(vars, 'ubar')) eq 0 then begin
      oini->VarPut, 'ubar', fltarr(dim[0:1]+[-1,0])
   endif

   if max(strmatch(vars, 'vbar')) eq 0 then begin
      oini->VarPut, 'vbar', fltarr(dim[0:1]+[0,-1])
   endif

   if max(strmatch(vars, 'u')) eq 0 then begin
      oini->VarPut, 'u', fltarr(dim+[-1,0,0])
   endif

   if max(strmatch(vars, 'v')) eq 0 then begin
      oini->VarPut, 'v', fltarr(dim+[0,-1,0])
   endif

   if max(strmatch(vars, 'temp')) eq 0 then begin
      oini->VarPut, 'temp', fltarr(dim)
   endif

   if max(strmatch(vars, 'salt')) eq 0 then begin
      oini->VarPut, 'salt', fltarr(dim)
   endif

   ;; Clean up

   obj_destroy, oini

   obj_destroy, oclm

end

