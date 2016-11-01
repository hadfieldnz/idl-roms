; svn $Id$
;+
; NAME:
;   MGH_ROMS_GRD_SUBSET
;
; PURPOSE:
;   Create a ROMS grid file containing a subset of another grid.
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_SUBSET, file_grd, file_src
;
; POSITIONAL PARAMETERS:
;   file_grd (input, scalar string)
;     The name of a ROMS grid file to be generated (write-only)
;
;   file_src (input, scalar string)
;     The name of the ROMS grid file that the subset is to be taken from (read-only)
;
; SIDE EFFECTS:
;   Writes (if necessary overwrites) the grid file.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2004-02:
;     Written
;-

pro mgh_roms_grd_subset, file_grd, file_src, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(file_grd) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file_grd'

   if n_elements(file_grd) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'file_grd'

   if n_elements(file_src) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file_src'

   if n_elements(file_src) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'file_src'

   if ~ file_test(file_src, /READ) then $
        message, 'Source file cannot be read'

   message, /INFORM, 'Creating grid file '+strtrim(file_grd,2)

   ogrd = obj_new('MGHncFile', file_grd, /CREATE, /CLOBBER)

   ogrd->AttAdd, /GLOBAL, 'type', "ROMS grid file"

   message, /INFORM, 'Reading data from source file '+strtrim(file_src,2)

   osrc = obj_new('MGHncReadFile', file_src)

   if n_elements(xi_range) eq 0 then $
        xi_range = [0,ogrd->DimInfo('xi_rho', /DIMSIZE)-1]

   if n_elements(eta_range) eq 0 then $
        eta_range = [0,ogrd->DimInfo('eta_rho', /DIMSIZE)-1]

   msg = ['Extracting subset',strtrim(xi_range, 2),strtrim(eta_range, 2)]

   message, /INFORM, strjoin(temporary(msg), ' ')

   ;; Define dimensions

   n_rho = [xi_range[1]-xi_range[0]+1,eta_range[1]-eta_range[0]+1]

   ogrd->DimAdd, 'xi_rho', n_rho[0]
   ogrd->DimAdd, 'xi_u', n_rho[0]-1
   ogrd->DimAdd, 'xi_v', n_rho[0]
   ogrd->DimAdd, 'xi_psi', n_rho[0]-1

   ogrd->DimAdd, 'eta_rho', n_rho[1]
   ogrd->DimAdd, 'eta_u', n_rho[1]
   ogrd->DimAdd, 'eta_v', n_rho[1]-1
   ogrd->DimAdd, 'eta_psi', n_rho[1]-1

   ;; Define variables & copy data

   osrc->GetProperty, VAR_NAMES=vars

   for i=0,n_elements(vars)-1 do begin

      var = vars[i]

      if var eq 'hraw' then continue

      dims = osrc->VarDimNames(var, COUNT=n_dims)

      osrc->VarInfo, var, DATATYPE=nctype
      ogrd->VarAdd, var, dims, NCTYPE=nctype

   endfor

   for i=0,n_elements(vars)-1 do begin

      var = vars[i]

      if var eq 'hraw' then continue

      dims = osrc->VarDimNames(var, COUNT=n_dims)

      case 1B of

         n_dims eq 0: begin
            ogrd->VarPut, var, osrc->VarGet(var)
         end

         n_dims eq 2 && dims[0] eq 'xi_rho' && dims[1] eq 'eta_rho': begin
            offset = [xi_range[0],eta_range[0]]
            count = [xi_range[1]-xi_range[0]+1,eta_range[1]-eta_range[0]+1]
            ogrd->VarPut, var, osrc->VarGet(var, OFFSET=offset, COUNT=count)
         end

         n_dims eq 2 && dims[0] eq 'xi_u' && dims[1] eq 'eta_u': begin
            offset = [xi_range[0],eta_range[0]]
            count = [xi_range[1]-xi_range[0],eta_range[1]-eta_range[0]+1]
            ogrd->VarPut, var, osrc->VarGet(var, OFFSET=offset, COUNT=count)
         end

         n_dims eq 2 && dims[0] eq 'xi_v' && dims[1] eq 'eta_v': begin
            offset = [xi_range[0],eta_range[0]]
            count = [xi_range[1]-xi_range[0]+1,eta_range[1]-eta_range[0]]
            ogrd->VarPut, var, osrc->VarGet(var, OFFSET=offset, COUNT=count)
         end

         n_dims eq 2 && dims[0] eq 'xi_psi' && dims[1] eq 'eta_psi': begin
            offset = [xi_range[0],eta_range[0]]
            count = [xi_range[1]-xi_range[0],eta_range[1]-eta_range[0]]
            ogrd->VarPut, var, osrc->VarGet(var, OFFSET=offset, COUNT=count)
         end

      endcase

   endfor

   goto, skipit

   v2d_rho = ['pm','pn','x_rho','y_rho','lon_rho','lat_rho','dmde','dndx','angle','f','h','mask_rho']
   v2d_u = ['mask_u']
   v2d_v = ['mask_v']
   v2d_psi = ['mask_psi']

   for i=0,n_elements(v2d)-1 do begin
      if osrc->HasVar(v2d[i]) then ogrd->VarAdd, v2d[i], ['xi_rho','eta_rho']
   endfor

   ;; Copy data

   ogrd->VarPut, 'spherical', osrc->VarGet('spherical')

   ogrd->VarPut, 'xl', 0.
   ogrd->VarPut, 'el', 0.

skipit:

   ;; Close files

   obj_destroy, ogrd

   obj_destroy, osrc

end

