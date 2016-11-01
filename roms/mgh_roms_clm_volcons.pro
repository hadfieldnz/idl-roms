; svn $Id$
;+
; NAME:
;   MGH_ROMS_CLM_VOLCONS
;
; PURPOSE:
;   Adjust velocities in a ROMS climate file (normal barotropic
;   velocities on boundary only) so that volume flux through
;   boundaries is zero.
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_VOLCONS, file_clm
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file to be modified.
;
; TO DO:
;   Allow selection of boundaries over which adjustment is to be done.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2005-09:
;     Written.
;-

pro mgh_roms_clm_volcons, file_clm

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

   ;; In future this procedure may be generalised to allow selection
   ;; of the boundaries over which the correction is made. Currently
   ;; the code to support this is incomplete, so we specify that all
   ;; boundaries will be adjusted.

   boundary = [0,1,2,3]

   ;; Open the file as a ROMS-history object: get grid data and
   ;; time-varying transport data and calculate imbalance vs time.

   message, /INFORM, string(FORMAT='(%"Reading data from climatology file %s")', file_clm)

   onc = obj_new('MGHromsHistory', file_clm)

   dim = [onc->DimInfo('xi_rho', /DIMSIZE), $
          onc->DimInfo('eta_rho', /DIMSIZE)]

   pm = onc->VarGet('pm')
   pn = onc->VarGet('pn')

   h = onc->VarGet('h')

   mask = onc->HasVar('mask_rho') ? onc->VarGet('mask_rho') : replicate(1, dim)

   ;; Assume that the time variable associated with the barotropic
   ;; velocities is specified via an attribute

   n_rec = onc->DimInfo(onc->AttGet('ubar', 'time'), /DIMSIZE)

   ;; Calculate correction, currently assuming that it is to be applied over all boundaries.

   ucorr = fltarr(n_rec)

   for r=0,n_rec-1 do begin

      tport = onc->GetTransportBox(RECORD=r)

      land = where(~finite(tport.depth), n_land)
      if n_land gt 0 then tport.depth[land] = 0
      mgh_undefine, land, n_land

      ucorr[r] = tport.transport[n_elements(tport.transport)-1] / $
                 total(tport.depth*mgh_diff(tport.distance))

   endfor

   obj_destroy, onc

   ;; Open the file again for modification

   message, /INFORM, string(FORMAT='(%"Opening climatology file %s")', file_clm)

   oclm = obj_new('MGHncFile', file_clm, /MODIFY)

   ;; Establish scaling factors for ubar and vbar data & whether values should be rounded
   ;; before writing.

   ubar_rnd = max(strmatch(['INT','SHORT','LONG'], oclm->VarInfo('ubar', /DATATYPE), /FOLD_CASE)) gt 0
   ubar_add = oclm->HasAtt('ubar', 'add_offset') ? oclm->AttGet('ubar', 'add_offset') : 0
   ubar_scl = oclm->HasAtt('ubar', 'scale_factor') ? oclm->AttGet('ubar', 'scale_factor') : 1

   vbar_rnd = max(strmatch(['INT','SHORT','LONG'], oclm->VarInfo('vbar', /DATATYPE), /FOLD_CASE)) gt 0
   vbar_add = oclm->HasAtt('vbar', 'add_offset') ? oclm->AttGet('vbar', 'add_offset') : 0
   vbar_scl = oclm->HasAtt('vbar', 'scale_factor') ? oclm->AttGet('vbar', 'scale_factor') : 1

   ;; Work through boundaries in turn, reading, adjusting & writing data.

   d0 = dim[0]
   d1 = dim[1]

   for b=0,n_elements(boundary)-1 do begin

      case boundary[b] of

         0: begin

            ;; Southern

            message, /INFORM, 'Adjusting vbar at southern boundary'
            for r=0,n_rec-1 do begin
               offset = [1,0,r]
               count = [d0-2,1,1]
               vbar = oclm->VarGet('vbar', OFFSET=offset, COUNT=count, /AUTOSCALE)
               vbar = (vbar-ucorr[r]-vbar_add)/vbar_scl
               if vbar_rnd then vbar = round(vbar)
               oclm->VarPut, 'vbar', vbar, OFFSET=offset, COUNT=count
            endfor

         end

         1: begin

            ;; Eastern

            message, /INFORM, 'Adjusting ubar at eastern boundary'
            for r=0,n_rec-1 do begin
               offset = [d0-2,1,r]
               count = [1,d1-2,1]
               ubar = oclm->VarGet('ubar', OFFSET=offset, COUNT=count, /AUTOSCALE)
               ubar = (ubar+ucorr[r]-ubar_add)/ubar_scl
               if ubar_rnd then ubar = round(ubar)
               oclm->VarPut, 'ubar', ubar, OFFSET=offset, COUNT=count
            endfor

         end

         2: begin

            ;; Northern

            message, /INFORM, 'Adjusting vbar at northern boundary'
            for r=0,n_rec-1 do begin
               offset = [1,d1-2,r]
               count = [d0-2,1,1]
               vbar = oclm->VarGet('vbar', OFFSET=offset, COUNT=count, /AUTOSCALE)
               vbar = (vbar+ucorr[r]-vbar_add)/vbar_scl
               if vbar_rnd then vbar = round(vbar)
               oclm->VarPut, 'vbar', vbar, OFFSET=offset, COUNT=count
            endfor

         end

         3: begin

            message, /INFORM, 'Adjusting ubar at western boundary'
            for r=0,n_rec-1 do begin
               offset = [0,1,r]
               count = [1,d1-2,1]
               ubar = oclm->VarGet('ubar', OFFSET=offset, COUNT=count, /AUTOSCALE)
               ubar = (ubar-ucorr[r]-ubar_add)/ubar_scl
               if ubar_rnd then ubar = round(ubar)
               oclm->VarPut, 'ubar', ubar, OFFSET=offset, COUNT=count
            endfor

         end

      endcase

   endfor

   ;; Close the output file

   obj_destroy, oclm

   mgh_toc

end

