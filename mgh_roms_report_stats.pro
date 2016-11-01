; svn $Id$
;+
; NAME:
;   MGH_ROMS_REPORT_STATS
;
; PURPOSE:
;   For a specified variable in a ROMS file, work through the records
;   and report statistics of the variable.
;
; CALLING SEQUENCE:
;   mgh_report_stats, file, variable
;
; POSITIONAL PARAMETERS:
;   file (input, scalar string or object reference)
;     A reference to a ROMS history-file or station-file object or a
;     list of file names.
;
;   variable (input, scalar string)
;     The variable to be processed. Default is "zeta"
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2011-03:
;       Written.
;-

pro mgh_roms_report_stats, file, variable

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process positional parameters

   case size(file, /TNAME) of
      'STRING': begin
         case mgh_roms_file_type(file) of
            'ROMS station file': $
                 ofile = obj_new('MGHromsStation', file)
            else: $
                 ofile = obj_new('MGHromsHistory', file)
         endcase
         ofile_destroy = 1B
      end
      'OBJREF': begin
         ofile = file
         ofile_destroy = 0B
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'file'
   endcase

   if n_elements(variable) eq 0 then variable = 'zeta'

   dim = ofile->VarDimNames(variable)

   n_dim = n_elements(dim)

   has_unlim = ofile->DimInfo(dim[n_dim-1], /IS_UNLIMITED)

   if has_unlim then begin

      count = replicate(0, n_dim-1)
      offset = replicate(0, n_dim-1)

      n_rec = ofile->DimInfo(dim[n_dim-1], /DIMSIZE)

      vmin = fltarr(n_rec)
      vmax = fltarr(n_rec)

      for r=0,n_rec-1 do begin

         var = ofile->VarGet(variable, COUNT=[count,1], OFFSET=[offset,r], /AUTOSCALE)

         vmin[r] = min(var, /NAN)
         vmax[r] = max(var, /NAN)

         print, r, mgh_format_float([vmin[r],vmax[r]]), FORMaT='(%"%d: %s %s")'

      endfor

      print, mgh_format_float([min(vmin, /NAN),max(vmax, /NAN)]), FORMaT='(%"Overall: %s %s")'

   endif else begin

      count = replicate(0, n_dim)
      offset = replicate(0, n_dim)

      var = ofile->VarGet(variable, COUNT=count, OFFSET=offset, /AUTOSCALE)

      print, r, min(var, /NAN), max(var, /NAN)

   endelse




   if keyword_set(ofile_destroy) then obj_destroy, ofile

end

