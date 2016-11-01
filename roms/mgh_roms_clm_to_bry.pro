;+
; NAME:
;   MGH_ROMS_CLM_TO_BRY
;
; PURPOSE:
;   Given a ROMS climatology file, create a corresponding boundary file
;   and load data into it.
;
; CALLING SEQUENCE:
;   MGH_ROMS_CLM_TO_BRY, file_clm, file_bry
;
; POSITIONAL PARAMETERS:
;   file_clm (input, scalar string)
;     The name of a ROMS climatology file to be read.
;
;   file_bry (input, scalar string)
;     The name of a ROMS boundary file to be written.
;
; KEYWORD PARAMETERS:
;   PACK_DATA (input, switch)
;     If PACK_DATA is in effect, packed (and unpacked) variables in the
;     input file are copied to the output file with the packing (and
;     unpacking) intact. Otherwise all output variables are
;     unpacked. Default is off.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-01:
;     Written.
;   Mark Hadfield, 2004-01:
;     Modified for new conventions for time-related dimensions & variables.
;   Mark Hadfield, 2007-07:
;     Revised code for specifying names of variables to be copied from
;     the climatology file to the boundary file.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;   Mark Hadfield, 2010-02:
;     PACK_DATA functionality added.
;   Mark Hadfield, 2011-08:
;     Now reports progress at each boundary.
;   Mark Hadfield, 2011-11:
;     Variable handling overhauled.
;-
function mgh_roms_clm_to_bry_match, clm_var, var, COUNT=count

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Utility function for MGH_ROMS_CLM_TO_BRY: given the list of
   ;; variable names from a climatology file and a variable-name pattern,
   ;; return the list of matches along with the number via the COUNT keyword.

   l_match = where(strmatch(clm_var, var), count)

   return, count gt 0 ? clm_var[l_match] : -1

end

function mgh_roms_clm_to_bry_fmt, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Utility function for MGH_ROMS_CLM_TO_BRY: given the name
   ;; of a variable in the climatology file, return a C-style
   ;; format string to be used in constructing the corresponding
   ;; variable names the boundary file.

   case 1B of
      strmatch(var, 'dye*'): begin
         result = mgh_str_subst(var, 'dye', 'dye_%s')
      end
      else: result = var+'_%s'
   endcase

   return, '(%"'+result+'")'

end

pro mgh_roms_clm_to_bry, file_clm, file_bry, $
     PACK_DATA=pack_data, VARIABLES=variables

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check file name etc

   if n_elements(file_clm) eq 0 then $
        message, 'Name for climatology file not supplied'

   if n_elements(file_bry) eq 0 then $
        message, 'Name for boundary file not supplied'

   ;; Specify variables whose boundary values are to copied to the
   ;; boundary file, if they exist in the climatology file. Variable
   ;; names in the boundary file are formed from the corresponding
   ;; names in the climatology file by adding "_west", _north",
   ;; etc. Note that the standard ROMS varinfo.dat file does not
   ;; follow this convention consistently, which is a problem.
   ;; I have used a modified version of the varinfo.dat to work
   ;; around it, but I'm not sure if the modifications are still
   ;; there.

   if n_elements(variables) eq 0 then $
        variables = ['zeta','ubar','vbar','u','v','temp','salt','dye_??']

   n_var = n_elements(variables)

   ;; Open the files

   if n_elements(file_clm) gt 1 then begin
      message, /INFORM, string(FORMAT='(%"Opening climate file sequence beginning ' + $
                               'with %s")', file_clm[0])
      oclm = obj_new('MGHncSequence', file_clm)
   endif else begin
      message, /INFORM, string(FORMAT='(%"Opening climate file %s")', file_clm)
      oclm = obj_new('MGHncReadFile', file_clm)
   endelse

   message, /INFORM, string(FORMAT='(%"Creating boundary file %s")', file_bry)

   obry = obj_new('MGHncFile', file_bry, /CREATE, /CLOBBER)

   ;; Write some global attributes

   obry->AttAdd, /GLOBAL, 'type', "Boundary forcing file"

   ;; Copy dimensions

   message, /INFORM, 'Copying dimensions to boundary file'

   dims = oclm->DimNames()

   dpat = ['xi_rho','eta_rho','xi_u','eta_u','xi_v','eta_v','s_rho','*time*']

   for i=0,n_elements(dims)-1 do begin
      for j=0,n_elements(dpat)-1 do begin
         if strmatch(dims[i], dpat[j]) then obry->DimCopy, oclm, dims[i]
      endfor
   endfor

   ;; Specify time variables

   clm_var = oclm->VarNames()

   l_tvar = where(strmatch(clm_var, '*time*'), n_tvar)

   tvar = n_tvar gt 0 ? clm_var[l_tvar] : -1

   ;; Copy variables

   message, /INFORM, 'Copying time variables to boundary file'

   for i=0,n_tvar-1 do begin
      obry->VarCopy, oclm, tvar[i], /DEFINITIONS
      obry->VarCopy, oclm, tvar[i], /ATTRIBUTES
   endfor

   message, /INFORM, 'Defining boundary variables in boundary file'

   if keyword_set(pack_data) then begin

      ;; If PACK_DATA is in effect, scaling attributes for each
      ;; variable are copied to the output file and data will be
      ;; copied over without unpacking

      autoscale = 0B

      for i=0,n_var-1 do begin

         match = mgh_roms_clm_to_bry_match(clm_var, variables[i], $
                                           COUNT=n_match)

         for j=0,n_match-1 do begin

            vdims = oclm->VarDimNames(match[j], COUNT=n_vdims)
            vtype = oclm->VarInfo(match[j], /DATATYPE)
            vtatt = oclm->AttGet(match[j], 'time')
            if oclm->HasAtt(match[j], 'add_offset') then $
                 vadd = oclm->AttGet(match[j], 'add_offset')
            if oclm->HasAtt(match[j], 'scale_factor') then $
                 vscl = oclm->AttGet(match[j], 'scale_factor')

            fmt = mgh_roms_clm_to_bry_fmt(match[j])

            obry->VarAdd, string(FORMAT=fmt, 'west'), $
                          [vdims[1],vdims[2:n_vdims-1]], $
                          _STRICT_EXTRA=create_struct(vtype, 1B)
            obry->AttAdd, string(FORMAT=fmt, 'west'), 'time', vtatt
            if n_elements(vadd) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'west'), 'add_offset', vadd
            if n_elements(vscl) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'west'), 'scale_factor', vscl

            obry->VarAdd, string(FORMAT=fmt, 'east'), [vdims[1],vdims[2:n_vdims-1]], $
                          _STRICT_EXTRA=create_struct(vtype, 1B)
            obry->AttAdd, string(FORMAT=fmt, 'east'), 'time', vtatt
            if n_elements(vadd) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'east'), 'add_offset', vadd
            if n_elements(vscl) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'east'), 'scale_factor', vscl

            obry->VarAdd, string(FORMAT=fmt, 'south'), [vdims[0],vdims[2:n_vdims-1]], $
                          _STRICT_EXTRA=create_struct(vtype, 1B)
            obry->AttAdd, string(FORMAT=fmt, 'south'), 'time', vtatt
            if n_elements(vadd) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'south'), 'add_offset', vadd
            if n_elements(vscl) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'south'), 'scale_factor', vscl

            obry->VarAdd, string(FORMAT=fmt, 'north'), [vdims[0],vdims[2:n_vdims-1]], $
                          _STRICT_EXTRA=create_struct(vtype, 1B)
            obry->AttAdd, string(FORMAT=fmt, 'north'), 'time', vtatt
            if n_elements(vadd) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'north'), 'add_offset', vadd
            if n_elements(vscl) gt 0 then $
                 obry->AttAdd, string(FORMAT=fmt, 'north'), 'scale_factor', vscl

            mgh_undefine, vadd, vscl

         endfor

      endfor

   endif else begin

      ;; If PACK_DATA is not in effect, output variables are all
      ;; unscaled floating point and data will be unpacked
      ;; automatically when read.

      autoscale = 1B

      for i=0,n_var-1 do begin

         match = mgh_roms_clm_to_bry_match(clm_var, variables[i], $
                                           COUNT=n_match)

         for j=0,n_match-1 do begin

            vdims = oclm->VarDimNames(match[j], COUNT=n_vdims)
            vtatt = oclm->AttGet(match[j], 'time')

            fmt = mgh_roms_clm_to_bry_fmt(match[j])

            obry->VarAdd, string(FORMAT=fmt, 'west'), [vdims[1],vdims[2:n_vdims-1]]
            obry->AttAdd, string(FORMAT=fmt, 'west'), 'time', vtatt

            obry->VarAdd, string(FORMAT=fmt, 'east'), [vdims[1],vdims[2:n_vdims-1]]
            obry->AttAdd, string(FORMAT=fmt, 'east'), 'time', vtatt

            obry->VarAdd, string(FORMAT=fmt, 'south'), [vdims[0],vdims[2:n_vdims-1]]
            obry->AttAdd, string(FORMAT=fmt, 'south'), 'time', vtatt

            obry->VarAdd, string(FORMAT=fmt, 'north'), [vdims[0],vdims[2:n_vdims-1]]
            obry->AttAdd, string(FORMAT=fmt, 'north'), 'time', vtatt

         endfor

      endfor

   endelse

   message, /INFORM, 'Loading time variables to boundary file'

   for i=0,n_tvar-1 do begin
      obry->VarCopy, oclm, tvar[i], /DATA
   endfor

   message, /INFORM, 'Loading boundary variables...'

   for i=0,n_var-1 do begin

      match = mgh_roms_clm_to_bry_match(clm_var, variables[i], COUNT=n_match)

      for j=0,n_match-1 do begin

         message, /INFORM, '...'+match[j]

         n_vdims = oclm->VarInfo(match[j], /N_DIMS)

         nn = replicate(0, n_vdims-2)

         fmt = mgh_roms_clm_to_bry_fmt(match[j])

         obry->VarPut, string(FORMAT=fmt, 'west'), $
                       reform(oclm->VarGet(match[j], AUTOSCALE=autoscale, $
                                           COUNT=[1,0,nn], OFFSET=[0,0,nn]))
         obry->VarPut, string(FORMAT=fmt, 'east'), $
                       reform(oclm->VarGet(match[j], AUTOSCALE=autoscale, $
                                           COUNT=[1,0,nn], OFFSET=[-1,0,nn]))
         obry->VarPut, string(FORMAT=fmt, 'south'), $
                       reform(oclm->VarGet(match[j], AUTOSCALE=autoscale, $
                                           COUNT=[0,1,nn], OFFSET=[0,0,nn]))
         obry->VarPut, string(FORMAT=fmt, 'north'), $
                       reform(oclm->VarGet(match[j], AUTOSCALE=autoscale, $
                                           COUNT=[0,1,nn], OFFSET=[0,-1,nn]))
      endfor

   endfor

   obj_destroy, [oclm,obry]

end

