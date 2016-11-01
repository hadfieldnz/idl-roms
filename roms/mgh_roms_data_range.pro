; svn $Id$
;+
; ROUTINE NAME:
;   MGH_ROMS_DATA_RANGE
;
; PURPOSE:
;   For a specified ROMS variable, this function returns a sensible
;   data range to be used as a default value.
;
; CALLING SEQUENCE:
;   result = mgh_roms_data_range(var)
;
; POSITIONAL PARAMETERS:
;   var (input, scalar string)
;     The name of a ROMS 2-D or 3-D variable.
;
; RETURN VALUE:
;   The function returns a 2-element, numeric vector equal to the
;   expected range of the variable in a typical ocean/coastal
;   simulation.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-11:
;     Written.
;   Mark Hadfield, 2003-12:
;     Superseded by MGH_ROMS_RESOLVE_DATA.
;-

function MGH_ROMS_DATA_RANGE, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   compile_opt OBSOLETE

   if size(var, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', var

   if n_elements(var) ne 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', var

   case var of
      'h': $
           result = [0,4000]
      'zeta': $
           result = [-1,1]
      'u': $
           result = [-1,1]*0.3
      'v': $
           result = [-1,1]*0.3
      'spd': $
           result = [0,0.4]
      'ubar': $
           result = [-1,1]*0.3
      'vbar': $
           result = [-1,1]*0.3
      'sbar': $
           result = [0,0.4]
      'w': $
           result = [-1,1]*1.E-3
      'omega': $
           result = [-1,1]*1.E-3
      'temp': $
           result = [0,25]
      'salt': $
           result = [34,36]
      'AKs': $
           result = [0,0.1]
      'AKt': $
           result = [0,0.1]
      'AKv': $
           result = [0,0.1]
      'rho': $
           result = [20.,50.]
      'NO3': $
           result = [0.,5.]
      'NH4': $
           result = [0.,1.]
      'TIC': $
           result = [2000.,2200.]
      'alkalinity': $
           result = [2200.,2400.]
      else: $
           result = [-1,1]
   endcase

   return, result

end
