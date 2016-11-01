; svn $Id$
;+
; NAME:
;   MGH_ROMS_STAGGER
;
; PURPOSE:
;   Given an array representing values on one of the ROMS horizontal grids
;   (rho, u, v, psi) interpolate to another such grid and return the result
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_STAGGER(a)
;
; POSITIONAL PARAMETERS:
;   A (input, 2D or 3D numeric array)
;     An array of values on FROM_GRID to be interpolated to TO_GRID
;
; KEYWORD PARAMETERS:
;   FROM_GRID (input, string)
;     The grid associated with the input values. Default is 'rho'.
;
;   TO_GRID (input, string)
;     The grid associated with the output values. Default is 'rho'.
;
; RETURN VALUE:
;   The function returns a numeric array with each dimension
;   in the input unchanged or expanded/contracted by 1.
;
; PROCEDURE:
;   Convert from input grid to RHO grid and then to output grid using MGH_STAGGER.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, Dec 2000:
;       Written.
;-
function MGH_ROMS_STAGGER, a, FROM_GRID=from_grid, TO_GRID=to_grid

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if size(a, /N_DIMENSIONS) lt 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'A'

   if n_elements(from_grid) eq 0 then from_grid = 'rho'

   if n_elements(to_grid) eq 0 then to_grid = 'rho'

   ;; Interpolate from initial grid to RHO grid

   delta = replicate(0, size(a, /N_DIMENSIONS))

   case strlowcase(from_grid) of
      'rho':
      'psi': delta[0:1] = [1,1]
      'u'  : delta[0:1] = [1,0]
      'v'  : delta[0:1] = [0,1]
   endcase

   a_rho = mgh_stagger(a, DELTA=delta)

   ;; Interpolate from RHO grid to final grid

   delta = replicate(0, size(a, /N_DIMENSIONS))

   case strlowcase(to_grid) of
      'rho':
      'psi': delta[0:1] = [-1,-1]
      'u'  : delta[0:1] = [-1,0]
      'v'  : delta[0:1] = [0,-1]
   endcase

   return, mgh_stagger(a_rho, DELTA=delta)

end
