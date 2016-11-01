; svn $Id$
;+
; NAME:
;   MGH_ROMS_GSVEL
;
; PURPOSE:
;   Calculate geostrophic surface velocity from a ROMS SSH field.
;
; POSITIONAL PARAMETERS:
;   zeta (input, numeric 2-D array)
;     Sea surface height, defined at rho points
;
;   pp (input, numeric 2-D array)
;     Grid metric coefficient pn (DIRECTION=0) or pm (DIRECTION=1), defined at rho points.
;
;   f (input, numeric 2-D array)
;     Coriolis parameter, defined at rho points
;
; KEYWORD PARAMETERS:
;   DIRECTION (input, integer scalar, default = 0)
;     Direction of velocity component to be calculated: 0 = xi
;     component and 1 = eta component.
;
;   MASK (input, numeric 2-D array)
;     Sea mask (0=land, 1=sea) defined at rho points.
;
; RETURN VALUE:
;   The function returns a 2D array containing geostrophic velocity
;   components at v (DIRECTION=0) or u points (DIRECTION=1)
;
; PROCEDURE:
;   Straightforward differencing.
;
;###########################################################################
;
; This software is provided subject to the following conditions:
;
; 1.  NIWA makes no representations or warranties regarding the
;     accuracy of the software, the use to which the software may
;     be put or the results to be obtained from the use of the
;     software.  Accordingly NIWA accepts no liability for any loss
;     or damage (whether direct of indirect) incurred by any person
;     through the use of or reliance on the software.
;
; 2.  NIWA is to be acknowledged as the original author of the
;     software where the software is used or presented in any form.
;
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2005-08:
;     Written.
;-
function mgh_roms_gsvel, zeta, p, f, $
     DIRECTION=direction, MASK=mask

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if size(zeta, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'zeta'

   if size(p, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'p'

   if size(f, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'f'

   if size(zeta, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'zeta'

   if size(p, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'p'

   if size(f, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'f'

   dim = size(zeta, /DIMENSIONS)

   d0 = dim[0]  &  d1 = dim[1]

   if ~ array_equal(size(p, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'zeta', 'p'

   if ~ array_equal(size(f, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'zeta', 'f'

   ;; Process DIRECTION keyword then specify constants for constructing the
   ;; finite-difference expression

   if n_elements(direction) eq 0 then direction = 0

   case direction of
      0: begin
         gdel = [0,-1]
         gdim = 2
         gmult = -1
      end
      1: begin
         gdel = [-1,0]
         gdim = 1
         gmult = 1
      end
   end

   ;; Calculate geostrophic surface velocity

   grav = 9.81

   result = gmult * mgh_diff(zeta, gdim) * grav * $
            mgh_stagger(p, DELTA=gdel) / mgh_stagger(f, DELTA=gdel)

   ;; Apply mask, if specified

   if n_elements(mask) gt 0 then begin

      if ~ array_equal(size(mask, /DIMENSIONS), dim) then $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'zeta', 'mask'

;      mask_psi = fix(round(mask[0:d0-2,0:d1-2] * mask[1:d0-1,0:d1-2] * $
;                           mask[0:d0-2,1:d1-1] * mask[1:d0-1,1:d1-1]))

      if min(mask) eq 0 then begin

         case direction of
            0: gmask = fix(round(mask[*,0:d1-2] * mask[*,1:d1-1]))
            1: gmask = fix(round(mask[0:d0-2,*] * mask[1:d0-1,*]))
         endcase

         result[where(gmask eq 0)] = !values.f_nan

      endif

   endif

   return, result

end
