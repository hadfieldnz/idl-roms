; svn $Id$
;+
; NAME:
;   MGH_ROMS_LAPLACIAN
;
; PURPOSE:
;   Calculate the Laplacian of a scalar on the ROMS grid. This scalar must be defined
;   on psi points.
;
; POSITIONAL PARAMETERS:
;   psi (input, numeric 2-D array)
;     Scalar, defined at psi points
;
;   pm (input, numeric 2-D array)
;     Grid metric in xi direction, defined at rho points.
;
;   pn (input, numeric 2-D array)
;     Grid metric in eta direction, defined at rho points.
;
; KEYWORD PARAMETERS:
;   mask (input, numeric 2-D array)
;     Sea mask (0=land, 1=sea) defined at rho points.
;
; RETURN VALUE:
;   The function returns a double-precision 2D array containing stream function
;   values, defined at psi points.
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
;   Mark Hadfield, 2003-08:
;     Written.
;-
function mgh_roms_laplacian, psi, pm, pn, MASK=mask

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if size(psi, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'psi'

   if size(pm, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pm'

   if size(pn, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pn'

   if size(psi, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'psi'

   if size(pm, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

   if size(pn, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

   dim = size(psi, /DIMENSIONS) + [1,1]

   d0 = dim[0]  &  d1 = dim[1]

   if ~ array_equal(size(pm, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'psi', 'pm'

   if ~ array_equal(size(pn, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'psi', 'pn'

   ;; Take first difference. Results are on u (dpsidx) & v (dpsidy) points:
   ;; interior values are calculated and extrapolated to the boundary

   pm_v = mgh_stagger(pm, DELTA=[0,-1])
   pn_u = mgh_stagger(pn, DELTA=[-1,0])

   dpsidx = fltarr(dim-[0,1])
   dpsidy = fltarr(dim-[1,0])

   dpsidx[1,0] = (psi[1:d0-2,*]-psi[0:d0-3,*]) * pm_v[1:d0-2,*]
   dpsidy[0,1] = (psi[*,1:d1-2]-psi[*,0:d1-3]) * pn_u[*,1:d1-2]

   dpsidx[0,*] = dpsidx[1,*]
   dpsidx[d0-1,*] = dpsidx[d0-2,*]

   dpsidy[*,0] = dpsidy[*,1]
   dpsidy[*,d1-1] = dpsidy[*,d1-2]

   ;; Second difference

   pm_p = mgh_stagger(pm, DELTA=[-1,-1])
   pn_p = mgh_stagger(pn, DELTA=[-1,-1])

   d2psidx2 = (dpsidx[1:d0-1,*]-dpsidx[0:d0-2,*]) * pm_p
   d2psidy2 = (dpsidy[*,1:d1-1]-dpsidy[*,0:d1-2]) * pn_p

   return, d2psidx2 + d2psidy2

end
