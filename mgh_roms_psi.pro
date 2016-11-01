; svn $Id$
;+
; NAME:
;   MGH_ROMS_PSI
;
; PURPOSE:
;   Calculate stream function for a ROMS vector field
;
; POSITIONAL PARAMETERS:
;   u (input, numeric 2-D array)
;     Vector component in xi direction, defined at u points.
;
;   v (input, numeric 2-D array)
;     Vector component in eta direction, defined at v points
;
;   pm (input, numeric 2-D array)
;     Grid metric in xi direction, defined at rho points.
;
;   pn (input, numeric 2-D array)
;     Grid metric in eta direction, defined at rho points.
;
;   Note that, to calculate stream function at all psi points,
;   all U & v points are required, including the tangential
;   boundary values that are not used in the computation.
;   (This is an anomaly, which I hope to correct some time.)
;
; KEYWORD PARAMETERS:
;   CONVERGENCE (output, numeric scalar)
;     Imbalance in stream function integrated around the boundary in
;     m^3 s^-1.
;
;   MASK (input, numeric 2-D array)
;     Sea mask (0=land, 1=sea) defined at rho points.
;
; RETURN VALUE:
;   The function returns a double-precision 2D array containing stream function
;   values, defined at psi points.
;
; PROCEDURE:
;   Naive integration, with allowance for weak divergence in the vector field.
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
;   Mark Hadfield, 2003-07:
;     Written.
;-
function mgh_roms_psi, u, v, pm, pn, $
     CONVERGENCE=convergence, MASK=mask

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if size(u, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'u'

   if size(v, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'v'

   if size(pm, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pm'

   if size(pn, /N_ELEMENTS) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pn'

   if size(u, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'u'

   if size(v, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'v'

   if size(pm, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

   if size(pn, /N_DIMENSIONS) ne 2 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

   dim = size(u, /DIMENSIONS) + [1,0]

   d0 = dim[0]  &  d1 = dim[1]

   if ~ array_equal(size(v, /DIMENSIONS), dim-[0,1]) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u', 'v'

   if ~ array_equal(size(pm, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u & v', 'pm'

   if ~ array_equal(size(pn, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u & v', 'pn'

   ;; Process keyword parameters

   if n_elements(mask) eq 0 then mask = replicate(1, dim)

   if ~ array_equal(size(mask, /DIMENSIONS), dim) then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u & v', 'mask'

   ;; Set up arrays defining the grid

   pm_v = mgh_stagger(pm, DELTA=[0,-1])
   pn_u = mgh_stagger(pn, DELTA=[-1,0])

   mask_u = fix(round(mask[0:d0-2,*] * mask[1:d0-1,*]))
   mask_v = fix(round(mask[*,0:d1-2] * mask[*,1:d1-1]))

   ;; Set up result array

   psi = dblarr(dim-[1,1])

   ;; Calculate perimeter values of stream function by integration
   ;; along boundaries (S, E, N, W)...

   ;; ...calculate the stream function increment at each normal
   ;; velocity point

   mskp = round([ mask_v[1:d0-2,0], $
                  reform(mask_u[d0-2,1:d1-2]), $
                  reverse(mask_v[1:d0-2,d1-2]), $
                  reverse(reform(mask_u[0,1:d1-2])) ])

   delp = mskp * $
          double([ + v[1:d0-2,0]/pm_v[1:d0-2,0], $
                   - reform(u[d0-2,1:d1-2]/pn_u[d0-2,1:d1-2]), $
                   - reverse(v[1:d0-2,d1-2]/pm_v[1:d0-2,d1-2]), $
                   + reverse(reform(u[0,1:d1-2]/pn_u[0,1:d1-2])) ])

   ;; ...adjust increment (at sea points only) so that the total is
   ;; zero

   l_mskp = where(mskp, n_mskp)
   case n_mskp gt 0 of
      0: convergence = 0
      1: begin
         convergence = total(delp*mskp, /DOUBLE)
         message, /INFORM, 'Adjusting for boundary convergence of '+ $
                  mgh_format_float(1.E-6*convergence)+' Sv'
         delp[l_mskp] -=  convergence / n_mskp
      end
   endcase

   ;; ...integrate to get stream function on perimeter; value is zero
   ;; at SW corner

   n_p = n_elements(delp)

   psip = dblarr(n_p)
   for i=1,n_p-1 do psip[i] = psip[i-1] + delp[i-1]*mskp[i-1]

   ;; ...load stream function into 2D array

   psi[0:d0-3,0] = psip[0:d0-3]
   psi[d0-2,0:d1-3] = psip[d0-2:d0+d1-5]
   psi[1:d0-2,d1-2] = reverse(psip[d0+d1-4:2*d0+d1-7])
   psi[0,1:d1-2] = reverse(psip[2*d0+d1-6:2*d0+2*d1-9])

   mgh_undefine, delp, mskp, n_mskp, l_mskp, n_p, psip

   ;; Estimate interior values by integrating in xi direction

   psix = psi

   for j=1,d1-3 do begin

      mskx = mask_v[1:d0-2,j]

      delx = double(mask_v[1:d0-2,j]*v[1:d0-2,j]/pm_v[1:d0-2,j])

      l_mskx = where(mskx, n_mskx)
      if n_mskx gt 0 then begin
         adjx = psix[d0-2,j] - psix[0,j] - total(delx)
         delx[l_mskx] +=  adjx / n_mskx
      endif

      for i=1,d0-3 do psix[i,j] = psix[i-1,j] + delx[i-1]

   endfor

   ;; Estimate interior values by integrating in eta direction

   psie = temporary(psi)

   for i=1,d0-3 do begin

      mske = mask_u[i,1:d1-2]

      dele = - double(mask_u[i,1:d1-2]*u[i,1:d1-2]/pn_u[i,1:d1-2])

      l_mske = where(mske, n_mske)
      if n_mske gt 0 then begin
         adje = psie[i,d1-2] - psie[i,0] - total(dele)
         dele[l_mske] +=  adje / n_mske
      endif

      for j=1,d1-3 do psie[i,j] = psie[i,j-1] + dele[j-1]

   endfor

   ;; Result is average of the two estimates

   message, /INFORM, 'Combining 2 stream function fields, maximum discrepancy is '+ $
            strjoin(mgh_format_float(1.E-6*mgh_minmax(psix-psie)), ' ')+' Sv'

   result = 0.5*(temporary(psix)+temporary(psie))

   return, result

end
