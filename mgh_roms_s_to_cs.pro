;+
; NAME:
;   MGH_ROMS_S_TO_CS
;
; PURPOSE:
;   This function calculates the s-coordinate stretching coefficient for a given set
;   of s-coordinate values.
;
; CATEGORY:
;   Ocean models.
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_S_TO_CS(S)
;
; POSITIONAL PARAMETERS:
;   S (input, numeric scalar or array)
;     S-coordinate values, 0 <= S <= 1.
;
; KEYWORD PARAMETERS:
;   THETA_S (input, numeric scalar)
;     Surface control parameter. Default is 0.
;
;   THETA_B (input, numeric scalar)
;      Bottom control parameter. Default is 0.
;
; EXPLANATION:
;   The s coordinate is defined in Song & Haidvogel, 1994, J Comp Phys
;   115, 228-244 and in the ROMS Users Manual.
;
; SEE ALSO:
;   MGH_ROMS_TO_Z
;
;###########################################################################
; Copyright (c) 1998-2012 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, ????:
;     Written.
;   Mark Hadfield, 1998-09:
;     Renamed MGH_ROMS_S_TO_CS.
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options.
;   Mark Hadfield, 2010-10:
;     Added a check for s-coordinate values outside the range [-1,0].
;   Mark Hadfield, 2011-05:
;     Updated for new stretching option 4: the Shchepetkin (2010) double
;     stretching function.
;-
function mgh_roms_s_to_cs, s, $
     THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(s) eq 0 then $
        message, 'One or more S values are required.'

   if min(s) lt -1 || max(s) gt 0 then $
        message, 'Invalid s-coordinate value'

   if n_elements(theta_s) eq 0 then theta_s = 0

   if n_elements(theta_b) eq 0 then theta_b = 0

   if n_elements(vstretch) eq 0 then vstretch = 1

   case vstretch of

      1: begin
         if theta_s gt 0 then begin
            pth = sinh(theta_s*s)/sinh(theta_s)
            rth = tanh(theta_s*(s+0.5D0))/(2*tanh(0.5D0*theta_s)) - 0.5D0
            cs = (1D0-theta_b)*pth + theta_b*rth
         endif else begin
            cs = s
         endelse
      end

      2: begin
         if theta_s gt 0 then begin
            alfa = 1.0D
            beta = 1.0D
            csur = (1-cosh(theta_s*s))/(cosh(theta_s)-1)
            if theta_b gt 0 then begin
               cbot = -1+sinh(theta_b*(s+1))/sinh(theta_b)
               w = (s+1)^alfa*(1+(alfa/beta)*(1-(s+1)^beta))
               cs = w*csur + (1-w)*cbot
            endif else begin
               cs = csur
            endelse
         endif else begin
            cs = s
         endelse
      end

      3: begin
         if theta_s gt 0 then begin
            exp_s = theta_s         ;;; surface stretching exponent
            exp_b = theta_b         ;;; bottom  stretching exponent
            alpha = 3.0D            ;;; scale factor for all hyperbolic functions
            cbot = alog(cosh(alpha*(s+1)^exp_b))/alog(cosh(alpha)) - 1
            csur = -alog(cosh(alpha*abs(s)^exp_s))/alog(cosh(alpha))
            w = (1-tanh(alpha*(s+0.5D)))/2.0D
            cs = w*Cbot + (1-w)*csur
         endif else begin
            cs = s
         endelse
      end

      4: begin
         if theta_s gt 0 then begin
            csur = (1-cosh(theta_s*s))/(cosh(theta_s)-1)
         endif else begin
            csur = -s^2
         endelse
         if theta_b gt 0 then begin
            cs = (exp(theta_b*csur)-1)/(1-exp(-theta_b))
         endif else begin
            cs = csur
         endelse
      end

   endcase

   return, cs

end
