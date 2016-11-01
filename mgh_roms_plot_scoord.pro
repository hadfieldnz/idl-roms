;+
; NAME:
;   MGH_ROMS_PLOT_SCOORD
;
; PURPOSE:
;   Plot or print vertical coordinate data
;
; POSITIONAL PARAMETERS:
;   option (input, integer scalar):
;     Specifies the type of plot/printout required.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2010-10:
;     Written.
;   Mark Hadfield, 2011-05:
;     Changed default parameters: now uses the Shchepetkin double stretching
;     function (VSTRETCH=4).
;   Mark Hadfield, 2012-10:
;     Option 0 now uses Motley graphics instead of IDL 8 new graphics.
;-
pro mgh_roms_plot_scoord, option, $
     H_MIN=h_min, H_MAX=h_max, HC=hc, N=n, THETA_S=theta_s, THETA_B=theta_b, $
     VSTRETCH=vstretch, VTRANSFORM=vtransform

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(option) eq 0 then option = 0

   if n_elements(h_max) eq 0 then h_max = 5000
   if n_elements(h_min) eq 0 then h_min = 0.02*h_max

   if n_elements(n) eq 0 then n = 30
   if n_elements(theta_s) eq 0 then theta_s = 4.0
   if n_elements(theta_b) eq 0 then theta_b = 2.0
   if n_elements(hc) eq 0 then hc = 50.0
   if n_elements(vstretch) eq 0 then vstretch = 4
   if n_elements(vtransform) eq 0 then vtransform = 2

   s_w = mgh_range(-1, 0, STRIdE=1./n)

   cs_w = mgh_roms_s_to_cs(s_w, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

   case option of

      0: begin

         ang = mgh_range(0, !pi/2, STRIDE=!pi/100)
         h = h_min + (h_max-h_min)*(sin(ang)^2)

         z_w = fltarr(n_elements(s_w), n_elements(h))

         for j=0,n_elements(h)-1 do begin
            z_w[*,j] = mgh_roms_s_to_z(s_w, h[j], $
                                       CS=cs_w, HC=hc, VTRANSFORM=vtransform)
         endfor

         ;; Plot with Motley graphics

         ograph = obj_new('MGHgrGraph2D')

         ograph->NewAxis, 0, RANGE=[0,n_elements(h)-1]
         ograph->NewAxis, 1, RANGE=[-h_max,0], TICKFORMAT='mgh_tf_negative'

         ograph->NewAtom, 'IDLgrPlot', -h
         for k=0,n_elements(s_w)-1 do $
              ograph->NewAtom, 'IDLgrPlot', z_w[k,*], COLOR=!color.blue

         mgh_new, 'MGH_Window', ograph

      end

      1: begin

         z_w_min = mgh_roms_s_to_z(s_w, h_min, $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         z_w_mid = mgh_roms_s_to_z(s_w, 0.5*(h_min+h_max), $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         z_w_max = mgh_roms_s_to_z(s_w, h_max, $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         print, 'Vertical S-coordinate System:'
         print, 'Level   S-coord     Cs-curve   at_hmin  at_hmid  at_hmax'
         fmt = '(I5,2F10.6,3F10.3)'
         for k=n_elements(s_w)-1,0,-1 do $
               print, FORMAT=fmt, k, s_w[k], cs_w[k], z_w_min[k], z_w_mid[k], z_w_max[k]

;         !null = plot(s_w)
;         !null = plot(cs_w)

      end

      2: begin

         ograph = obj_new('MGHgrGraph2D', ASPECT=0.7, NAME='S coordinate stretching function')

         ograph->NewFont

         ograph->NewAxis, 0, RANgE=[-1,0], TITLE='S'
         ograph->NewAxis, 1, RANgE=[-1,0], TITLE='Cs'

         ograph->NewAtom, 'IDLgrPlot', s_w, cs_w, COLOR=!color.black

         mgh_new, 'mgh_window', ograph

      end

      3: begin

         z_w_min = mgh_roms_s_to_z(s_w, h_min, $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         z_w_mid = mgh_roms_s_to_z(s_w, 0.5*(h_min+h_max), $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         z_w_max = mgh_roms_s_to_z(s_w, h_max, $
                                   CS=cs_w, HC=hc, VTRANSFORM=vtransform)

         ograph = obj_new('MGHgrGraph2D', ASPECT=0.7, NAME='S coordinate transformation')

         ograph->NewFont

         ograph->NewAxis, 0, RANgE=[-1,0], TITLE='Cs'
         ograph->NewAxis, 1, RANgE=[-h_max,0], /EXACT, /EXTEND, TITLE='Depth (m)', $
              TICKFORMAT='mgh_tf_negative'

         ograph->NewAtom, 'IDLgrPlot', cs_w, z_w_min, COLOR=!color.black
         ograph->NewAtom, 'IDLgrPlot', cs_w, z_w_mid, COLOR=!color.black
         ograph->NewAtom, 'IDLgrPlot', cs_w, z_w_max, COLOR=!color.black

         mgh_new, 'mgh_window', ograph

      end

   endcase

end
