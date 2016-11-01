;+
; NAME:
;   MGH_ROMS_PLOT_LAYOUT
;
; PURPOSE:
;   Generate & display a figure showing the staggering of variables on
;   the ROMS grid e.g. see
;
;   http://hadfieldm.greta.niwa.co.nz/boundaries/figures/ROMS_horizontal_grid.png
;
; POSITIONAL PARAMETERS:
;   option (input, integer scalar):
;     Specifies the type of plot required: 0 for horizontal grid; 1 for vertical
;     grid. Default is 0.
;
;###########################################################################
; Copyright (c) 2000-2012 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, ????:
;     Written as MGH_ROMS_GRID_LAYOUT or some such.
;   Mark Hadfield, ????:
;     Renamed MGH_ROMS_GRD_LAYOUT.
;   Mark Hadfield, 2010-10:
;     Renamed MGH_ROMS_PLOT_LAYOUT.
;-
pro mgh_roms_plot_layout, option, LM=lm, MM=mm, N=n

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(option) eq 0 then option = 0

  if n_elements(lm) eq 0 then lm = 5
  if n_elements(mm) eq 0 then mm = 4
  if n_elements(n) eq 0 then n = 6

  case option of

    0: begin

      mgh_graph_default, FONTSIZE=fsize, SCALE=scale

      margin = 0.15

      ograph = obj_new('MGHgrGraph', SCALE=1.3*scale, $
        VIEWPLANE_RECT=[-margin,-margin,1+2*margin,1+2*margin])

      ograph->GetProperty, DELTAZ=deltaz

      ograph->SetProperty, NAME='ROMS horizontal grid'

      ograph->NewFont, SIZE=1.1*fsize
      ograph->NewFont, SIZE=fsize, NAME='Helvetica*Italic'

      ;; Define grids

      x_rho = rebin(mgh_range(0,1,N_ELEMENTS=lm+2),lm+2,mm+2)
      y_rho = rebin(reform(mgh_range(0,1,N_ELEMENTS=mm+2),1,mm+2),lm+2,mm+2)

      x_u = mgh_roms_stagger(x_rho, TO='u')
      y_u = mgh_roms_stagger(y_rho, TO='u')

      x_v = mgh_roms_stagger(x_rho, TO='v')
      y_v = mgh_roms_stagger(y_rho, TO='v')

      x_psi = mgh_roms_stagger(x_rho, TO='psi')
      y_psi = mgh_roms_stagger(y_rho, TO='psi')

      ;; Rho points
      ograph->NewSymbol, 0, /FILL, NORM_SIZE=0.01, $
        COLOR=mgh_color('black'), RESULT=osym
      ograph->NewAtom, 'IDLgrPlot', x_rho[*], y_rho[*], LINESTYLE=6, SYMBOL=osym

      ;; U points
      ograph->NewSymbol, 0, NORM_SIZE=0.01, $
        COLOR=mgh_color('blue'), RESULT=osym
      ograph->NewAtom, 'IDLgrPlot', x_u[*], y_u[*], LINESTYLE=6, SYMBOL=osym

      ;; V points
      ograph->NewSymbol, 0, NORM_SIZE=0.01, $
        COLOR=mgh_color('(0,159,0)'), RESULT=osym
      ograph->NewAtom, 'IDLgrPlot', x_v[*], y_v[*], LINESTYLE=6, SYMBOL=osym

      ;; Draw a line around the physical boundary
      ograph->NewAtom, 'IDLgrPlot', COLOR=mgh_color('red'), $
        [x_psi[*,0],reform(x_psi[lm,*]), $
        reverse(x_psi[*,mm]),reverse(reform(x_psi[0,*]))], $
        [y_psi[*,0],reform(y_psi[lm,*]), $
        reverse(y_psi[*,mm]),reverse(reform(y_psi[0,*]))]

      ;; Labels

      del = [0.015,0.015]
      ograph->NewText, '!9r', LOCATIONS=[x_rho[2,2],y_rho[2,2]]+del, $
        /ENABLE_FORMATTING
      ograph->NewText, '!8u', LOCATIONS=[x_u[2,2],y_u[2,2]]+del, $
        /ENABLE_FORMATTING
      ograph->NewText, '!8v', LOCATIONS=[x_v[2,2],y_v[2,2]]+del, $
        /ENABLE_FORMATTING

      del = [0,-0.07]
      ograph->NewText, 'i=0', LOCATIONS=[x_rho[0,0],y_rho[0,0]]+del, ALIGN=0.5, $
        /ENABLE_FORMATTING, FONT=ograph->GetFont(POSITION=1)
      ograph->NewText, 'i=L', LOCATIONS=[x_rho[lm+1,0],y_rho[lm+1,0]]+del, $
        ALIGN=0.5, /ENABLE_FORMATTING, $
        FONT=ograph->GetFont(POSITION=1)

      del = [-0.03,0]
      ograph->NewText, 'j=0', LOCATIONS=[x_rho[0,0],y_rho[0,0]]+del, ALIGN=1, $
        VERTICAL=0.5, /ENABLE_FORMATTING, $
        FONT=ograph->GetFont(POSITION=1)
      ograph->NewText, 'j=M', LOCATIONS=[x_rho[0,mm+1],y_rho[0,mm+1]]+del, $
        ALIGN=1, VERTICAL=0.5, /ENABLE_FORMATTING, $
        FONT=ograph->GetFont(POSITION=1)

      mgh_new, 'MGH_Window', ograph

    end

    1: begin

      mgh_graph_default, FONTSIZE=fsize, SCALE=scale

      margin = 0.15

      ograph = obj_new('MGHgrGraph', SCALE=1.3*scale, $
        VIEWPLANE_RECT=[-margin,-margin,1+2*margin,1+2*margin])

      ograph->GetProperty, DELTAZ=deltaz

      ograph->SetProperty, NAME='ROMS vertical grid'

      ograph->NewFont, SIZE=fsize
      ograph->NewFont, SIZE=0.9*fsize, NAME='Helvetica*Italic'

      ;; Define grids

      s_w = mgh_range(0, 1, N_ELEMENTS=n+1)

      s_rho = mgh_stagger(s_w, DELTA=[-1])

      xx = [0.25,0.75]

      ;; W levels
      ograph->NewAtom, 'IDLgrPlot', xx, [s_w[0],s_w[0]], THICK=2, $
        COLOR=mgh_color('brown')
      for k=1,n-1 do begin
        ograph->NewAtom, 'IDLgrPlot', xx, [s_w[k],s_w[k]]
      endfor
      ograph->NewAtom, 'IDLgrPlot', xx, [s_w[n],s_w[n]], THICK=2, $
        COLOR=mgh_color('blue')

      ;; Rho levels
      for k=0,n-1 do begin
        ograph->NewAtom, 'IDLgrPlot', xx, [s_rho[k],s_rho[k]], LINESTYLE=2
      endfor

      ;; Labels

      del = 0.015
      ograph->NewText, '!9r', LOCATIONS=[xx[1]+del,s_rho[2]], $
        ALIGNMENT=0, VERTICAL_ALIGN=0.5, /ENABLE_FORMATTING
      ograph->NewText, '!8w', LOCATIONS=[xx[1]+del,s_w[2]], $
        ALIGNMENT=0, VERTICAL_ALIGN=0.5, /ENABLE_FORMATTING

      mgh_new, 'MGH_Window', ograph

    end

  endcase

end
