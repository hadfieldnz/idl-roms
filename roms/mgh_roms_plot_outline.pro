;+
; PROCEDURE NAME:
;   MGH_ROMS_PLOT_OUTLINE
;
; PURPOSE:
;   For one or more ROMS grid files, plot the (normal velocity) boundary
;   on top of a map with bathymetry
;
; CALLING SEQUENCE
;   mgh_roms_plot_outline, ogrd
;
; POSITIONAL PARAMETERS:
;   ogrd (input, object scalar or vector)
;     One or more references to ROMS grid files
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-05:
;     Written, generalising mgh_hikurangi_fig_outline and friends.
;   Mark Hadfield, 2009-08:
;     Now passes out a reference to the window object via the RESULT
;     keyword.
;   Mark Hadfield, 2009-10:
;     Perimeter (normal-velocity boundary) of each grid now calculated
;     by MGH_PERIM.
;   Mark Hadfield, 2014-10:
;     - Deleted GERRIS_FILE functionality.
;     - Added SHOW_INTERIOR functionality.
;-
pro mgh_roms_plot_outline, ogrd, RESULT=result, $
     DATA_RANGE=data_range, MARGIN=margin, $
     SHOW_INTERIOR=show_interior, INTERIOR_STRIDE=interior_stride, $
     BATHYMETRY_PROPERTIES=bathymetry_properties, $
     PALETTE_PROPERTIES=palette_properties, $
     INTERIOR_PROPERTIES=interior_properties, $
     MAP_PROPERTIES=map_properties, $
     OUTLINE_PROPERTIES=outline_properties, $
     _REF_EXTRA=extra

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(ogrd) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ogrd'

  if n_elements(margin) eq 0 then margin = {lon: [-1.5,1.5], lat: [-1,1]}

  if n_elements(interior_stride) eq 0 then interior_stride = 1

  n_grd = n_elements(ogrd)

  perim = ptrarr(n_grd)

  if keyword_set(show_interior) then grd_psi = ptrarr(n_grd)

  for s=0,n_grd-1 do begin

    grd = ogrd[s]->Retrieve(['lon_rho','lat_rho'])

    ;; Determine points on perimeter (normal-velocity boundary)

    lon_psi = mgh_stagger(grd.lon_rho, DELTA=[-1,-1])
    lat_psi = mgh_stagger(grd.lat_rho, DELTA=[-1,-1])

    if keyword_set(show_interior) then $
      grd_psi[s] = ptr_new({lon: lon_psi, lat: lat_psi})

    perim[s] = ptr_new(mgh_perim(lon_psi, lat_psi))

    mgh_undefine, lon_psi, lat_psi

  endfor

  ;; Determine boundaries for background map

  lon_range = mgh_minmax((*perim[0])[0,*])
  lat_range = mgh_minmax((*perim[0])[1,*])

  for s=1,n_grd-1 do begin
    lon_range = [lon_range[0] < min((*perim[s])[0,*]), $
      lon_range[1] > max((*perim[s])[0,*])]
    lat_range = [lat_range[0] < min((*perim[s])[1,*]), $
      lat_range[1] > max((*perim[s])[1,*])]
  endfor

  bound = {lon: lon_range+margin.lon, lat: lat_range+margin.lat}

  ;; Plot it

  bp = mgh_struct_merge({dataset: 'GEBCO', stride: 1, halo: 2}, $
    bathymetry_properties)

  mp = mgh_struct_merge({name: 'ROMS grid outline', fill: 1, resolution: 1}, $
    map_properties)

  pp = mgh_struct_merge({table: 'Blue/White'}, $
    palette_properties)

  nzregion_plot_bathymetry_density, RESULT=result, $
    BOUNDARIES=bound, DATA_RANGE=data_range, $
    BATHYMETRY_PROPERTIES=bp, $
    MAP_PROPERTIeS=mp, $
    PALETTE_PROPERTIES=pp, $
    _STRICT_EXTRA=extra

  result->GetProperty, GRAPHICS_TREE=ograph

  ograph->GetProperty, DELTAZ=deltaz

  for s=0,n_grd-1 do begin

    ograph->NewAtom, 'IDLgrPlot', (*perim[s])[0,*], (*perim[s])[1,*], $
      ZVALUE=10*deltaz, /USE_ZVALUE, COLOR=!color.red, THICK=2, $
      _STRIcT_EXTRA=outline_properties

    if keyword_set(show_interior) then begin
      gg = *grd_psi[s]
      dim = size(gg.lon, /DIMENSIONS)
      for i=0,dim[0]-1,interior_stride do begin
        ograph->NewAtom, 'IDLgrPlot', gg.lon[i,*], gg.lat[i,*], $
          ZVALUE=10*deltaz, /USE_ZVALUE, COLOR=!color.red, THICK=1, $
          _STRIcT_EXTRA=interior_properties
      endfor
      for j=0,dim[1]-1,interior_stride do begin
        ograph->NewAtom, 'IDLgrPlot', gg.lon[*,j], gg.lat[*,j], $
          ZVALUE=10*deltaz, /USE_ZVALUE, COLOR=!color.red, THICK=1, $
          _STRIcT_EXTRA=interior_properties
      endfor
    endif

  endfor

  result->Update

  heap_free, perim
  if keyword_set(show_interior) then heap_free, grd_psi

end
