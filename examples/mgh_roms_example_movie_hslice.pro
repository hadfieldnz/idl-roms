;+
; NAME:
;   MGH_ROMS_EXAMPLE_MOVIE_HSLICE
;
; PURPOSE:
;   Hslice movie examples
;
;###########################################################################
; Copyright (c) 2016 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-11:
;     Written.
;-
pro mgh_roms_example_movie_hslice, option

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(option) eq 0 then option = 0

   root = file_dirname(routine_filepath('mgh_roms_example_movie_hslice'))

   case option of

      0: begin

         name = 'Lake Jersey SSH'

         file = filepath('lake_jersey_his.nc', ROOT=root, SUBDIR=['data','lake_jersey','Forward'])

         mgh_new, 'mgh_roms_movie_hslice', file, 'zeta', $
            DATA_RANGE=[-0.04,0.04], $
            GRAPH_PROPERTIES={name: name}

      end

      1: begin

         name = 'Lake Jersey surface speed'

         file = filepath('lake_jersey_his.nc', ROOT=root, SUBDIR=['data','lake_jersey','Forward'])

         mgh_new, 'mgh_roms_movie_hslice', file, 'abs(u,v)', $
            DATA_RANGE=[0,0.5], $
            GRAPH_PROPERTIES={name: name}

      end

   endcase


end
