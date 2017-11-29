;+
; NAME:
;   MGH_SAN_PATH
;
; PURPOSE:
;   Given base names and locations of one or more files on a remote
;   data store (SAN), this routine constructs fully qualified file
;   names. Depending on the setting of the MIRROR keyword, these file names
;   point either to the master data store (MIRROR=0) or to a local mirror (MIRROR=1).
;   Several different master-mirror pairs may be available; they are selected via the
;   volume keyword
;
; CALLING SEQUENCE:
;   result = MGH_SAN_PATH(name, MIRROR=mirror, SUBDIRECTORY=subdirectory)
;
; POSITIONAL PARAMETERS:
;   name (input, optional, string scalar)
;     The file name with no directory portion. Default is ''.
;
; KEYWORD PARAMETERS:
;   MIRROR (input, switch)
;     Set this keyword to return a path on the mirror area on my PC.
;
;   SUBDIRECTORY (input, string array or scalar)
;     Subdirectories relative to the root of the SAN volume.
;
;   VOLUME (input, optional, string scalar)
;     The volume label. Default is '/hpcf/working/hadfield'
;
; RETURN VALUE:
;   The function returns the full path name of the file.
;
; DEPENDENCIES:
;   Highly specific to my setup. Requires a system variable !MGH_SAN
;   to be defined as an array of structures with master-mirror pairs.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-09:
;     Written.
;   Mark Hadfield, 2009-05:
;     In call to strsplit, changed PRESERVE_NULL keyword from 1 to 0
;     to avoid problems with leading or trailing directory separators.
;     I'm not sure of the consequences of this, so let's see.
;   Mark Hadfield, 2015-07:
;     The !MGH_SAN system variable is now assumed to be a hash rather
;     than an array. Volumes are selected via a hash key string rather
;     than an index number.
;   Mark Hadfield, 2015-07:
;     The routine now allows an 'access' field in each structure
;     in !MGH_SAN. This field can have values 'mounted' (the default) or
;     'ssh'. If the latter then the master is assumed to be accessed
;     by ssh and the result when MIRROR is not set is formed accordingly.
;-
function mgh_san_path, name, $
     MIRROR=mirror, SUBDIRECTORY=subdirectory, VOLUME=volume

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(name) eq 0 then name = ''

  n_name = n_elements(name)

  if n_elements(mirror) eq 0 then mirror = 0B

   if n_elements(volume) eq 0 then volume = mgh_san_default()

  if ~ isa(volume, 'STRING') then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'volume'

  vol = (!mgh_san)[volume]
  vol_access = mgh_struct_has_tag(vol, 'access') ? vol.access : 'mounted'

  case n_elements(subdirectory) of
    0:
    1: begin
      my_subdir = strsplit(subdirectory, '[\\/]', /REGEX, PRESERVE_NULL=0, /EXTRACT)
    end
    else: begin
      my_subdir = subdirectory
    end
  endcase

  if keyword_set(mirror) then begin
    result = filepath(name, ROOT=vol.mirror, SUBDIRECTORY=my_subdir)
  endif else begin
    case vol_access of
      'mounted': begin
        result = filepath(name, ROOT=vol.master, SUBDIRECTORY=my_subdir)
      end
      'ssh': begin
        if n_elements(my_subdir) gt 0 then begin
          result = strarr(n_name)
          for n=0,n_name-1 do $
            result[n] = strjoin([vol.master,my_subdir,name[n]], '/')
        endif else begin
          result = strarr(n_name)
          for n=0,n_name-1 do $
            result[n] = strjoin([vol.master,name[n]], '/')
        endelse
      end
    endcase
  endelse

  return, result

end
