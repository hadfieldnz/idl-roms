;+
; NAME:
;   MGH_SAN_FILE
;
; PURPOSE:
;   Given base names and locations of one or more files on a remote
;   data store (SAN), this routine checks that each file exists and
;   is readable. (If any of the files does not exist or is not
;   readable, the function raises an error.) Depending on the setting
;   of the MIRROR keyword, it will either return the fully qualified
;   path names of the master files, or make local copies (as
;   necessary) and return the names of those.
;
; CALLING SEQUENCE:
;   result = mgh_san_file(name)
;
; POSITIONAL ARGUMENTS:
;   name (input, string scalar or vector)
;     The base name, i.e. with no directory portion, of the file(s).
;
; KEYWORD ARGUMENTS:
;   MIRROR (input, switch)
;     Specifies where to look for the file and what action to
;     take. Valid values are:
;       0   Return path name of files on remote data store.
;       1   Make local copies (if necessary) and return the names of
;           those.
;     Default is 1 if a mirror is specified for that volume, otherwise 0.
;
;   MTIME (output, long or long64 scalar)
;     If this argument is present, return the latest last-modified
;     time stamp for the files in the sequence.
;
;   SUBDIRECTORY (input, string)
;     Subdirectories relative to the SAN volume.
;
;   VOLUME (input, optional, string scalar)
;     The volume name. The default is returned by MGH_SAN_DEFAULT.
;
; RETURN VALUE:
;   The function returns the full path name of the file(s).
;
; DEPENDENCIES:
;   Requires MGH_SAN_PATH to form path names. This requires a system variable,
;   !MGH_SAN, to be defined, as described in the documentation.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-09:
;     Written.
;   Mark Hadfield, 2011-04:
;     Now prints elapsed time if a copy is made.
;   Mark Hadfield, 2011-08:
;     - Added GUNZIP and MTIME keywords.
;   Mark Hadfield, 2011-08:
;     - GUNZIP keyword renamed GZIP and option 2 added, provided
;       automagical handling of remote gzipped files.
;   Mark Hadfield, 2012-01:
;     - Automagical handling of bzipped files on the master
;       was added.
;     - GZIP keyword deleted: the routine always searches for
;       bzipped, gzipped and non-packed files.
;   Mark Hadfield, 2013-02:
;     - Added a check for zero-sized mirror files, assumed to be
;       the result of an unsuccessful previous attempt to retrieve
;       the file.
;   Mark Hadfield, 2013-06:
;     - Copying of non-compressed files now being done by
;       MGH_FILE_COPY rather than FILE_COPY, in an attempt to reduce
;       the impact on other applications accessing the same network
;       volume.
;   Mark Hadfield, 2014-08:
;     - Changed preference order for packed vs unpacked files. An
;       unpacked is now preferred, if it exists, followed by packed
;       files with suffixes .gz and .bz2
;   Mark Hadfield, 2015-07:
;     The !MGH_SAN system variable is now assumed to be a hash rather
;     than an array. Volumes are selected via a hash key string rather
;     than an index number.
;   Mark Hadfield, 2015-07:
;     Now supports access to remote volumes via mounted file systems
;     (internal routine MGH_SAN_FILE_MOUNTED) or scp (internal routine MGH_SAN_FILE_SSH).
;   Mark Hadfield, 2016-01
;     Re-indented source code.
;   Mark Hadfield, 2019-05
;     In checking for the existence of the master file, this routine used
;     to call file_test with the READ keyword set. However this test gives
;     incorrect results on some Maui volumes (the file exists and can be 
;     read or copied, but the test returns !false). So now we just call
;     file_test.
;-
function mgh_san_file_mounted, name, $
     MIRROR=mirror, MTIME=mtime, SUBDIRECTORY=subdirectory, $
     VOLUME=volume

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   compile_opt HIDDEN

   if n_elements(volume) eq 0 then volume = mgh_san_default()

   if n_elements(mirror) eq 0 then begin
      msv = (!mgh_san)[volume]
      mirror = mgh_struct_has_tag(msv, 'mirror') && strlen(msv.mirror) gt 0
   endif

   n_file = n_elements(name)

   ;; Determine fully qualified names of master and mirror files.
   ;; Use of the MGH_REPRODUCE function in setting up the variables
   ;; ensures they match the input in size and shape.

   file_master = mgh_reproduce('', name)
   file_mirror = mgh_reproduce('', name)
   file_packed = mgh_reproduce(0B, name)

   for f=0,n_file-1 do begin
      file_mirror[f] = mgh_san_path(name[f], MIRROR=mirror, SUBDIRECTORY=subdirectory, VOLUME=volume)
      file_master[f] = mgh_san_path(name[f], MIRROR=0, SUBDIRECTORY=subdirectory, VOLUME=volume)
      case !true of
         file_test(file_master[f]): begin
            file_packed[f] = 0B
         end
         file_test(file_master[f]+'.gz'): begin
            file_master[f] += '.gz'
            file_packed[f] = 1B
         end
         file_test(file_master[f]+'.bz2'): begin
            file_master[f] += '.bz2'
            file_packed[f] = 2B
         end
         else: $
            message, 'File does not exist or is unreadable: '+file_master[f]
      endcase
   endfor

   ;; Copy files if necessary

   for f=0,n_file-1 do begin

      if (file_packed[f] gt 0 || mirror) then begin

         info_master = file_info(file_master[f])
         info_mirror = file_info(file_mirror[f])
         copy = (~ info_mirror.read) || (info_mirror.mtime lt info_master.mtime) || (info_mirror.size eq 0)
         if copy then begin
            message, /INFORM, string(FORMAT='(%"Copying %s to %s")', file_master[f], file_mirror[f])
            file_mkdir, file_dirname(file_mirror[f])
            if file_packed[f] gt 0 then begin
               t0 = systime(1)
               mgh_file_copy, file_master[f], file_mirror[f], $
                  BUNZIP2=file_packed[f] eq 2, GUNZIP=file_packed[f] eq 1
               t1 = systime(1)
            endif else begin
               t0 = systime(1)
               mgh_file_copy, file_master[f], file_mirror[f]
               t1 = systime(1)
            endelse
            size_mb = info_master.size/1048576.0
            if file_packed[f] gt 0 then begin
               info_mirror = file_info(file_mirror[f])
               size_mbu = info_mirror.size/1048576.0
               fmt = '(%"Copied %0.1f(%0.1f) MiB in %0.2f s, %0.2f(%0.2f) MiB/s")'
               message, /INFORM, $
                  string(FORMAT=fmt, size_mb, size_mbu, t1-t0, size_mb/(t1-t0), size_mbu/(t1-t0))
            endif else begin
               fmt = '(%"Copied %0.1f MiB in %0.2f s, %0.2f MiB/s")'
               message, /INFORM, $
                  string(FORMAT=fmt, size_mb, t1-t0, size_mb/(t1-t0))
            endelse
         endif
      endif

   endfor

   if arg_present(mtime) then begin
      mtime = 0
      for f=0,n_file-1 do begin
         info = file_info(file_mirror[f])
         mtime = mtime > info.mtime
      endfor
   endif

   return, file_mirror

end

function mgh_san_file_ssh, name, $
     MIRROR=mirror, MTIME=mtime, SUBDIRECTORY=subdirectory, $
     VOLUME=volume

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   compile_opt HIDDEN

   if n_elements(volume) eq 0 then volume = mgh_san_default()

   if n_elements(mirror) eq 0 then mirror = 1B

   n_file = n_elements(name)

   ;; Determine fully qualified names of master and mirror files.
   ;; Use of the MGH_REPRODUCE function in setting up the variables
   ;; ensures they match the input in size and shape.

   file_master = mgh_reproduce('', name)
   file_mirror = mgh_reproduce('', name)
   file_packed = mgh_reproduce(0B, name)

   for f=0,n_file-1 do begin
      file_mirror[f] = mgh_san_path(name[f], MIRROR=mirror, SUBDIRECTORY=subdirectory, VOLUME=volume)
      file_master[f] = mgh_san_path(name[f], MIRROR=0, SUBDIRECTORY=subdirectory, VOLUME=volume)
      case 1B of
         mgh_ssh_file_test(file_master[f]): begin
            file_packed[f] = 0B
         end
         mgh_ssh_file_test(file_master[f]+'.gz'): begin
            file_master[f] += '.gz'
            file_packed[f] = 1B
         end
         mgh_ssh_file_test(file_master[f]+'.bz2'): begin
            file_master[f] += '.bz2'
            file_packed[f] = 2B
         end
         else: $
            message, 'File does not exist or is unreadable: '+file_master[f]
      endcase
   endfor

   ;; Copy files if necessary

   for f=0,n_file-1 do begin

      if (file_packed[f] gt 0 || mirror) then begin

         info_master = mgh_ssh_file_info(file_master[f])
         info_mirror = file_info(file_mirror[f])
         copy = (~ info_mirror.read) || (info_mirror.mtime lt info_master.mtime) || (info_mirror.size eq 0)
         if copy then begin
            message, /INFORM, string(FORMAT='(%"Copying %s to %s")', file_master[f], file_mirror[f])
            file_mkdir, file_dirname(file_mirror[f])
            if file_packed[f] gt 0 then begin
               t0 = systime(1)
               mgh_ssh_file_copy, file_master[f], file_mirror[f], $
                  BUNZIP2=file_packed[f] eq 2, GUNZIP=file_packed[f] eq 1
               t1 = systime(1)
            endif else begin
               t0 = systime(1)
               mgh_ssh_file_copy, file_master[f], file_mirror[f]
               t1 = systime(1)
            endelse
            size_mb = info_master.size/1048576.0
            if file_packed[f] gt 0 then begin
               info_mirror = file_info(file_mirror[f])
               size_mbu = info_mirror.size/1048576.0
               fmt = '(%"Copied %0.1f(%0.1f) MiB in %0.2f s, %0.2f(%0.2f) MiB/s")'
               message, /INFORM, $
                  string(FORMAT=fmt, size_mb, size_mbu, t1-t0, size_mb/(t1-t0), size_mbu/(t1-t0))
            endif else begin
               fmt = '(%"Copied %0.1f MiB in %0.2f s, %0.2f MiB/s")'
               message, /INFORM, $
                  string(FORMAT=fmt, size_mb, t1-t0, size_mb/(t1-t0))
            endelse
         endif
      endif

   endfor

   if arg_present(mtime) then begin
      mtime = 0
      for f=0,n_file-1 do begin
         info = file_info(file_mirror[f])
         mtime = mtime > info.mtime
      endfor
   endif

   return, file_mirror

end

function mgh_san_file, name, VOLUME=volume, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(name) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'name'

   if n_elements(volume) eq 0 then volume = mgh_san_default()

   vol = (!mgh_san)[volume]
   vol_access = mgh_struct_has_tag(vol, 'access') ? vol.access : 'mounted'

   case vol_access of
      'mounted': $
         return, mgh_san_file_mounted(name, VOLUME=volume, _STRICT_EXTRA=extra)
      'ssh': $
         return, mgh_san_file_ssh(name, VOLUME=volume, _STRICT_EXTRA=extra)
   endcase

end
