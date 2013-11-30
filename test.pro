pro test

  a=readfits('CLEAN_INTEN_GONG_X_T.fits',hdr)
  c=readfits('data/tdizi100212t1033.fits',hdr2)
 
  ss = size(a)

  b=fltarr(ss[1],ss[2],ss[3]-1)

  for i=0,ss[3]-2 do begin
     b[*,*,i]=a[*,*,i+1]-a[*,*,i]
     writefits,'RAW_DATA/raw_data'+STRTRIM(i,2)+'.fits',b[*,*,i],hdr2
  endfor


end
