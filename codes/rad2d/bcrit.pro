function get_beta, bcool  

result = newton(bcool,'newtfunc')

return, result(0) 
end 

function newtfunc, bcool 

common share, qsmall, theta, bulkv, alphafix, viscg, Qval, gmma, sstar 

bigQ_marg = 1d0/sqrt(2d0*gmma*(2d0-qsmall))

if(alphafix gt 0d0) then begin 
alpha = alphafix 
endif else begin 
alpha = 1d0 - theta 
alpha/= (gmma-1d0)*bcool*qsmall^2
endelse

if(viscg gt 0d0) then begin
bigQ = Qval*bigQ_marg/sqrt(alpha) 
endif else begin
bigQ = Qval*bigQ_marg 
endelse 

res = 6.^(3./4.)*(alpha/(2d0-qsmall))^(1./4.)*(4.*alpha*(1d0+3d0*bulkv/4d0) + 3.*gmma*bcool)^(1./4.) - 3.*theta
res/= bigQ*(4.*alpha*(1d0+3d0*bulkv/4d0) + 3.*gmma*bcool) 
res-= sstar 

return, res 
end


pro bcrit, tgrow=tgrow, smallq=smallq, fixalpha=fixalpha, gvisc=gvisc, bigQ=bigQ, tirr=tirr, bvisc = bvisc $
          , label=label, legend=legend, title=title, xrange=xrange, yrange=yrange  

!p.font = 0 

common share, qsmall, theta, bulkv, alphafix, viscg, Qval, gmma, sstar 

if not keyword_set(smallq) then smallq = 1.5
if not keyword_set(tirr) then tirr = 0d0 
if not keyword_set(bvisc) then bvisc = 0d0 
if not keyword_set(fixalpha) then fixalpha = 0d0 

qsmall = smallq 
theta  = tirr 
bulkv = bvisc 
alphafix = fixalpha 
if keyword_set(gvisc) then begin
viscg = 1d0
endif else viscg = -1d0 
if not keyword_set(bigQ) then bigQ = 1d0 
Qval = bigQ 

gmma_min = 1.01 
gmma_max = 3.0
ng       = 128  
gmma_arr = gmma_min + dindgen(ng)*(gmma_max-gmma_min)/(ng-1d0)

ncases   = n_elements(tgrow) 
bvals    = dblarr(ncases, ng) 
rates    = 1d0/(2d0*!dpi*tgrow) 

;bval_theory = 1d0/(sqrt(gmma_arr)-1d0)^1.5 

bvals_theory = dblarr(ncases,ng) 


bvals(*,*) = 1d0 

for i=0, ncases-1 do begin
sstar = rates(i) 
for j=0, ng-1 do begin
gmma = gmma_arr(j) 
bvals(i,j) = get_beta(bvals(i,j))


bigQ_marg = 1d0/sqrt(2d0*gmma*(2d0-qsmall))
bcool = bvals(i,j) 
if(alphafix gt 0d0) then begin
alpha = alphafix
endif else begin
alpha = 1d0 - theta
alpha/= (gmma-1d0)*bcool*qsmall^2
endelse

if(viscg gt 0d0) then begin
bigQ = Qval*bigQ_marg/sqrt(alpha)
endif else begin
bigQ = Qval*bigQ_marg
endelse

bigX = sstar^2 + 2d0*(2d0-smallq) 
bigX*= bigQ^2 

bvals_theory(i,j) = (bigX*tirr - 1d0)/(sstar*(1d0-bigX*gmma)) 

endfor
endfor


     loadct, 6,/silent
     color_arr = dindgen(ncases)*256d0/ncases


     if not keyword_set(yrange) then begin
        yrange1=[min(bvals),max(bvals)]
     endif else begin
        yrange1=yrange
     endelse

     xtitle = textoidl('\gamma') + '!X'
     ytitle = textoidl('\beta') + '!X'
     if not keyword_set(title) then title = ''

     set_plot,'ps'
     file = strcompress('bcrit.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
     plot, gmma_arr, bvals(0,*),xmargin=[8.,2],ymargin=[3.2,1.8], ystyle=1, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange1,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, color=color_arr(0), /ylog  $
     ,title = textoidl(title)

     oplot, gmma_arr, bvals_theory(0,*), thick=4, linestyle=1 , color=color_arr(0) 
     for i=1,ncases-1 do begin
           oplot, gmma_arr, bvals(i,*), thick=4,color=color_arr(i)
           oplot, gmma_arr, bvals_theory(i,*), thick=4, linestyle=1 , color=color_arr(i) 
     endfor

     if keyword_set(legend) then begin
        x0=legend(0)
        x1=legend(1)
        y0=legend(2)
        dy=legend(3)
        for j=0, n_elements(label)-1 do begin
           ;ynew = 10.^(alog10(y0) - dy*j)
           ynew = y0 - dy*j 
           oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j)
           xyouts, x1, ynew, textoidl(label(j)),charsize=1.5
        endfor
     endif
   
     device,/close 


end
