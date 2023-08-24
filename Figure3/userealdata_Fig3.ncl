;*************************************************
; NCL Graphics: lsm_3.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


;;;; how to run this script
;;;; ncl 'period="end_century"' 'season="meanAnnual"' 'format="x11"' Fig3.ncl
;;;; ncl 'period="end_century"' 'season="DJF"' 'format="x11"' Fig3.ncl
;;;; ncl 'period="end_century"' 'season="JJA"' 'format="x11"' Fig3.ncl


undef ("read_sigtest")
function read_sigtest(patho,varname,rcp,period,fname_string)
begin
	data_filename = patho+varname+"_"+rcp+"_"+period+"_siglvl"+fname_string+".nc"
	a_sigtest = addfile(data_filename,"r")
	sigtest=(a_sigtest->sigtest)
	
	return sigtest
end

undef ("cal_diff")
function cal_diff(var_p3,var_p0)
begin
	vardiff = var_p3-var_p0
	copy_VarCoords(var_p3, vardiff)
	vardiff@units = "gC m~S~-2~N~d~S~-1"
	vardiff@cell_methods = "time: mean"
	vardiff@_FillValue = var_p3@_FillValue
	vardiff@missing_value = var_p3@missing_value
	vardiff = lonFlip(vardiff)  ;;; convert 0~360 longitude to -180~180
	return vardiff
end

undef ("cal_diff_lat")
function cal_diff_lat(var_p3,var_p0,unit)
begin
	vardiff = var_p3-var_p0
	copy_VarCoords(var_p3, vardiff)
	vardiff@long_name = unit
	vardiff@units = "gC m~S~-2~N~d~S~-1"
	vardiff@cell_methods = "time: mean"
	vardiff@_FillValue = var_p3@_FillValue
	vardiff@missing_value = var_p3@missing_value
	vardiff = lonFlip(vardiff)  ;;; convert 0~360 longitude to -180~180
	return vardiff
end

undef ("set_boundary_bio_ref")
function set_boundary_bio_ref(varname)
begin	
	if (varname .eq. "TSA") then
		upper_boundary=4.0
		lower_boundary=(-1)*upper_boundary
		space=upper_boundary/10
		LeftString="T-2m"
		unit=LeftString+" (~F34~0~F~C)"
		colors=(/10,15,20,30,40,50,60,70,80,90,100,100,110,120,130,140,150,160,170,180,190,200/) ; from blue to red
	end if
	
	if (varname .eq. "TREFMXAV" ) then
		upper_boundary=4.0
		lower_boundary=(-1)*upper_boundary
		space=upper_boundary/10
		LeftString="Max T-2m"
		unit=LeftString+" (~F34~0~F~C)"
		colors=(/10,15,20,30,40,50,60,70,80,90,100,100,110,120,130,140,150,160,170,180,190,200/) ; from blue to red
	end if
	
	if (varname .eq. "TREFMNAV") then
		upper_boundary=4.0
		lower_boundary=(-1)*upper_boundary
		space=upper_boundary/10
		LeftString="Min T-2m"
		unit=LeftString+" (~F34~0~F~C)"
		colors=(/10,15,20,30,40,50,60,70,80,90,100,100,110,120,130,140,150,160,170,180,190,200/) ; from blue to red
	end if
		
	return [/upper_boundary,lower_boundary,space,unit,colors/]
end


undef ("read_lat_mean")
function read_lat_mean(varname,lats,lons,rcp_bases,rcps,period,fname_string,patho,unit)
begin
	forests = new((/3,dimsizes(lats),dimsizes(lons)/),float)
	bioenergys = new((/3,dimsizes(lats),dimsizes(lons)/),float)
	
	do ensi=0,2
		fname=patho+varname+"_"+rcp_bases(ensi)+"_"+period+fname_string+".nc"
		print((/fname/))
		addfilercp = addfile(fname,"r")
		reforest = (addfilercp->diff)
		printVarSummary(reforest)
		forests(ensi,:,:) = reforest
	
		fname=patho+varname+"_"+rcps(ensi)+"_"+period+fname_string+".nc"
		print((/fname/))
		addfilercp = addfile(fname,"r")
		bioenergy = (addfilercp->diff)
		bioenergys(ensi,:,:) = bioenergy
	
		delete(bioenergy)
		delete(reforest)
	end do
	
	forest =  dim_avg_n_Wrap(forests,0)
	value =  dim_avg_n_Wrap(bioenergys,0)
	
	meanann_value=cal_diff_lat(value,forest,unit)
	plotvalue1 = dim_avg_n_Wrap(meanann_value,1)
	
	return plotvalue1
end

;************************************************
;;;;****************** main function
;************************************************
begin
	wks2 = gsn_open_wks(format,"Fig3"+"_"+season)
	
	gsn_define_colormap(wks2,"testcmap")
	plot2 = new(8,graphic)

	res                      = True                ; plot mods desired
	res@cnFillMode            = "RasterFill"       ; raster plot
	res@cnFillOn              = True               ; color on
	res@lbLabelBarOn          = False              ; no label bar
	res@cnLinesOn             = False              ; no contour lines
	res@cnLineLabelsOn        = False             ; no line labels
	res@cnLevelSelectionMode = "ManualLevels"
	res@gsnDraw               = False             ; don't draw
	res@gsnFrame              = False             ; don't advance frame
	res@gsnAddCyclic          = False
	res@mpFillOn                = True               ; Turn on map fill
	res@mpProjection="Robinson" ;"Mercator" ;LambertEqualArea, Stereographic
	res@mpFillAreaSpecifiers    = (/"Land", "Antarctica"/)
	res@mpSpecifiedFillColors   = (/"white","LightGrey"/)
	res@mpOutlineOn            = True
	res@mpPerimOn              = False
	res@mpGridAndLimbOn        = False
	res@mpGridMaskMode         = "MaskLand"
	res@mpOceanFillColor       = "grey89" ;"LightGrey"	
	res@tmYRBorderOn = False
	res@tmXTBorderOn = False
	res@tmYLBorderOn = False
	res@tmXBBorderOn = False
	res@tmYROn       = False              ; no y right tickmarks
	res@tmXTOn       = False              ; no x top tickmarks
	res@tmYLOn       = False              ; no y left tickmarks
	res@tmXBOn       = False              ; no bottom tickmarks
	res@gsnPanelLabelBar   = False
	res@pmLabelBarWidthF   = 0.5 ; 0.5 
	res@pmLabelBarHeightF  = 0.06
	res@gsnCenterStringOrthogonalPosF = -0.06           ; move whole label bar down
	res@gsnLeftStringOrthogonalPosF = -0.05           ; move whole label bar down
	res@lbTitleString       = ""
	res@lbTitlePosition     = "Top"              ; title position
	res@lbTitleFontHeightF  = 0.03                ; make title smaller
	res@lbTitleDirection    = "Across"             ; title direction
	res@lbTitleExtentF=0.07
	res@lbTitleJust="TopCenter"
	res@tmXBLabelFontHeightF = 0.03
	res@tmYLLabelFontHeightF = 0.03
	res@lbLabelFontHeightF   = 0.03
	res@lbBoxLinesOn=False
	res@gsnMaximize          = True              ; Maximize plot
	res@cnInfoLabelOn       = False       ; turn off info label
	res@cnConstFLabelOn      = False
	res@cnConstFEnableFill   = True
	res@cnLevelSelectionMode = "ManualLevels"
	res@lbTitlePosition     = "Right"              ; title position
	res@lbTitleFontHeightF  = 0.03                ; make title smaller
	res@lbTitleDirection    = "Across"             ; title direction
	res@tmXBLabelFontHeightF = 0.03
	res@tmYLLabelFontHeightF = 0.03
	res@lbLabelFontHeightF   = 0.03	
	res@gsnCenterStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF   = 0.03
	res@gsnCenterStringOrthogonalPosF = 0.06          ; move center string down
	res@gsnLeftStringOrthogonalPosF = -0.05 
	res@gsnCenterString      = ""
	res@gsnRightString       = ""
	res@lbLabelBarOn          = True              ; no label bar 
	res@cnFillOn              = True               ; color on
	res@cnLinesOn             = False              ; no contour lines
	res@cnLineLabelsOn        = False              ; no line labels
	
	archivedir="./"
	patho=archivedir+"Figure3_Data/"
	
	;;; plot significant level
	res2 = True                            ; res2 probability plots
	res2@gsnDraw             = False       ; Do not draw plot
	res2@gsnFrame            = False       ; Do not advance frome
	res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
	res2@cnMinLevelValF      = 0.00        ; set min contour level
	res2@cnMaxLevelValF      = 1.05        ; set max contour level
	res2@cnLevelSpacingF     = 0.05        ; set contour spacing
	res2@cnInfoLabelOn       = False       ; turn off info label
	res2@cnLinesOn           = False       ; do not draw contour lines
	res2@cnLineLabelsOn      = False       ; do not draw contour labels
	res2@cnFillScaleF        = 0.5         ; add extra density, 0.6
	sig_lvl = 0.05
	opt     = True
	opt@gsnShadeFillType = "pattern"
	opt@gsnShadeLow = 10 ; 17 is dot, 10 is hatch
	opt@gsnShadeFillDotSizeF       = 0.003
	
	fname_string="_" + season	
	temperature_unit = "~F34~0~F~C"
	
	plotunits=(/temperature_unit,temperature_unit,temperature_unit, temperature_unit/)
	plotvars=(/"T-2m","T-2m","Max T-2m","Min T-2m"/)
	vars=(/"TSA","TSA","TREFMXAV","TREFMNAV" /)
	
	dsizes = dimsizes(vars)
 	nvars = dsizes(0)
	
	rcp_bases = (/"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP1-2.6Lu.101_reforestation",\
			"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP1-2.6Lu.102_reforestation",\
			"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP1-2.6Lu.103_reforestation"/)
			
	rcps = (/"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP2-2.6Lu.101_withbioenergy",\
			"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP2-2.6Lu.102_withbioenergy",\
			"b.e21.BSSP126cmip6.f09_g17.CMIP6-SSP1-2.6-SSP2-2.6Lu.103_withbioenergy"/)
	
	fname=patho+"TSA"+"_"+rcp_bases(0)+"_"+period+fname_string+".nc"
	print((/fname/))
	addfilercp = addfile(fname,"r")
	lats   = addfilercp->lat
	lons   = addfilercp->lon
	delete(addfilercp)
	
	ii = [/0,2,4,6/]
	do varID=0,nvars-1
		ploti = ii[varID]
		varname =vars(varID)
		
		forests = new((/3,dimsizes(lats),dimsizes(lons)/),float)
		bioenergys = new((/3,dimsizes(lats),dimsizes(lons)/),float)
		
		do ensi=0,2
			fname=patho+varname+"_"+rcp_bases(ensi)+"_"+period+fname_string+".nc"
			print((/fname/))
			addfilercp = addfile(fname,"r")
			reforest = (addfilercp->diff)
			forests(ensi,:,:) = reforest
		
			fname=patho+varname+"_"+rcps(ensi)+"_"+period+fname_string+".nc"
			print((/fname/))
			addfilercp = addfile(fname,"r")
			bioenergy = (addfilercp->diff)
			bioenergys(ensi,:,:) = bioenergy

			delete(bioenergy)
			delete(reforest)
		end do
		
		forest =  dim_avg_n_Wrap(forests,0)
		value =  dim_avg_n_Wrap(bioenergys,0)
		meanann_value=cal_diff(value,forest) 	
		boundary_array = set_boundary_bio_ref(varname) 
		res@cnMinLevelValF 	  = boundary_array[1]
		res@cnMaxLevelValF 	  = boundary_array[0]
		res@cnLevelSpacingF   = boundary_array[2]
		res@lbTitleString     = boundary_array[3]
		delete(res@cnFillColors)
		res@cnFillColors       = boundary_array[4]
		res@gsnCenterString        = ""
		res@gsnLeftString          = ""
		plot2(ploti)=gsn_csm_contour_map_ce(wks2,meanann_value, res)
		
		sigtests = new((/3,dimsizes(lats),dimsizes(lons)/),float)
		do ensi=0,2
			sigtest_one_ens=read_sigtest(patho,varname,rcps(ensi),period,fname_string)
			sigtests(ensi,:,:) = sigtest_one_ens
			delete(sigtest_one_ens)
		end do
		sigtest =  dim_avg_n_Wrap(sigtests,0)
		
		plot_prob   = gsn_csm_contour(wks2,sigtest, res2)
		plot_prob   = gsn_contour_shade(plot_prob, sig_lvl, 30, opt)  ; shade all areas less than the 0.05 contour level
		overlay (plot2(ploti), plot_prob)
		
		delete(meanann_value)
		delete(forest)
		delete(value)
		delete(sigtest)
		delete(forests)
		delete(bioenergys)
		delete(sigtests)
	end do

		
	delete(res)
	res                      =  True                           ;-- resource object
	res@xyLineThicknessF     =  3.0                            ;-- set line thickness
	res@xyLineColor          =  "black"                          ;-- set line color
	res@xyDashPatterns       =  0                              ;-- solid line
	res@gsnDraw               = False             ; don't draw
	res@gsnFrame              = False             ; don't advance frame
	res@vpWidthF         = 0.2            ; Change the aspect ratio, but 
	res@vpHeightF        = 0.49            ; make plot as large as possible.
	res@gsnXRefLine           = 0.0
	res@gsnPanelXWhiteSpacePercent = 5
	res@gsnPanelYWhiteSpacePercent = 3
	res@tmXTOn               = False                ; Turn off x-axis (top) major ticks
	res@tmYRMinorOn          = False                ; Turn off y-axis (right) minor ticks
	res@tmYROn               = False                ; Turn off y-axis (right) major ticks
	res@trYMinF 	  = -90
	res@trYMaxF 	  = 90
	res@tmYRLabelFontHeightF = 0.03
	res@tmXBLabelFontHeightF = 0.03
	res@tiXAxisFontHeightF =0.035
 	res@tmXBLabelStride = 3  ; set the x bottom tick interval
 	
	maxt=2.0
	mint=-2.0
	lower_boundarys=(/mint,mint,mint,mint/)
 	upper_boundarys=(/maxt,maxt,maxt,maxt/)
	ii= [/1,3,5,7/]
	do varID=0,3
		ploti = ii[varID]
		varname = vars(varID)
		unit=plotunits(varID)
		
		;;; mean annual time series
		fname_string="_meanAnnual"
		plotvalue1 = read_lat_mean(varname,lats,lons,rcp_bases,rcps,period,fname_string,patho,plotvars(varID) + " ("+unit+")")
		
		;;; JJA time series
		fname_string="_JJA"
		plotvalue2 = read_lat_mean(varname,lats,lons,rcp_bases,rcps,period,fname_string,patho,plotvars(varID) + " ("+unit+")")
		
		;;; DFJ time series
		fname_string="_DJF"
		plotvalue3 = read_lat_mean(varname,lats,lons,rcp_bases,rcps,period,fname_string,patho,plotvars(varID) + " ("+unit+")")
	
		res@xyLineColors = (/"black","orange","blue"/)
		plotvalues = new((/3,dimsizes(lats)/),float)
		plotvalues(0,:) = plotvalue1
		plotvalues(1,:) = plotvalue2
		plotvalues(2,:) = plotvalue3
	
		res@gsnCenterString           = ""
		res@gsnLeftString          = ""
		res@trXMinF 	  = lower_boundarys(varID)
		res@trXMaxF 	  = upper_boundarys(varID)
		plot2(ploti) = gsn_csm_xy(wks2, plotvalues,lats,  res)
		
		print("min annual, JJA, DJF for boreal")
		begi=158-30
		endi=158
		print(min(plotvalue1(begi:endi)))
		print(min(plotvalue2(begi:endi)))
		print(min(plotvalue3(begi:endi)))
		print("max")
		print(max(plotvalue1(begi:endi)))
		print(max(plotvalue2(begi:endi)))
		print(max(plotvalue3(begi:endi)))
		print("avg")
		print(avg(plotvalue1(begi:endi)))
		print(avg(plotvalue2(begi:endi)))
		print(avg(plotvalue3(begi:endi)))
		
		delete(plotvalues)
		delete(plotvalue1)
		delete(plotvalue2)
		delete(plotvalue3)
	end do
	
	
;************************************************
; create panel plot
;************************************************
    resPanel                    = True
	resPanel@gsnMaximize        = True                       ; use full page
	resPanel@gsnPanelLabelBar   = False                       ; add color bar
	resPanel@gsnPanelYWhiteSpacePercent = 20
	resPanel@gsnPanelXWhiteSpacePercent = 10
	resPanel@txString         = ""
	
	row=4
	col=2
	
	gsn_panel(wks2,plot2,(/row,col/),resPanel)
 end



