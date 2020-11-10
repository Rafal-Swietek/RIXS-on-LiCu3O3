#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//Function Extract_MetaData(basename, suffix, first, last)//Basename="EA_", suffix="avgy" or ""
	//string Basename, suffix
	//variable first, last
	//setdatafolder root:
	//make/o/n=(abs(last-first+1)) $("ScanNum_"+num2str(first)+"_"+num2str(last)), $("TEY_"+num2str(first)+"_"+num2str(last)), $("hv_"+num2str(first)+"_"+num2str(last))
	//wave ScanNum=$("ScanNum_"+num2str(first)+"_"+num2str(last))
	//wave TEY=$("TEY_"+num2str(first)+"_"+num2str(last))
	//wave hv=$("hv_"+num2str(first)+"_"+num2str(last))
	//ScanNum=first+p
	//ncNote2wave("Attr_TEY", NameofWave(TEY), basename,suffix,NameofWave(ScanNum))
	//ncNote2wave("Attr_ActualPhotonEnergy", NameofWave(hv), basename,suffix,NameofWave(scannum))
	//end
	
//Function Stack_Dither(basename, suffix, first,last)
	//string basename, suffix
	//variable first, last
	//variable i
	//for(i=0;i<abs(last-first+1);i+=1)
		//wave wv=$WaveNamewithNum(basename,first+i,suffix)
		//if(i==0)
			//duplicate/o wv Dither
			//redimension/n=(-1,dimsize(Dither,1)+last-first+1) Dither
			
		//else
			//wave dither= Dither
			//dither[][]+=wv[p][q+i]
		//endif
	//endfor
	
//End


function joe(first,last,angle1,angle2,E1,deltaE,diode,inc)

	variable first, last, angle1, angle2, E1, deltaE,inc
	
	string diode 
	string filename = diode+"_"+num2str(first)+"to"+num2str(last)
	
		
			if(first<1000)
			setdatafolder("root:mda_0"+num2str(first))
			else
			setdatafolder("root:mda_"+num2str(first))
			//print "root:mda_"+num2str(first)
			endif
			
			make/O/N=(1001,(last-first)/inc+1) $filename
			setscale/p x, angle1, (angle2-angle1)/1000, $filename
			setscale/p y, E1, deltaE, $filename
			wave jochen =$filename
	
	variable i
	
	for(i=first;i<=last;i+=inc)
			
			if(i<1000)
			setdatafolder("root:mda_0"+num2str(i))
			//print "root:mda_0"+num2str(i)
			else
			setdatafolder("root:mda_"+num2str(i))
			//print "root:mda_0"+num2str(i)
			endif
			
			make/o/n=1001 xrange
			setscale/p x, angle1, (angle2-angle1)/1000, xrange
			xrange[]=x
			
			//Interpolate2/T=3/N=1001/F=0/Y=$"c29idd_ca3_read_SS" $"c29idd_m7_VAL", $"c29idd_ca3_read"
			Interpolate2/T=3/N=1001/I=3/F=0/Y=$diode+"_SS"/X=xrange $"c29idd_m7_VAL", $diode

			
			wave interpolate =$diode+"_SS"
			
			
			//root:mda_0498:jochen[][i-first]=interpolate[p]
			jochen[][(i-first)/inc]=interpolate[p]

	endfor
	
	
		
			if(first<1000)
			setdatafolder("root:mda_0"+num2str(first))
			else
			setdatafolder("root:mda_"+num2str(first))
			//print "root:mda_"+num2str(first)
			endif
		
end




function loadRSXS(first,last,angle,dangle,E1,deltaE)

	variable first, last, angle, dangle, E1, deltaE
	
	
	string filename = "D3_"+num2str(first)+"to"+num2str(last)
	
	
			wave c29idd_ca3_read
			
			setdatafolder("root:mda_"+num2str(first))	
			make/O/N=(dimsize($"c29idd_ca3_read",0),(last-first)+1) $filename
			//print dimsize($"c29idd_ca3_read",0)
			setscale/i x, angle, dangle, $filename
			setscale/p y, E1, deltaE, $filename
			wave jochen =$filename
			//print dimsize(jochen,0)
			jochen[][] = 0
	variable i
	
	for(i=first;i<=last;i+=1)
			
			if(i<1000)
			setdatafolder("root:mda_0"+num2str(i))
			else
			setdatafolder("root:mda_"+num2str(i))
			
			endif
			
			
			wave joe = $"c29idd_ca3_read"
			wave normwave = $"c29idb_ca14_read"
			
			//print dimsize(normwave,0)
			
			//if (dimsize(joe,0)!=80)
			//print "root:mda_0"+num2str(i)
			//print dimsize(joe,0)
			//endif
			//Interpolate2/T=3/N=1001/I=3/F=0/Y=$diode+"_SS"/X=xrange $"c29idd_m7_VAL", $diode

			jochen[][i-first]= joe[p]/normwave[p]
			//print i
	endfor
	

		
	setdatafolder("root:mda_"+num2str(first))
end

function loadRSXS2(first,last,angle,dangle,E1,deltaE)

	variable first, last, angle, dangle, E1, deltaE
	
	
	string filename = "D14_"+num2str(first)+"to"+num2str(last)
	
	
			wave c29idd_ca14_read
			
			setdatafolder("root:mda_"+num2str(first))	
			make/O/N=(dimsize($"c29idb_ca14_read",0),(last-first)+1) $filename
			//print dimsize($"c29idd_ca3_read",0)
			setscale/i x, angle, dangle, $filename
			setscale/p y, E1, deltaE, $filename
			wave jochen =$filename
			//print dimsize(jochen,0)
			jochen[][] = 0
	variable i
	
	for(i=first;i<=last;i+=1)
			
			if(i<1000)
			setdatafolder("root:mda_0"+num2str(i))
			else
			setdatafolder("root:mda_"+num2str(i))
			
			endif
			
			
			wave joe = $"c29idb_ca14_read"
			//wave normwave = $"c29idb_ca14_read"
			
			//print dimsize(normwave,0)
			
			//if (dimsize(joe,0)!=80)
			//print "root:mda_0"+num2str(i)
			//print dimsize(joe,0)
			//endif
			//Interpolate2/T=3/N=1001/I=3/F=0/Y=$diode+"_SS"/X=xrange $"c29idd_m7_VAL", $diode

			jochen[][i-first]= joe[p]///normwave[p]
			//print i
	endfor
	

		
	setdatafolder("root:mda_"+num2str(first))
end



function loadRSXSall(first,last,angle,dangle,E1,deltaE)

	variable first, last, angle, dangle, E1, deltaE
	
	
	string filename = "D3_"+num2str(first)+"to"+num2str(last)
	string filename2 = "D2_"+num2str(first)+"to"+num2str(last)
	string filename_cur = "Dcur_"+num2str(first)+"to"+num2str(last)
	string filename2norm = "D2_"+num2str(first)+"to"+num2str(last)+"_norm"
	string filenamenorm = "D3_"+num2str(first)+"to"+num2str(last)+"_norm"
	string filename14 = "D14_"+num2str(first)+"to"+num2str(last)+"_norm"
	
	
			wave c29idd_ca2_read, c29idd_ca3_read, cS_SRcurrentAI_VAL, c29idb_ca14_read
			variable i
			
			if(first<1000)
				setdatafolder("root:mda_0"+num2str(first))
			else
				setdatafolder("root:mda_"+num2str(first))	
			endif	
			make/O/N=(dimsize($"c29idd_ca3_read",0),(last-first)+1) $filename
			make/O/N=(dimsize($"c29idd_ca2_read",0),(last-first)+1) $filename2
			make/O/N=(dimsize($"cS_SRcurrentAI_VAL",0),(last-first)+1) $filename_cur
			make/O/N=(dimsize($"c29idd_ca3_read",0),(last-first)+1) $filenamenorm
			make/O/N=(dimsize($"c29idd_ca2_read",0),(last-first)+1) $filename2norm
			make/O/N=(dimsize($"c29idb_ca14_read",0),(last-first)+1) $filename14
			//print dimsize($"c29idd_ca3_read",0)
			setscale/i x, angle, dangle, $filename
			setscale/p y, E1, deltaE, $filename
			setscale/i x, angle, dangle, $filename2
			setscale/p y, E1, deltaE, $filename2
			setscale/i x, angle, dangle, $filenamenorm
			setscale/p y, E1, deltaE, $filenamenorm
			setscale/i x, angle, dangle, $filename2norm
			setscale/p y, E1, deltaE, $filename2norm
			setscale/i x, angle, dangle, $filename_cur
			setscale/p y, E1, deltaE, $filename_cur
			setscale/i x, angle, dangle, $filename14
			setscale/p y, E1, deltaE, $filename14
			wave jochen1 =$filename
			wave jochen2 =$filename2
			wave jochen3 =$filename_cur
			wave jochen4 =$filenamenorm
			wave jochen5 =$filename2norm
			wave jochen6 =$filename14
			//print dimsize(jochen,0)
			jochen1[][] = 0
			jochen2[][] = 0
			jochen3[][] = 0
			jochen4[][] = 0
			jochen5[][] = 0
			jochen6[][] = 0
	
	for(i=first;i<=last;i+=1)
			
			if(i<1000)
			setdatafolder("root:mda_0"+num2str(i))
			else
			setdatafolder("root:mda_"+num2str(i))
			
			endif
			
			
			wave joe1 = $"c29idd_ca3_read"
			wave joe2 = $"c29idd_ca2_read"
			wave current = $"cS_SRcurrentAI_VAL"
			wave normwave = $"c29idb_ca14_read"
			
			//print dimsize(normwave,0)
			
			//if (dimsize(joe,0)!=80)
			//print "root:mda_0"+num2str(i)
			//print dimsize(joe,0)
			//endif
			//Interpolate2/T=3/N=1001/I=3/F=0/Y=$diode+"_SS"/X=xrange $"c29idd_m7_VAL", $diode

			jochen1[][i-first]= joe1[p]/normwave[p]
			jochen2[][i-first]= joe2[p]/normwave[p]
			jochen3[][i-first]= current[p]
			jochen4[][i-first]= joe1[p]/normwave[p]/current[p]
			jochen5[][i-first]= joe2[p]/normwave[p]/current[p]
			jochen6[][i-first]= normwave[p]
			//print i
	endfor
		
	if(first<1000)
		setdatafolder("root:mda_0"+num2str(first))
	else
		setdatafolder("root:mda_"+num2str(first))
	endif
	 //NewImageTool5(filename)
	 //NewImageTool5(filename2)
	 //NewImageTool5(filenamenorm)
	 //NewImageTool5(filename2norm)
	 //NewImageTool5(filename_cur)
	 //NewImageTool5(filename14)
	 killwaves/z c29idd_ca2_read, c29idd_ca3_read, cS_SRcurrentAI_VAL, c29idb_ca14_read
	 killwaves/z joe1, joe2
	 setdatafolder("root")
end