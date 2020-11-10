#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function CreateGraph()
	setdatafolder("root:DielectricTensor:")
	string name = "xy"
	string nameRe = "EpsRe_"+name, nameIm = "EpsIm_"+name
	wave joe = $nameRe	
	display joe
	AppendToGraph/R=imag :EpsIm_xy
	ModifyGraph lblPosMode(imag)=1,freePos(imag)=0,axRGB(imag)=(65280,16384,16384)
	Legend/C/N=text0/J/A=RT/X=0.00/Y=0.00 "Optical Conductivity\r\\s("+nameRe+") Re \F'Symbol'e\F'Geneva'\B"+name+"\M\r\\s("+nameIm+") Im \F'Symbol'e\F'Geneva'\B"+name+"\M"
	Legend/C/N=text0/J/X=0./Y=0.00
	ModifyGraph lsize=1.15,rgb(EpsRe_xy)=(0,0,0)
end function

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	Here we check if CMXD is exhibited in the LiCu3O3 crystal

function MagneticDichroism(p1,p2,E1,p3,p4,E2)
	variable p1,p2,p3,p4,E1,E2
	
	//	Dichroism from imaginary part of refractzive index: change to TEY if necessary
	string dataname = "root:mda_" + num2str(p1)+"to"+num2str(p2)
	setdatafolder(dataname)
	//string filename = "XAS_"+num2str(p1)+"to"+num2str(p2)+"nobackground"
	string filename = "RefIndex_Im" + num2str(p1) + "to"+num2str(p2)
	wave joe = $filename
	
	dataname = "root:mda_" + num2str(p3)+"to"+num2str(p4)
	setdatafolder(dataname)
	//filename = "XAS_"+num2str(p3)+"to"+num2str(p4)+"nobackground"
	filename = "RefIndex_Im" + num2str(p3) + "to"+num2str(p4)
	wave jochen = $filename
	
	setdatafolder("root:")
	filename = "Dichroism"+num2str(p1)+"_"+num2str(p3)
	Duplicate/O joe, dummy1, test1
	Duplicate/O jochen, dummy2, test2
	Redimension/N=1000 dummy1, dummy2
	interpolate2/F=1 /N=1000 /Y=dummy1 test1
	interpolate2/F=1 /N=1000 /Y=dummy2 test2
	Duplicate/O dummy1, energy
	Duplicate/O dummy1, $filename
	wave dichro = $filename
	dichro[] = 0
	variable Eend = min( pnt2x(jochen,dimsize(jochen,0)-1), pnt2x(joe,dimsize(joe,0)-1))
	variable E = max(E1,E2), dE = (Eend-E)/999, i
	energy[] = E + dE*p
	//kappa to absorption
		variable h = 4.14e-15/2/pi, c = 3e8
		dummy1[] =2*energy[p]*dummy1[p]/c/h
		dummy2[] =2*energy[p]*dummy2[p]/c/h
	//----------------
	setscale/i  x,E, Eend, dichro
	//variable TEYarea = area(dummy1)
	//dummy1[] /= TEYarea
	//TEYarea = area(dummy2)
	//dummy2[] /= TEYarea
	for(i=0;i<dimsize(dichro,0);i+=1)
		dichro[i] = -dummy1[x2pnt(dummy1,energy[i])] + dummy2[x2pnt(dummy2,energy[i])]
	endfor
	display dummy1, dummy2, dichro
	ModifyGraph lsize=1.2,rgb(dummy2)=(0,0,0)
	Legend/C/N=text0/J/A=RT/X=0.00/Y=0.0 "\Z12\s(dummy1) \F'Symbol'm\F'Geneva' \BLCP\M \r\Z12\s(dummy2) \F'Symbol'm\F'Geneva' \BRCP\M ("+num2str(p3)+")"
	//Label left "Total Electron Yield (A.u.)";DelayUpdate
	Label left "Linear Absorption Coefficient \r\F'symbol'm\F'geneva' [1/m]";DelayUpdate
	Label bottom "photon energy[eV]"
	TextBox/C/N=text2/A=MT/X=0.00/Y=0.00 "\\Z14XMCD"
	TextBox/C/N=text1/A=MT/X=0.00/Y=0.00/E "\Z12Absorption from imaginary part of refractive index"
	//SavePICT/E=-4/B=576
	killwaves/z test1, test2, dummy1, dummy2
end
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

function KKanalysis(p1, p2, E1, dE)
	variable p1, p2, E1, dE
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	
	//Real part----------------------------------------------------------
	filename = "RefIndex_Im"+file
	wave joe = $filename
	Duplicate/O joe, dummy1, dummyIm
	dummy1[] = joe[p]
	
	KKtransform_extrapolated(dummy1,p1,p2,1,10000)
	//Kramers_Kronig(dummy1, p1, p2, 1)
	
	filename = "dummy1_KK"
	wave joe = $filename
	Duplicate/O joe, dummy1KK
	dummy1KK[] = joe[p]
	filename = "RefIndex_Re"+file+"_unwiggled"
	wave jochen = $filename
	dummy1[] = jochen[p]
	
	variable dummyVar = dummy1[0]
	//dummy1[] += 1 - dummyVar
	//dummyVar= dummy1KK[0]
	//dummy1KK[] += 1 - dummyVar
	display dummy1, dummy1KK
	ModifyGraph lsize=1.2,rgb(dummy1KK)=(0,0,0)
	ModifyGraph nticks(left)=10
	Legend/C/N=text0/J/A=RT/X=0.00/Y=0.0 "\Z12\\s(dummy1KK) Real part from Kramers-Kronig\r\\s(dummy1) Real part from peak position"
	Label left "n";DelayUpdate
	Label bottom "photon energy[eV]"
	TextBox/C/N=text2/A=LT/X=0.00/Y=0.00 "\\Z12LV("+num2str(p1)+")"
	//-------------------------------------------------------------------------
	//Imaginary part-----------------------------------------------------
	filename = "RefIndex_Re"+file
	wave joe = $filename
	Duplicate/O joe, dummy2
	dummy2[] = joe[p]
	
	IKKtransform_extrapolated(dummy2,p1,p2,1,10000)
	//Kramers_Kronig_inverse(dummy2, p1, p2, 1)
	
	filename = "dummy2_KK"
	wave joe = $filename
	Duplicate/O joe, dummy2KK
	dummy2[] = dummyIm[p]
	//dummyVar= dummy2[10]
	//dummy2[] -= dummyVar
	//dummyVar= dummy2KK[10]
	//dummy2KK[] -= dummyVar
	display dummy2, dummy2KK
	ModifyGraph lsize=1.2,rgb(dummy2KK)=(0,0,0)
	ModifyGraph nticks(left)=10
	Legend/C/N=text0/J/A=RT/X=0.00/Y=0.0 "\\Z12\s(dummy2KK) Imaginary part from Kramers-Kronig\r\\s(dummy2) Imaginary part from peak width"
	Label left "\F'Symbol'k";DelayUpdate
	Label bottom "photon energy[eV]"
	TextBox/C/N=text2/A=LT/X=0.00/Y=0.00 "\\Z12LV("+num2str(p1)+")"
	//-------------------------------------------------------------------------
	
	killwaves/z energy, dummyIm, dummy1, dummy1KK, dummy2, dummy2KK
	setdatafolder("root:")
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

function sth()
//	string filename = "Fit_"+num2str(575)+"to"+num2str(815)+"_abs"
//	wave joe =$filename
//	Duplicate/O joe, dummy1
//	SubstractBack(dummy1,925.15,0.15)
//	variable dummyVar= wavemax(dummy1)
//	dummy1[] /= dummyVar
//	
//	filename = "Fit_"+num2str(960)+"to"+num2str(1300)+"_abs"
//	wave joe =$filename
//	Duplicate/O joe, dummy2
//	SubstractBack(dummy2,925,0.15)
//	dummyVar= wavemax(dummy2)
//	dummy2[] /= dummyVar
//	
//	filename = "Fit_"+num2str(1324)+"to"+num2str(1689)+"_abs"
//	wave joe =$filename
//	Duplicate/O joe, dummy3
//	SubstractBack(dummy3,925.15,0.15)
//	dummyVar= wavemax(dummy3)
//	dummy3[] /= dummyVar
//	
//	filename = "Fit_"+num2str(1713)+"to"+num2str(2044)+"_abs"
//	wave joe =$filename
//	Duplicate/O joe, dummy4
//	SubstractBack(dummy4,926.65,0.15)
//	dummyVar= wavemax(dummy4)
//	dummy4[] /= dummyVar

	string filename = "Fit_"+num2str(1324)+"to"+num2str(1689)+"_abs"
	wave joe =$filename
	Duplicate/O joe, dummy1
	SubstractBack(dummy1,925.15,0.15)
	variable dummyVar= wavemax(dummy1)
	dummy1[] /= dummyVar
	
	filename = "Fit_"+num2str(2057)+"to"+num2str(2590)+"_abs"
	wave joe =$filename
	Duplicate/O joe, dummy2
	SubstractBack(dummy2,925,0.15)
	dummyVar= wavemax(dummy2)
	dummy2[] /= dummyVar
	//display dummy1, dummy2
	
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	This function evaluates the XAS spectrum for a given polarization (defined by the two variables) 
// as an average over the scattering angles

function XAS(p1, p2, E1, dE)
	variable p1, p2, E1, dE
	
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	filename = "D2_"+file
	wave joe = $filename
	filename = "Xas"+file
	Make/O/N=(dimsize(joe,1)) test
	Duplicate/O test, $filename
	wave xas = $filename
	setscale/p x, E1, dE, xas
	variable i
	Redimension/N=(dimsize(joe,0)) test
	for(i=0;i<dimsize(joe,1);i+=1)
		test[] = joe[p][i]
		xas[i] = faverage(test)
	endfor
	setdatafolder("root:")
	killwaves/z test
	
end function

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	This function creates a background to your TEY data,
// E1 is first energy, E2 is last (of the inverted TEY, so the not the very last but almost)
// This algorithm creates the background making on 3 different subsets of the xrange which are divided by Em1 and Em2
function CreateShirleyBackground(dummy,back, dE, E1, E2, Em1, Em2, Em3)
	wave dummy, back
	variable dE, E1, E2, Em1, Em2, Em3
	Duplicate/O dummy, energy
	energy[] = E1 + dE*p
	variable i, j, k, dummyVar= 0, B, kn
	// Iterative Shirley Algorithm : do more iterations!
		back[] = 0
		Variable num_of_iter = 200 //number of shirley iterations
		variable El = Em1, Er = E2 //energy range
		variable Ir = dummy[x2pnt(dummy,Er)]
		for(k=0;k<=num_of_iter;k+=1)
			dummyVar= 0
			for(j=x2pnt(dummy,El);j<x2pnt(dummy,Er);j+=1)
				dummyVar+= (dummy[j] - Ir - back[j])*dE //B0 = 0
			endfor
			kn = ( dummy[x2pnt(dummy,El)] - dummy[x2pnt(dummy,Er)] )/dummyVar//k_k
			for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
				B = 0
				for(j=i;j<x2pnt(dummy,Er);j+=1)
					B += (dummy[j] - Ir - back[j])*dE
				endfor
				back[i] = kn*B //B_k
			endfor
		endfor
		if(Em2!=0)
			Ir = back[x2pnt(back,El)] //changed dummy to back here !!!
			Er = El //energy range
			El = Em2 
			for(k=0;k<=num_of_iter;k+=1)
				dummyVar= 0
				for(j=x2pnt(dummy,El);j<x2pnt(dummy,Er);j+=1)
					dummyVar+= (dummy[j] - Ir - back[j])*dE //B0 = 0
				endfor
				kn = ( dummy[x2pnt(dummy,El)] - dummy[x2pnt(dummy,Er)] )/dummyVar//k_k
				for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
					B = 0
					for(j=i;j<x2pnt(dummy,Er);j+=1)
						B += (dummy[j] - Ir - back[j])*dE
					endfor
					back[i] = kn*B //B_k
				endfor
			endfor
			for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
				back[i] += Ir //B_k
			endfor
		endif
		if(Em3!=0)
			Ir = back[x2pnt(back,El)] //changed dummy to back here !!!
			Er = El //energy range
			El = Em3
			for(k=0;k<=num_of_iter;k+=1)
				dummyVar= 0
				for(j=x2pnt(dummy,El);j<x2pnt(dummy,Er);j+=1)
					dummyVar+= (dummy[j] - Ir - back[j])*dE //B0 = 0
				endfor
				kn = ( dummy[x2pnt(dummy,El)] - dummy[x2pnt(dummy,Er)] )/dummyVar//k_k
				for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
					B = 0
					for(j=i;j<x2pnt(dummy,Er);j+=1)
						B += (dummy[j] - Ir - back[j])*dE
					endfor
					back[i] = kn*B //B_k
				endfor
			endfor
			for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
				back[i] += Ir //B_k
			endfor
		endif
		Make/O/N=2 test1
		Curvefit/Q=1 line kwCWave=test1 dummy[x2pnt(energy,E1),x2pnt(energy,El)] /X=energy
		for(i=0;i<x2pnt(dummy,El);i+=1)
			back[i] = test1[0] + test1[1]*energy[i] //back[x2pnt(back,El)] //B_k
		endfor
		Duplicate/O back, test
		back[] = test[dimsize(test,0)-p-1]
	//-----------------------------------------
	killwaves/z test, test1, energy
end function

function CreateShirleyBackground_1peak(dummy,back, dE, E1, E2)
	wave dummy, back
	variable dE, E1, E2
	Duplicate/O dummy, energy
	energy[] = E1 + dE*p
	variable i, j, k, dummyVar= 0, B, kn
	// Iterative Shirley Algorithm : do more iterations!
		back[] = 0
		Variable num_of_iter = 200 //number of shirley iterations
		variable El = E1, Er = E2 //energy range
		variable Ir = dummy[x2pnt(dummy,Er)]
		for(k=0;k<=num_of_iter;k+=1)
			dummyVar= 0
			for(j=x2pnt(dummy,El);j<x2pnt(dummy,Er);j+=1)
				dummyVar+= (dummy[j] - Ir - back[j])*dE //B0 = 0
			endfor
			kn = ( dummy[x2pnt(dummy,El)] - dummy[x2pnt(dummy,Er)] )/dummyVar//k_k
			for(i=x2pnt(dummy,El);i<x2pnt(dummy,Er);i+=1)
				B = 0
				for(j=i;j<x2pnt(dummy,Er);j+=1)
					B += (dummy[j] - Ir - back[j])*dE
				endfor
				back[i] = kn*B //B_k
			endfor
		endfor
	//-----------------------------------------
	killwaves/z test, test1, energy
end function

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	This function returns the polarization of the incoming, and outgoing wave wave
function switchPolarization(polarization, theta, p1)
	wave/c polarization
	variable p1, theta
	switch(p1)
		case 575: //LCP
			polarization[0] = 1
			polarization[1] = cmplx(0,1)*sin(theta) //[x2pnt(theta,omega[p])]
			polarization[2] = cmplx(0,1)*cos(theta)
			polarization[] /= sqrt(2)
			break;
		case 960: //LV
			polarization[0] = 1
			break;	
		case 1324: //RCP	
			polarization[0] = 1
			polarization[1] = -cmplx(0,1)*sin(theta)
			polarization[2] = -cmplx(0,1)*cos(theta)
			polarization[] /= sqrt(2)
			break;
		case 1713: //LH
			polarization[1] = sin(theta)
			polarization[2] = cos(theta)
			break;
		case 2057: //RCP	
			polarization[0] = 1
			polarization[1] = -cmplx(0,1)*sin(theta)
			polarization[2] = -cmplx(0,1)*cos(theta)
			polarization[] /= sqrt(2)
			break;
		case 2595: //LV
			polarization[0] = 1
			break;					
	endswitch
	killwaves/z energy, num
end function

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Creating RIXS intensity along Bragg peak
function RIXSonpeak(p1,p2)
	variable p1, p2
	string file = num2str(p1) + "to" + num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	
	filename = "Fit_" + file + "_"
	wave x0 = $filename+"x0"
	wave RIXS = $"D3_"+file
	filename = "RIXS_onpeak_"+file
	Duplicate/O x0, $filename, energy
	wave rixs_onpeak = $filename
	variable num = dimsize(x0,0), i
	variable E1 = pnt2x(x0,0), dE = pnt2x(x0,1) - pnt2x(x0,0)
	energy[] = E1 + dE*p
	for(i=0;i<num;i+=1)
		rixs_onpeak[i] = RIXS[x2pnt(RIXS,x0[i])][i]
	endfor
	setdatafolder("root:")
end function
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
