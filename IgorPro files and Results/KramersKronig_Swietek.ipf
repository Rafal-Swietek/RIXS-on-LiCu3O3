#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//	This function does the Kramer-Kronig transformation on given set of points defined in 'xdat'
// the wave 'inwave' is the wave on want to transform
//	The transformation is done over the given energy range, thus it not trustworthy
function Kramers_Kronig(inwave, p1, p2, C)
	wave inwave
	variable  p1, p2, C //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	//string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	//setdatafolder(filename)
		
		//filename = "RefIndex_Im"+num2str(p1)+"to"+num2str(p2)
		//wave inputwave = $filename //wave to transform and energy axis scaling
		
	variable E1 = pnt2x(inwave,0), dE =  pnt2x(inwave,1) - pnt2x(inwave,0)
	String newName = NameOfWave(inwave) + "_KK" // Compute output wave name
	Duplicate/O inwave, $newName // Create output wave
	Wave outwave = $newName // Create wave reference for output wave
	Duplicate/O inwave, xdat
	xdat[] = E1+dE*p
	
	Variable i, j, delx, dummyVar, eps
	eps = 1e-8
	//KKtrasform
	i = 0
	do
		dummyVar = 0
		//KK integral with cauchy principal value
		j = 1
		do
			if( abs(xdat[i] - xdat[j]) > eps)
				delx = xdat[j] - xdat[j-1]
				dummyVar += (xdat[j]*inwave[j]/(xdat[j]^2 - xdat[i]^2))*delx
			endif
			j += 1
		while(j < dimsize(inwave,0) )
		//-------------------------------------
		outwave[i] = C+2/pi*dummyVar  
		i += 1	
	while(i < dimsize(inwave,0) )
	
	//setdatafolder("root:")
	
end function


//inverse Kramers-Kronig transform
function Kramers_Kronig_inverse(inwave, p1, p2, C)
	wave inwave
	variable  p1, p2, C //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	//string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	//setdatafolder(filename)
		
		//filename = "RefIndex_Im"+num2str(p1)+"to"+num2str(p2)
		//wave inputwave = $filename //wave to transform and energy axis scaling
		
	variable E1 = pnt2x(inwave,0), dE =  pnt2x(inwave,1) - pnt2x(inwave,0)
	String newName = NameOfWave(inwave) + "_KK" // Compute output wave name
	Duplicate/O inwave, $newName // Create output wave
	Wave outwave = $newName // Create wave reference for output wave
	Duplicate/O inwave, xdat
	xdat[] = E1+dE*p
	
	Variable i, j, delx, dummyVar, eps
	eps = 1e-8

	//KKtrasform
	for(i=0; i <  dimsize(inwave,0); i+=1)
		dummyVar = 0	
		//KK integral with cauchy principal value
		for(j=1;j < dimsize(inwave,0); j+=1)
			if( abs(xdat[i] - xdat[j]) > eps)
				delx = xdat[j] - xdat[j-1]
				dummyVar += xdat[i]*(inwave[j]-C)/(xdat[j]^2 - xdat[i]^2) *delx
			endif
		endfor
		//-------------------------------------
		outwave[i] = -2/pi*dummyVar  	
	endfor
	
	//setdatafolder("root:")
end function

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

//This function does the Kramers-Kronig transform including extrapolation of the data at both ends
//the inputwave is the one to transform, the xdat wave, which includes the energy scaling is extrapolated with an uniform scaling
function KKtransform_extrapolated(inputwave, p1, p2, C,length_ext)
	wave inputwave
	variable  p1, p2, C, length_ext //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	//string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	//setdatafolder(filename)
		
		//filename = "RefIndex_Im"+num2str(p1)+"to"+num2str(p2)
		//wave inputwave = $filename //wave to transform and energy axis scaling
		
	variable E1 = pnt2x(inputwave,0), dE =  pnt2x(inputwave,1) - pnt2x(inputwave,0)
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave outwave = $newName // Create wave reference for output wave
	Duplicate/O inputwave, xdat, dummy
	xdat[] = E1+dE*p
	
	Variable i, j, k, delx, dummyVar, eps = 1e-8, num
	delx = xdat[1] - xdat[0]
	//Fitting data off range-------------------
		Variable length,length1,length2
		length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
		length2 = length_ext
		length = length1 + length2
		Make/O /N=4 left_fit_coef, left_fit_coef2 
		Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
		
		xdata[0,length1] = delx*p
		for(i=0;i<length2;i+=1)
			xdata[dimsize(xdat,0)+length1,dimsize(xdat,0)+length1+i] = xdat[dimsize(xdat,0)-1] + delx*(i+1)
		endfor
		setscale/i x 0, xdata[dimsize(xdata,0)-1], inwave, xdata
		
		CurveFit/Q=1 lor kwCWave=left_fit_coef,  inputwave[0,x2pnt(inputwave,939.5)] /X=xdat
		for(k=0;k<length1;k+=1)
			xdata[k] = k*delx
			inwave[k] = left_fit_coef[0] + left_fit_coef[1]/((xdata[k] - left_fit_coef[2])^2+left_fit_coef[3] )
			//inwave[k] = inwave[length1]
		endfor
		CurveFit/Q=1 lor kwCWave=left_fit_coef,  inputwave[x2pnt(inputwave,957.5),dimsize(xdat,0)-1] /X=xdat
		for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
			xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
			//inwave[k] = inputwave[dimsize(xdat,0)-1]
			inwave[k] = left_fit_coef[0] + left_fit_coef[1]/((xdata[k] - left_fit_coef[2])^2+left_fit_coef[3] )
		endfor
		for(k=length1; k<length1+dimsize(inputwave,0);k+=1)
			inwave[k] = inputwave[k - length1]
			xdata[k] = xdat[k-length1]
		endfor
	//display inwave
	//KKtrasform
	for(i=length1; i < length1 + dimsize(inputwave,0); i+=1)
		dummyVar = 0
		//Cauchy principal value
		for(j=1;j < dimsize(inwave,0); j+=1)
			if( abs(xdata[i] - xdata[j]) > eps)
				delx = xdata[j] - xdata[j-1]
				//dummyVar += (xdata[j]*inwave[j]-xdata[i]*inwave[i])/(xdata[j]^2 - xdata[i]^2)*delx // similiar KKT
				dummyVar += (xdata[j]*inwave[j])/(xdata[j]^2 - xdata[i]^2)*delx
			endif
		endfor
		dummy[i-length1] = C+2/pi*dummyVar
	endfor
	outwave[] = dummy[p]
	outwave[0] = outwave[1]
	killwaves/z inwave, right_fit_coef, left_fit_coef, xdata,W_sigma,W_fitConstants, xdat, W_coef, dummy

end function

//--------------------------------------------------------------------------------
//	This one does the inverse KKT in the same manner as above
function IKKtransform_extrapolated(inputwave, p1, p2, C, length_ext)  
	
	wave inputwave
	variable  p1, p2, C, length_ext //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	//string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	//setdatafolder(filename)
	
	//filename = "RefIndex_Re"+num2str(p1)+"to"+num2str(p2)+"_unwiggled"
	//wave inputwave = $filename //wave to transform and energy axis scaling
	
	variable E1 = pnt2x(inputwave,0), dE =  pnt2x(inputwave,1) - pnt2x(inputwave,0)
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave outwave= $newName // Create wave reference for output wave
	Duplicate/O inputwave, energy
	energy[] = E1+dE*p
	
	Duplicate/O inputwave, dummy, xdat
	dummy[] = 0
	interpolate2/F=1 /N=(dimsize(inputwave,0)) /Y = xdat energy
	
	Variable i, j, k, delx, dummyVar, eps =1e-8
	delx = xdat[1] - xdat[0]
	//Fitting data off range-------------------
	Variable length,length1,length2
	length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
	length2 = length_ext
	length = length1 + length2		 //right fit is exponential decay and left is Hillequation
	
	Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
	Make/O /N=4 fit_coef
	Curvefit/Q=1/NTHR=0 Sigmoid, kwCWave = fit_coef, inputwave[0,x2pnt(inputwave,936.6)] /X=xdat 
	for(k=0;k<length1;k+=1)
		xdata[k] = k*delx
		//inwave[k] = inputwave[0]
		inwave[k] = fit_coef[0] + fit_coef[1]/(1 + exp((fit_coef[2]-xdata[k])/fit_coef[3]) )
	endfor
	Make/O /N=2 fit_coef
	//Curvefit/Q=1/NTHR=0 line, kwCWave = fit_coef, inputwave[x2pnt(inputwave,967),dimsize(xdat,0)-1] /X=xdat //for LCP 575 commentize
	for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
		xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
		inwave[k] = inputwave[dimsize(inputwave,0)-1] //for LCP 575
		//inwave[k] = fit_coef[0] + fit_coef[1]*xdata[k]
	endfor
	for(k=length1; k<length1+dimsize(xdat,0);k+=1)
		inwave[k] = inputwave[k - length1]
		xdata[k] = xdat[k-length1]
	endfor
	setscale/i x 0, xdata[dimsize(xdata,0)-1], inwave, xdata
	//display inwave
	//KKtrasform
	for(i=length1; i < length1 + dimsize(inputwave,0); i+=1)
		dummyVar = 0
		//Cauchy principal value
		for(j=1;j < dimsize(inwave,0); j+=1)
			if( abs(xdata[i] - xdata[j]) > eps)
				delx = xdata[j] - xdata[j-1]
				dummyVar += (xdata[i]*(inwave[j]-C)/(xdata[j]^2 - xdata[i]^2))*delx 
			endif
		endfor
		dummy[i-length1] =  -2/pi*dummyVar 
	endfor
	outwave[] = dummy[p]
	outwave[0] = outwave[1]
	killwaves/z  inwave, right_fit_coef, left_fit_coef, xdata, ,W_sigma,W_fitConstants, xdat
	
end function
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This function fits the dielectric permittivity using Kramers-Kronig constrainted methods

function Comparer(p1,p2)
	variable p1, p2
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	filename = "Refindex_Im"+file
	
	
end function

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
function linearize_to1(inwave,E1,E2,var)
	wave inwave
	variable E1, E2, var
	Duplicate/O inwave, xdat
	variable i	
	for(i=0;i<dimsize(inwave,0);i+=1)
		xdat[i] = E1 + (E2-E1)/(dimsize(inwave,0)-1)*i
	endfor
	Duplicate/O inwave, inwave_line
	Curvefit/Q=1 line inwave /X=xdat /D=inwave_line
	variable sth = inwave[0]
	inwave[] = inwave[p]-inwave_line[p]+var
	killwaves/z  xdat, inwave_line
end function



function SubstractBack(dummy,E1,dE)
	wave dummy
	variable E1, dE
	Make/O/N = 2, test1, test2, test3
	Duplicate/O dummy, energy
	energy[] = E1 + dE*p
	variable L = dimsize(dummy,0)
	Curvefit/Q=1 line kwCWave=test1 dummy[0,x2pnt(energy,935)] /X=energy
	Curvefit/Q=1 line kwCWave=test2 dummy[x2pnt(energy,944),x2pnt(energy,956)] /X=energy
	Curvefit/Q=1 line kwCWave=test3 dummy[x2pnt(energy,962),L-1]  /X=energy 
	dummy[0,x2pnt(energy,939.4)] -= test1[1]*energy[p]+test1[0]
	dummy[x2pnt(energy,939.4)+1,x2pnt(energy,958.15)] -= test2[1]*energy[p]+test2[0]
	dummy[x2pnt(energy,958.15)+1,L-1] -= test3[1]*energy[p]+test3[0]
end function




Function Eps12(wavetotransform,energy)
	
	Wave wavetotransform, energy// Reference to the input wave received as parameter
	
	String thetaName = NameOfWave( wavetotransform) + "_theta" // Compute output wave name
	Duplicate/O wavetotransform, $thetaName // Create output wave
	Wave/c thetawave = $thetaName // Create wave reference for output wave
	redimension/c $thetaName
	

	String rName = NameOfWave( wavetotransform) + "_r" // Compute output wave name
	Duplicate/O wavetotransform, $rName // Create output wave
	Wave/c rwave = $rName // Create wave reference for output wave
	redimension/c $rName

	String eps1Name = NameOfWave( wavetotransform) + "_eps1" // Compute output wave name
	Duplicate/O wavetotransform, $eps1Name // Create output wave
	Wave eps1wave = $eps1Name // Create wave reference for output wave
	
	String eps2Name = NameOfWave( wavetotransform) + "_eps2" // Compute output wave name
	Duplicate/O wavetotransform, $eps2Name // Create output wave
	Wave eps2wave = $eps2Name // Create wave reference for output wave
	
	String ref1Name = NameOfWave( wavetotransform) + "_ref1" // Compute output wave name
	Duplicate/O wavetotransform, $ref1Name // Create output wave
	Wave ref1wave = $ref1Name // Create wave reference for output wave
	
	String ref2Name = NameOfWave( wavetotransform) + "_ref2" // Compute output wave name
	Duplicate/O wavetotransform, $ref2Name // Create output wave
	Wave ref2wave = $ref2Name // Create wave reference for output wave

	//Duplicate/O wavetotransform, dummy

	thetawave[]=cmplx(0.6,0)
	rwave[]=cmplx(0,0)
	eps1wave[]=0
	eps2wave[]=0
	
	variable i, k, delx
	
	for(k=0;k<Dimsize(wavetotransform,0);k+=1)
		for(i=1;i<dimsize(energy,0);i+=1)
			if(abs(energy[i]-energy[k])>0.000001)
				delx = energy[i] - energy[i-1]
				//check this formula!
				thetawave[k]+=-2*energy[k]/pi*ln(  sqrt(  wavetotransform[i]  ))/((energy[i])^2-(energy[k])^2)*delx
			endif
		endfor
	endfor
	
	//here you can play with the sign: -1 and +1 give totally different results
	thetawave[]*=-1
	rwave[0,dimsize(energy,0)-1]=sqrt(wavetotransform[p])*exp(cmplx(0,1)*thetawave[p])
	
	eps1wave[]=real((1-rwave[p])^2 / (1+rwave[p])^2)
	eps2wave[]=imag((1-rwave[p])^2 / (1+rwave[p])^2)
	ref1wave[] = sqrt( 0.5*( sqrt(eps1wave[p]^2+eps2wave[p]^2) + eps1wave[p]) )
	ref2wave[] = sqrt( 0.5*( sqrt(eps1wave[p]^2+eps2wave[p]^2) - eps1wave[p]) )
	
End

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This function analyses the influcence of length of the extrapolated region in KKT on the transform itself

function KK_vs_length(p1,p2)

	variable p1, p2
	string filename = "root:KKanalysis"
	setdatafolder(filename)	
	Make/O/N=12 Ls 
	Ls[] = {2,4,8,16,32,64,128,256,512,1024,2048,4096}
	//----------KKT	
		filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
		setdatafolder(filename)
		filename = "RefIndex_Im"+num2str(p1)+"to"+num2str(p2)
		wave inputwave = $filename //wave to transform and energy axis scaling
	
	filename = "root:KKanalysis"
	setdatafolder(filename)
	Duplicate/O inputwave, dummyR
	variable i
	for(i=1;i<10001;i+=100)
		String newName = "dummyR_"+num2str(i)//Ls[i]) // Compute output wave name
		Duplicate/O dummyR, $newName // Create output wave
		Wave outwave= $newName
		KKtransform_extrapolated(outwave, p1, p2, 1,i)//Ls[i])
		killwaves/z outwave
	endfor
	//----------IKKT
end



//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	 KK transform via double fast fourier transform with extrapolation; O(N)=2Nln(N), while trivial KKT is O(N)=N^2

function KKT_viaFFT(p1,p2,C,length_ext)
	variable  p1, p2, C, length_ext //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	setdatafolder(filename)
		
		filename = "RefIndex_Im"+num2str(p1)+"to"+num2str(p2)
		wave inputwave = $filename //wave to transform and energy axis scaling
		
	variable E1 = pnt2x(inputwave,0), dE =  pnt2x(inputwave,1) - pnt2x(inputwave,0)
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave outwave = $newName // Create wave reference for output wave
	Duplicate/O inputwave, xdat, dummy
	xdat[] = E1+dE*p
	
	Variable i, j, k, delx, dummyVar, eps = 1e-8, num
	delx = xdat[1] - xdat[0]
	//Fitting data off range-------------------
		Variable length,length1,length2
		length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
		length2 = length_ext
		length = length1 + length2
		Make/O /N=4 left_fit_coef, left_fit_coef2 
		Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
		
		xdata[0,length1] = delx*p
		for(i=0;i<length2;i+=1)
			xdata[dimsize(xdat,0)+length1,dimsize(xdat,0)+length1+i] = xdat[dimsize(xdat,0)-1] + delx*(i+1)
		endfor
		setscale/i x 0, xdata[dimsize(xdata,0)-1], inwave, xdata
		
		CurveFit/Q=1 lor kwCWave=left_fit_coef,  inputwave[0,x2pnt(inputwave,939.5)] /X=xdat
		for(k=0;k<length1;k+=1)
			xdata[k] = k*delx
			inwave[k] = left_fit_coef[0] + left_fit_coef[1]/((xdata[k] - left_fit_coef[2])^2+left_fit_coef[3] )
			//inwave[k] = inwave[length1]
		endfor
		CurveFit/Q=1 lor kwCWave=left_fit_coef,  inputwave[x2pnt(inputwave,957.5),dimsize(xdat,0)-1] /X=xdat
		for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
			xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
			//inwave[k] = inputwave[dimsize(xdat,0)-1]
			inwave[k] = left_fit_coef[0] + left_fit_coef[1]/((xdata[k] - left_fit_coef[2])^2+left_fit_coef[3] )
		endfor
		for(k=length1; k<length1+dimsize(inputwave,0);k+=1)
			inwave[k] = inputwave[k - length1]
			xdata[k] = xdat[k-length1]
		endfor
	//Algorithm
		filename = "root:KKT_viaFFT"
		setdatafolder(filename)
		Make/o/C/N=(dimsize(inwave,0)) dummy1
		Duplicate/c/O inwave, dummy2
		FFT/DEST = dummy1 inwave
		dummy1[] = -pi*cmplx(0,1)*dummy1[p] //change from initial formula
		IFFT/DEST = dummy2 dummy1
		dummy2[] = dummy2[p]/pi+C  //C just the same as for KKT normal
		KKtransform_extrapolated(inwave, p1, p2, 1,length_ext)
		setscale/i x 0, xdata[dimsize(xdata,0)-1], dummy2
		//killwaves/z dummy1, inwave
end

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	 IKK transform via double fast fourier transform with extrapolation; O(N)=2Nln(N), while trivial KKT is O(N)=N^2

function IKKT_viaFFT(p1,p2,C,length_ext)
	variable  p1, p2, C, length_ext //C=1 for Refreactive Index and C=0 for Optical Conductivity 
	string filename = "root:mda_"+num2str(p1)+"to"+num2str(p2)
	setdatafolder(filename)
		
		filename = "RefIndex_Re"+num2str(p1)+"to"+num2str(p2)
		wave inputwave = $filename //wave to transform and energy axis scaling
		
	variable E1 = pnt2x(inputwave,0), dE =  pnt2x(inputwave,1) - pnt2x(inputwave,0)
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave outwave = $newName // Create wave reference for output wave
	Duplicate/O inputwave, xdat, dummy
	xdat[] = E1+dE*p
	
	Variable i, j, k, delx, dummyVar, eps = 1e-8, num
	delx = xdat[1] - xdat[0]
	//Fitting data off range-------------------
		Variable length,length1,length2
		length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
		length2 = length_ext
		length = length1 + length2		 //right fit is exponential decay and left is Hillequation
		
		Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
		Make/O /N=4 fit_coef
		Curvefit/Q=1/NTHR=0 Sigmoid, kwCWave = fit_coef, inputwave[0,x2pnt(inputwave,936.6)] /X=xdat 
		for(k=0;k<length1;k+=1)
			xdata[k] = k*delx
			//inwave[k] = inputwave[0]
			inwave[k] = fit_coef[0] + fit_coef[1]/(1 + exp((fit_coef[2]-xdata[k])/fit_coef[3]) )
		endfor
		Make/O /N=2 fit_coef
		//Curvefit/Q=1/NTHR=0 line, kwCWave = fit_coef, inputwave[x2pnt(inputwave,967),dimsize(xdat,0)-1] /X=xdat //for LCP 575 commentize
		for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
			xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
			inwave[k] = inputwave[dimsize(inputwave,0)-1] //for LCP 575
			//inwave[k] = fit_coef[0] + fit_coef[1]*xdata[k]
		endfor
		for(k=length1; k<length1+dimsize(xdat,0);k+=1)
			inwave[k] = inputwave[k - length1]
			xdata[k] = xdat[k-length1]
		endfor
	//Algorithm
		filename = "root:IKKT_viaFFT"
		setdatafolder(filename)
		IKKtransform_extrapolated(inwave, p1, p2, 1,length_ext)
		Make/o/C/N=(dimsize(inwave,0)) dummy1
		Duplicate/c/O inwave, dummy2
		inwave[] = pi*(inwave[p] - C)
		FFT/DEST = dummy1 inwave
		dummy1[] = -dummy1[p]/(pi*cmplx(0,1)) //change from initial formula
		IFFT/DEST = dummy2 dummy1
		setscale/i x 0, xdata[dimsize(xdata,0)-1], dummy2
		//killwaves/z dummy1, inwave
end