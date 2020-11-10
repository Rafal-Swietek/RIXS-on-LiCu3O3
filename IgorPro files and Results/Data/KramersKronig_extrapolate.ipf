#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//This function does the Kramers-Kronig transform including extrapolation of the data at both ends
//the inputwave is the one to transform, the xdat wave, which includes the energy scaling is extrapolated with an uniform scaling
function KKtransform_extrapolated(inputwave, xdat)

	wave inputwave, xdat 	//wave to transform and energy axis scaling
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave inputwave = $newName // Create wave reference for output wave

	Duplicate/O inputwave, outwave
	outwave[] = 0
	
	Variable i, j, k, delx, temp, start, stop, eps = 1e-8
	delx = xdat[1] - xdat[0]
	start = DateTime
	//Fitting data off range-------------------
	Variable length,length1,length2
	length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
	length2 = 100000
	length = length1 + length2
	Make/O /N=3 right_fit_coef 	//right fit is exponential decay and left is Hillequation
	Make/O /N=4 left_fit_coef 
	
	Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
	
	CurveFit/NTHR=0 HillEquation kwCWave=left_fit_coef,  inputwave[0,50] /X=xdat 
	for(k=0;k<length;k+=1)
		xdata[k] = k*delx
		inwave[k] = left_fit_coef[0] + (left_fit_coef[1]-left_fit_coef[0])/(1+(left_fit_coef[3]/xdata[k])^left_fit_coef[2])		
	endfor
	
	Curvefit/NTHR=0 exp_XOffset, kwCWave = right_fit_coef, inputwave[dimsize(xdat,0)-29,dimsize(xdat,0)] /X=xdat 
	Variable x0 = xdat[dimsize(inputwave,0)-29]
	for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
		xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
		inwave[k] = right_fit_coef[0] + right_fit_coef[1]*exp(-(xdata[k] - x0)/right_fit_coef[2] )
	endfor
	for(k=length1; k<length1+dimsize(xdat,0);k+=1)
		inwave[k] = inputwave[k - length1]
		xdata[k] = xdat[k-length1]
	endfor
	//----------------------------------------------
	//KKtrasform
	for(i=length1; i < length1 + dimsize(inputwave,0); i+=1)
		temp = 0
		//Cauchy principal value
		for(j=1;j < dimsize(inwave,0); j+=1)
			if( abs(xdata[i] - xdata[j]) > eps)
				delx = xdata[j] - xdata[j-1]
				temp += (xdata[j]*inwave[j]/(xdata[j]^2 - xdata[i]^2))*delx //rectangular integration
				// for use traezoid integration add anothe condition xdat[j-1]!=xdat[i]
				//temp +=0.5*( xdata[j]*inwave[j]/(xdata[j]^2 - xdata[i]^2) + xdata[j-1]*inwave[j-1]/(xdata[j-1]^2 - xdata[i]^2) )*delx //trapezoidal integration
			endif
		endfor
		outwave[i-length1] = 1+2/pi*temp
	endfor
	inputwave[] = outwave[p]
	killwaves/z outwave, inwave, right_fit_coef, left_fit_coef, xdata,W_sigma,W_fitConstants
	stop = DateTime
	print "Execution time: " + num2str(stop - start) +" seconds"
	
end function



function IKKtransform_extrapolated(inputwave, xdat)

	wave inputwave, xdat 	//wave to transform and energy axis scaling
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Wave inputwave = $newName // Create wave reference for output wave

	Duplicate/O inputwave, outwave
	outwave[] = 0
	
	Variable i, j, k, delx, temp, start, stop, eps =1e-8
	delx = xdat[1] - xdat[0]
	start = DateTime
	//Fitting data off range-------------------
	Variable length,length1,length2
	length1 = ceil(xdat[0]/delx)		// length of inputwave extension to zero
	length2 = 100000
	length = length1 + length2	
	Make/O /N=4 right_fit_coef 	//right fit is sigmoid and left is lorentzian
	Make/O /N=4 left_fit_coef 
	
	Make/O /N=(dimsize(inputwave,0)+length) inwave, xdata 	 //resized input wave and energy data wave
	
	CurveFit lor kwCWave=left_fit_coef,  inputwave[0,50] /X=xdat 
	Variable x0 = xdat[0]
	for(k=0;k<length1;k+=1)
		xdata[k] = k*delx
		inwave[k] = left_fit_coef[0] + left_fit_coef[1]/((xdata[k] - left_fit_coef[2])^2+left_fit_coef[3] )
	endfor
	
	Curvefit Sigmoid, kwCWave = right_fit_coef, inputwave[dimsize(inputwave,0)-29,dimsize(inputwave,0)] /X=xdat 
	for(k=length1+dimsize(xdat,0);k<dimsize(inwave,0);k+=1)
		xdata[k] = xdat[dimsize(xdat,0)-1] + (k+1-length1-dimsize(xdat,0))*delx
		inwave[k] = right_fit_coef[0] + right_fit_coef[1]/(1+exp( (right_fit_coef[2]-xdata[k])/right_fit_coef[3] ) )	
	endfor
	for(k=length1; k<length1+dimsize(inputwave,0);k+=1)
		inwave[k] = inputwave[k - length1]
		xdata[k] = xdat[k-length1]
	endfor
	//----------------------------------------------
	//KKtrasform
	for(i=length1; i < length1 + dimsize(inputwave,0); i+=1)
		temp = 0
		//Cauchy principal value
		for(j=1;j < dimsize(inwave,0); j+=1)
			if( abs(xdata[i] - xdata[j]) > eps)
				delx = xdata[j] - xdata[j-1]
				temp += (xdata[i]*inwave[j]/(xdata[j]^2 - xdata[i]^2))*delx  //rectangular integration
			endif
		endfor
		outwave[i-length1] =  -2/pi*temp 
	endfor
	inputwave[] = outwave[p]
	killwaves/z outwave, inwave, right_fit_coef, left_fit_coef, xdata, ,W_sigma,W_fitConstants
	stop = DateTime
	print "Execution time: " + num2str(stop - start) +" seconds"
	
end function



function linearize_to1(inwave,xdata)
	wave inwave, xdata
	Duplicate/O inwave, xdat
	variable i	
	for(i=0;i<dimsize(inwave,0);i+=1)
		xdat[i] = xdata[i]
	endfor
	Duplicate/O inwave, inwave_line
	Curvefit line inwave /X=xdat /D=inwave_line

	for(i=0;i<dimsize(inwave,0);i+=1)
		inwave[i] = inwave[i]/inwave_line[i]
	endfor
	killwaves/z inwave_line, xdata, xdat
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


function RotatebyAngle2D(theta, inwave)
//This function rotates a wave counterclockwise by an angle theta
	wave inwave
	variable theta
	Duplicate/O inwave, dumm
	Make/O/N=(2,2) RotationMatrix
	RotationMatrix[0][0] = cos(theta)
	RotationMatrix[1][0] = -sin(theta)
	RotationMatrix[0][1] = sin(theta)
	RotationMatrix[1][1] = cos(theta)
	Variable i, j
	for(i=0;i<2;i+=1)
		for(j=0;j<2;j+=1)
			dumm[i] = RotationMatrix[i][j]*inwave[j]
		endfor
	endfor
	inwave[] = dumm[p]
	killwaves/z dumm, RotationMatrix
	return inwave
end function



function CreateFunEq(intensity, eps1, eps2, polarization, theta, sxx, sxy, sxz, szz, k)	
	wave intensity, eps1, eps2, theta
	string polarization
	variable sxx, sxy, sxz, szz, k // optical conductivity independent tensor elements
	
	Make/C/N=(dimsize(eps1,0),3) E_in, E_out //0,1,2 is x,y,z
	
	string newName = NameOfWave(intensity) + "_fun"
	Duplicate/O intensity, $newName
	wave intensity = $newName //reference for output wave
	E_in[][] = cmplx(0,0)
	E_out[][] = cmplx(0,0)
	variable dummy
	strswitch(polarization) //Incoming wave is in the yz - plane, propagating in the z-direction
		case "s": // E = |E|*[1,0,0]
			E_in[][0] = eps1[p] + cmplx(0,1)*eps2[p]
			E_out[][0] = E_in[p][0]//*reflectivity[p]
				dummy = intensity[k]
				dummy -= E_out[k][0]*sxx*E_in[k][0] 
				return dummy
			break;
		case "p":  // E = |E|*[0,1,0]
			E_in[][1] = eps1[p]  + cmplx(0,1)*eps2[p] 
			E_out[][1] = E_in[p][1]
				dummy = intensity[k]
				Make/N=3 E_in2
				E_in2[] = E_in[k][p]
				dummy -= E_out[k][1]*sxx*RotatebyAngle2D(2*theta[k],E_in2) //p-pol wave turns by 2theta after reflection
				return dummy
			break;
		case " RCP": //E = |E|*[0,1,-i]
			E_in[][1] = cmplx(eps1[p] ,eps2[p] )
			E_in[][2] = -cmplx(0,1)*cmplx(eps1[p] ,eps2[p] )
			E_out[][1] = E_in[p][1]
			E_out[][2] = E_in[p][2]
				dummy = intensity[k]
				dummy -= E_out[k][1]*sxx*E_in[k][1] + E_out[k][1]*sxz*E_in[k][2]+E_out[k][2]*sxz*E_in[k][1]+E_out[k][2]*szz*E_in[k][2]
				return dummy
			break;
		case "LCP": //E = |E|*[0,1,i]
			E_in[][1] = cmplx(eps1[p] ,eps2[p] )
			E_in[][2] = cmplx(0,1)*cmplx(eps1[p] ,eps2[p] )
			E_out[][1] = E_in[p][1]
			E_out[][2] = E_in[p][2]
				dummy = intensity[k]
				dummy -= E_out[k][1]*sxx*E_in[k][1] + E_out[k][1]*sxz*E_in[k][2]+E_out[k][2]*sxz*E_in[k][1]+E_out[k][2]*szz*E_in[k][2]
				return dummy
			break;	
	endswitch
	//for(k=0;j<dimsize(intensity,0);j+=1)
	//endfor
	print 33
	killwaves/z E_in, E_out, dummy
end function


function OpticalCondExtract(intensity, eps1, eps2, theta,energy)
	wave intensity, eps1, eps2, theta, energy
	variable sxx, sxy, sxz, szz
	variable i, k, fun1, fun2, fun3, fun4
	Make/N=(dimsize(intensity,1))  epsre1, epsre2, epsre3, epsre4, epsim1, epsim2, epsim3, epsim4
	Make/N=(dimsize(intensity,1)) intensity1, intensity2, intensity3, intensity4
	epsre1[] = eps1[p][0]
	epsre2[] = eps1[p][1]
	epsre3[] = eps1[p][2]
	epsre4[] = eps1[p][3]
	epsim1[] = eps2[p][0]
	epsim2[] = eps2[p][1]
	epsim3[] = eps2[p][2]
	epsim4[] = eps2[p][3]
	intensity1[] = intensity[p][0]
	intensity2[] = intensity[p][1]
	intensity3[] = intensity[p][2]
	intensity4[] = intensity[p][3]
	Make/N=3 cw1, cw2, cw3
	cw1[] = 0
	cw2[] = 0
	cw3[] = 0
	for(i=0;i<dimsize(energy,0);i+=1)
		fun1 = CreateFunEq(intensity1, epsre1, epsim1, "s", theta, sxx, sxy, sxz, szz, i)
		fun2 = CreateFunEq(intensity2, epsre2, epsim2, "RCP", theta, sxx, sxy, sxz, szz, i)
		fun3 = CreateFunEq(intensity3, epsre3, epsim3, "LCP", theta, sxx, sxy, sxz, szz, i)
		fun4 = CreateFunEq(intensity4, epsre4, epsim4, "RCP", theta, sxx, sxy, sxz, szz, i)
		FindRoots/Q fun1, cw1,fun2, cw2,fun3, cw3
		
	endfor
	 
	 
	 killwaves/z  epsre1, epsre2, epsre3, epsre4, epsim1, epsim2, epsim3, epsim4, intensity1, intensity2, intensity3, intensity4
end function

function RemoveWiggle(wavetoremove, wiggle, energy)
	wave wavetoremove, wiggle, energy
	
	string filename = NameOfWave(wavetoremove) + "_unwiggled"
	Make/O/N=(dimsize(wavetoremove,0)) dummy, data
	wave output = $filename
	Duplicate/O energy, energy2
	Duplicate/O wiggle, energy00, Emess //,line
	//Curvefit line wiggle /D=line
	//wiggle[] = wiggle[p]/line[p] - 1
	Emess[] = 917.5 + 0.15*p + wiggle[p]
	energy00[] = 917.5 + 0.15*p
	
	Duplicate/O wavetoremove, waveline
	linearize_to1(waveline, energy)
	filename = NameOfWave(waveline)
	wave waveline = $filename
	variable i, j, chi2, chi2previous=1e5, temp, delE, Ereal, E
	//Interpolation of data
	Variable num = 3000 //number of points for interpolation
	interpolate2 /F=1 /N=(num) /X = energy2 /Y = dummy energy, wavetoremove
	//interpolate2 /F=1 /N=(num)  /Y = energy2 energy
	//
	for(delE=0;delE<dimsize(wiggle,0);delE+=1)	
		chi2 = 0
		for(i=0;i<dimsize(wavetoremove,0);i+=1)
			E = energy[i]+delE
			Ereal = Emess[x2pnt(Emess,E)]
			data[i] = wavetoremove[x2pnt(wavetoremove,Ereal)]
			chi2 += (wavetoremove[i]-waveline[i])^2/waveline[i]
		endfor
		if(chi2<chi2previous)
			temp = delE
		endif
		chi2previous = chi2
	endfor
	for(i=0;i<dimsize(wavetoremove,0);i+=1)
		E = energy[i]+delE
		Ereal = Emess[x2pnt(Emess,E)]
		data[i] = wavetoremove[x2pnt(wavetoremove,Ereal)]
	endfor
	
	output[] = data[p]
	filename = "dummy_CS"
	wave delete = $filename	
	killwaves/z dummy, energy2, data, energy00, Emess
end function