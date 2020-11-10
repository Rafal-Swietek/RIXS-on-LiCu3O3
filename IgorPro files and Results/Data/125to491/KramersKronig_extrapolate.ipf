#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//This function does the Kramers-Kronig transform including extrapolation of the data at both ends
//the inputwave is the one to transform, the xdat wave, which includes the energy scaling is extrapolated with an uniform scaling
function KKtransform_extended(inputwave, xdat)

	wave inputwave, xdat 	//wave to transform and energy axis scaling
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Display inputwave vs xdat
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
	Display inputwave vs xdat
	killwaves/z outwave, inwave, right_fit_coef, left_fit_coef, xdata,W_sigma,W_fitConstants
	stop = DateTime
	print "Execution time: " + num2str(stop - start) +" seconds"
	
end function



function IKKtransform_extended(inputwave, xdat)

	wave inputwave, xdat 	//wave to transform and energy axis scaling
	String newName = NameOfWave(inputwave) + "_KK" // Compute output wave name
	Duplicate/O inputwave, $newName // Create output wave
	Display inputwave vs xdat
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
	Display inputwave vs xdat
	killwaves/z outwave, inwave, right_fit_coef, left_fit_coef, xdata, ,W_sigma,W_fitConstants
	stop = DateTime
	print "Execution time: " + num2str(stop - start) +" seconds"
	
end function

function resize_to1(inwave)
	wave inwave
	variable i, k
	k = inwave[0]
	for(i=0;i<dimsize(inwave,0);i+=1)
		inwave[i] = inwave[i]/k
	endfor

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