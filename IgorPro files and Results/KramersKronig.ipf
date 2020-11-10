#pragma rtGlobals=1		// Use modern global access method.

Function KKtransform(wavetotransform,emin,emax)

	variable emin, emax
	
	Wave wavetotransform// Reference to the input wave received as parameter
	String newName = NameOfWave(wavetotransform) + "_KK" // Compute output wave name
	Duplicate/O wavetotransform, $newName // Create output wave
	Wave wavetotransform = $newName // Create wave reference for output wave
	
	Duplicate/O wavetotransform, dummy

	dummy[][]=0
	variable i, j,k
	
	for(k=0;k<Dimsize(wavetotransform,0);k+=1)
	
	for(i=emin;i<emax;i+=1)
	for(j=emin;j<emax;j+=1)
	
	if(abs(i-j)>0.000001)
	dummy[k][i]+=wavetotransform[k][j]/(j-i)
	endif
	
	
	endfor
	endfor
	endfor

	wavetotransform[][]=dummy[p][q]
      killwaves/z maxima, dummy
End






Function KKtransform2(wavetotransform,emin,emax)

	variable emin, emax
	
	Wave wavetotransform// Reference to the input wave received as parameter
	
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

	//Duplicate/O wavetotransform, dummy

	thetawave[][]=cmplx(0,0)
	rwave[][]=cmplx(0,0)
	eps1wave[][]=0
	eps2wave[][]=0
	
	variable i, j,k
	
	for(k=0;k<Dimsize(wavetotransform,0);k+=1)
		for(i=emin;i<emax;i+=1)
			for(j=emin;j<emax;j+=1)
				if(abs(i-j)>0.000001)
					//check this formula!
					thetawave[k][i]+=ln(  sqrt(  wavetotransform[k][j]  ))/((1/j)^2-(1/i)^2)/(j^2)
				endif
			endfor
			
			thetawave[k][i]/=i
		endfor
	endfor
	
	//here you can play with the sign: -1 and +1 give totally different results
	thetawave[][]*=1
	rwave[][emin,emax]=sqrt(wavetotransform[p][q])*exp(cmplx(0,1)*thetawave[p][q])
	
	eps1wave[][]=real((1-rwave[p][q])^2 / (1+rwave[p][q])^2)
	eps2wave[][]=imag((1-rwave[p][q])^2 / (1+rwave[p][q])^2)

	
End

