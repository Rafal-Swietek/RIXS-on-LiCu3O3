#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//-----------------------------------------------------------------------------------------------------------------------------------------------
//test function using FFT, problem: not known position of frequencies creating wiggle
function testfunc(inwave)
	wave inwave
	Duplicate/O inwave, dummy
	Duplicate/c/O inwave, inwave_fft
	FFT/OUT=1/DEST=inwave_fft inwave
	Duplicate/c/O inwave_fft, wave10
	Duplicate/O wave10, reali, imagi
	redimension/r/N=(dimsize(wave10,0)) reali, imagi
	variable i,u=dimsize(wave10,0)-27,v
	Make/O/N=(dimsize(inwave,0),dimsize(wave10,0)) out
	for(v=0;v<dimsize(wave10,0);v+=1)
		wave10[] = inwave_fft[p]
		reali[] = real(wave10[p])
		imagi[] = imag(wave10[p])
		for(i=70;i<dimsize(wave10,0);i+=1)
			reali[i] = 0
			imagi[i] = 0
		endfor
		for(i=v;i<dimsize(wave10,0);i+=1)
			reali[i] = 0
			imagi[i] = 0
		endfor
		//reali[18,23] = 0
		//imagi[18,23] = 0
		wave10[] = reali[p] + cmplx(0,1)*imagi[p]
		IFFT/DEST=dummy wave10
		out[][v]= dummy[p]
	endfor
	killwaves/z dummy, reali, imagi, wave10, inwave_fft
end function


//Below function to create all permutations
//Use: test all possible FFT with cutoff frequencies
function createpermutations(inwave, left, right,index)
	wave inwave
	variable left, right, index
	Duplicate/O inwave, dummy
	variable i, dummyVar
	if(left==right)
		string filename = "permutation"+num2str(index)
		duplicate/O inwave, $filename
		wave joe = $filename
		index += 1
	else
		for(i=left;i<=right;i+=1)
			dummyVar= dummy[i]
			dummy[i] = dummy[left]
			dummy[left] = dummyVar
			
			createpermutations(dummy,left+1,right,i)
			
			dummyVar= dummy[i]
			dummy[i] = dummy[left]
			dummy[left] = dummyVar
		endfor
	endif
end

function aaaaa()
	//setdatafolder("root:125to491:D3 fit:permutations")
	make/O/N=3 rol
	rol[] = {1,2,3}
	variable index = 1
	createpermutations(rol,0,2,index)
end


//--------------------------------------------------------------------------------------
//Non used stuff saved for later	
	//E = energy2[i] + wiggle2[mod(x2pnt(wiggle2,energy2[i])+delE,dimsize(wiggle2,0))]//gemesene energy, offset on x energy axis
	//E = wiggle2[mod(x2pnt(wiggle2,energy[i]),dimsize(wiggle2,0))]+energy[i] + delE*0.001 //offset on mess energy axis
	//if(E<energy2[0])
	//	E = E + energy2[dimsize(energy2,0)-1] - energy2[0]
	//endif
	//if(E>energy[dimsize(energy,0)-1])
	//	E = E - energy2[dimsize(energy2,0)-1] + energy2[0]
	//endif
	//m = x2pnt(data,energy2[i])
	//data[m] = dummy[x2pnt(dummy,E)]
//-------------------------------------------------------------------
	//Interpolation of data
		//Variable num = 3000 //number of points for interpolation
		//Make/O/N=(num) dummy, energy2
		//interpolate2 /F=2 /N=(num) /Y = dummy waveline
		//interpolate2 /F=2 /N=(num) /Y = energy2 energy
		//interpolate2/F=2 /N=(2000) /Y = wiggle2 wig
		//setscale/i x E1, energy[dimsize(data,0)-1], dummy
//-----------------------------------------------------
	//--------------
	//Fitted model function
	//KKtransform_extrapolated(imaginalis,energy)
	//filename = NameOfWave(imaginalis)+"_KK"
	//wave wave_fit = $filename
	//wave_fit[] = 65*wave_fit[p] - 64.0003 //adjusting values to refindex aka waveline
//-----------------------------------------------------------------------------
//	redimension/c $filename
//	FitFunc[] = 1
//	Make/O/N=4 wp, w0, g
//	//Coefficients taken from RefFit
//		w0[] = {937.96, 940.94, 958.11, 961.47}
//		wp[] = {0.075757, 0.03862, 0.036529, 0.024543}
//		g[] = {0.61277, 2.0743, 0.84953, 2.6148}
//		for(i=0;i<=2;i+=1)
//			//Complex dielectric model function Lörentzian type
//			FitFunc[] -=  wp[i]^2/( energy[p]^2 - w0[i]^2 + cmplx(0,1)*g[i]*energy[p]) //fitted function as reference for chi2 test
//		endfor
//		Make/O/N=(dimsize(waveline,0)) wave_fit, Imaginalis_fit
//		wave_fit[] = real(sqrt(FitFunc[p]))
//		imaginalis_fit[] = imag(sqrt(FitFunc[p]))
//		linearize_to1(wave_fit, 925,979.9,1)
//		Duplicate/O wave_fit, wave_fitRef
//		wave_fit[] = asin(sqrt( (h*c/2/latticepar/energy[p])^2 + 2*(1-wave_fit[p]) ))
//		//wave_fit[] = acos(sqrt( wave_fit[p]^2 - (h*c/2/latticepar/energy[p])^2))
//		killwaves/z FitFunc, wp, w0, g, eps0
//-----------------------------------------------------------------------------
	//filename = "Refractive_RefFit" //function to compare in chi2 test taken from RefFit analysis of Imaginary part
	//wave FitFunc = $filename
	//Duplicate/O waveline, test
	//interpolate2/F=2 /N=(dimsize(waveline,0)) /Y = test FitFunc
	//redimension/N=(dimsize(waveline,0)) FitFunc
	//FitFunc = test
	//setscale/i x 925, 979.9, FitFunc
	//killwaves/z test
	
	function RotatebyAngle2D(theta, inwave)
//This function rotates a wave counterclockwise by an angle theta, around x-axis
	wave/c inwave
	wave theta
	Duplicate/c/O inwave, dumm
	Make/O/N=(2,2) RotationMatrix
	Variable i, j, k
	for(i=0;i<2;i+=1)
		for(j=0;j<2;j+=1)
			for(k=0;k<dimsize(inwave,0);k+=1)
				RotationMatrix[0][0] = cos(theta[k])
				RotationMatrix[1][0] = -sin(theta[k])
				RotationMatrix[0][1] = sin(theta[k])
				RotationMatrix[1][1] = cos(theta[k])
				dumm[k][i] = RotationMatrix[i][j]*inwave[k][j]
			endfor
		endfor
	endfor
	inwave[][] = dumm[p][q]
	killwaves/z dumm, RotationMatrix
	return inwave
end function

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This function extracts the imaginary part of the conductivity tensor for tetragonal symmetry 
// assuming no CMXD (Circular Magnetic X-ray Dichroism)
function OpticalCondExtract_nonMagn()
	variable i, j, k
	setdatafolder("root:")
	variable E2 = 961 //961 if LCP is taken into account
	variable num = (E2-926.8)/(0.15)+1
	string filename
	Make/O/N=(num) omega
	omega[] = 926.8 + p*0.15
	Make/O/N=(num,3) tey
	Make/O/N=(num) ang1, ang2, ang3
	variable h = 4.14e-15, c=3e8
	
		setdatafolder("root:mda_"+num2str(960)+"to"+num2str(1300))
		filename = "RefIndex_Im"+num2str(960)+"to"+num2str(1300) //LV polarisation
		//filename = "Fit_"+num2str(960)+"to"+num2str(1300)+"_abs"   //LV polarisation
		wave joe = $filename		
		tey[][0] = 2*omega[p]/h/c*joe[x2pnt(joe,omega[p])]
		//tey[][0] = joe[x2pnt(joe,omega[p])]
		setdatafolder("root:")
		filename = "ScatAng"+num2str(960)+"to"+num2str(1300)
		wave angles = $filename
		ang1[] = angles[x2pnt(joe,omega[p])]
		ang1[] *= pi/180
		
		setdatafolder("root:mda_"+num2str(1713)+"to"+num2str(2044))
		filename = "RefIndex_Im"+num2str(1713)+"to"+num2str(2044) //LH polarisation
		//filename = "Fit_"+num2str(1713)+"to"+num2str(2044)+"_abs"   //LH polarisation
		wave joe = $filename
		tey[][1] = 2*omega[p]/h/c*joe[x2pnt(joe,omega[p])]
		//tey[][1] = joe[x2pnt(joe,omega[p])]
		setdatafolder("root:")
		filename = "ScatAng"+num2str(1713)+"to"+num2str(2044)
		wave angles = $filename
		ang2[] = angles[x2pnt(angles,omega[p])]
		ang2[] *= pi/180
		
		setdatafolder("root:mda_"+num2str(575)+"to"+num2str(815))
		filename = "RefIndex_Im"+num2str(575)+"to"+num2str(815) //LCP polarisation
		//filename = "Fit_"+num2str(575)+"to"+num2str(939)+"_abs"   //LCP polarisation
		wave joe = $filename
		tey[][2] = 2*omega[p]/h/c*joe[x2pnt(joe,omega[p])]
		//tey[][2] = joe[x2pnt(joe,omega[p])]
		setdatafolder("root:")
		filename = "ScatAng"+num2str(575)+"to"+num2str(815)
		wave angles = $filename
		ang3[] = angles[x2pnt(angles,omega[p])]
		ang3[] *= pi/180
		
	filename = "ImSigma_xx"
	Duplicate/O ang1, $filename
	wave sxx = $filename
	filename = "ImSigma_zz"
	Duplicate/O ang1, $filename
	wave szz = $filename
	//filename = "ImSigma_zz_fromLCP"
	//Duplicate/O ang1, $filename
	//wave szz2 = $filename
	setscale/i x 926.8, E2, sxx, szz
	//Inverting the matrix equation µ_TEY ~ -Im( €° * Sigma * € )
		sxx[] = -tey[p][0]
		//sxx[] = tey[p][1] - 2*tey[p][2]
		szz[] = tey[p][0]*tan(ang2[p])^2 - tey[p][1]/cos(ang2[p])^2
		//szz2[] = (tey[p][0]*(1+sin(ang3[p])^2) - 2*tey[p][2])/cos(ang3[p])^2
	ModifyGraph lsize=1.2,rgb(ImSigma_zz)=(0,0,0)
	Legend/C/N=text0/J/A=RT/X=0.00/Y=0.00 "\\s("+NameOfWave(sxx)+") Sigma_xx\r\\s("+NameOfWave(szz)+") Sigma_zz"
	Label left "Optical Conductivity - imaginary part";DelayUpdate
	Label bottom "photon energy[eV]"
	 //setdatafolder("root:")		 
	 killwaves/z  omega, E1, tey, ang1, ang2, ang3
end function