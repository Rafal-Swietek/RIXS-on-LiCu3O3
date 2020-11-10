#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//	This function removes the oscilating behaviour in the 
// peak position by smoothing the data to mak a reference function and using the chi-square test
// substracting the reference oscilations with the best shift

function RemoveWiggle(wiggledummyVar, E1, dE,latticepar,p1,p2)
	wave  wiggledummyVar
	variable E1, dE, latticepar, p1, p2
	wave wavetoremove = $"Fit_"+num2str(p1)+"to"+num2str(p2)+"_x0"
	wave error = $"Fit_"+num2str(p1)+"to"+num2str(p2)+"_x0Sigma"
	linearize_to1(wiggledummyVar,917.5,979.9,0)
	
	Duplicate/O wavetoremove, data, energy, wave_fit
	string filename = Nameofwave(wavetoremove)+"_unwiggled"
	Duplicate/O wavetoremove, $filename
	wave output = $filename
	
	energy[] = E1 + 0.15*p
	setscale/p x E1, dE, energy	
	Duplicate/O wavetoremove, waveline
	//linearize_to1(waveline, 925,979.9,1)
	variable h=4.14e-15, c=3e8, i, j
		Duplicate/O waveline, angles //angles in radians
		angles[] = asin(h*c/2/energy[p]*waveline[p]/latticepar) //Seve [1]
		setscale/p x E1, dE, angles
	variable a = h*c/2/latticepar
	
	data[] = 1+0.5*(a/energy[p])^2-0.5*sin(angles[p])^2
	setscale/p x E1, dE, data
	Duplicate/O wavetoremove,FitFunc;
	Smooth 20, FitFunc
	//wave_fit[] = asin(sqrt( (h*c/2/latticepar/energy[p])^2 + 2*(1-FitFunc[p]) ))
	wave_fit[] = FitFunc[p]
	display waveline, wave_fit
	//---------------------------------------------
	variable chi2, chi2previous=1e20, dummyVar, delE, E, eps, m
	waveline[] = wavetoremove[p]
	Duplicate/O wiggledummyVar, wig
	Duplicate/O wig, wiggle2, wig2
	wiggle2[] = 0
	setscale/i x 917.5, energy[dimsize(energy,0)-1], wiggle2

	setscale/i x 917.5, 979.9, wiggle2
	setscale/p x E1, dE, data
	Duplicate/O angles, ang, angdummyVar
	variable L, tmp
	if(dimsize(wiggledummyVar,0) < dimsize(angles,0))
		Redimension/N=2000 wig, ang, wiggle2
		interpolate2/F=1/N=2000 /Y = wig wig2
		L = (979.9-E1)/(979.9-917.5)*599 + 1
		Redimension/N=(L) angdummyVar
		interpolate2/F=1 /N=(L) /Y=angdummyVar angles
		L = 2000 - L
		setscale/i x 917.5, 979.9, ang, wiggle2
		for(i=0;i<dimsize(ang,0);i+=1)
			if( i < L )
				tmp = (angles[x2pnt(angles,930)]-angles[1])/(930-E1)*(pnt2x(ang,i)-E1) + angles[1]
				//wiggle2[i] = sin(tmp)^2/cos(tmp)/a*wig[i] //wiggle in angles
				wiggle2[i] = sin(tmp)/a*wig[i] //wiggle in peak position
			else
				ang[i] = angdummyVar[i-L]
				//wiggle2[i] = sin(ang[i])^2/cos(ang[i])/a*wig[i] //wiggle in angles
				wiggle2[i] = sin(ang[i])/a*wig[i] //wiggle in peak position
			endif
		endfor
	else
		Redimension/N=(dimsize(wiggle2,0)) ang
		ang[] = 0
		L = dimsize(wiggle2,0)-dimsize(angles,0)
		Make/O/N=2 test1
		Curvefit/Q line kWCwave = test1, angles
		ang[0,L-1] = angles[0]
		for(i=L;i<dimsize(ang,0);i+=1)
			ang[i] = angles[i-L]
		endfor
		//wiggle2[] = sin(ang[p])^2/cos(ang[p])/a*wig[p] //wiggle in angles
		wiggle2[] = sin(ang[p])/a*wig[p] //wiggle in peak position
	endif
	Duplicate/O data, Ref2
	make/O/N=(dimsize(wiggle2,0)) chi
	make/O/N=(dimsize(wiggle2,0),dimsize(data,0)) wave1
	for(delE=0;delE<dimsize(wiggle2,0);delE+=1)	
		chi2 = 0
		data[] = waveline[p]
		for(i=0;i<dimsize(energy,0);i+=1)	
			data[i] -= wiggle2[mod(x2pnt(wiggle2,pnt2x(data,i))+delE,dimsize(wiggle2,0))] //for last two datasets
			//data[i] -= wiggle2[mod(i+delE,dimsize(wiggle2,0))] 
			chi2 += (data[i]-wave_fit[i])^2/error[i]
			wave1[delE][i] = data[i]
		endfor
		if(chi2<chi2previous)
			dummyVar = delE
			chi2previous = chi2
		endif
		chi[delE] =  chi2
	endfor
	variable dw = pnt2x(wig,1) - pnt2x(wig,0)
	setscale/p x 0, dw, chi, wave1
	data[] = waveline[p]
	dummyVar = 405
	for(i=0;i<dimsize(energy,0);i+=1)
			data[i] -= wiggle2[mod(x2pnt(wiggle2,pnt2x(data,i))+dummyVar,dimsize(wiggle2,0))] //for last two datasets
			//data[i] -= wiggle2[mod(i+dummyVar,dimsize(wiggle2,0))]
			//Ref2[i] = 1+0.5*(a/energy[i])^2-0.5*(h*c/2/energy[i]*data[i]/latticepar)^2 //From Angles to Refractive Index
	endfor
	output[] = data[p]
	killwaves/z  data, energy, Ref, wiggle2, angles
	killwaves/z  wig, dummy, W_sigma, W_coef, Ref2, waveline, wave_fit, test, FitFunc, wig2
end function
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This function takes peak position and width and converses them to the real and
// imaginary part of the refractive index using the formulas from [1] and [2]
Function QuantitiesFromFit(p1,p2)
	variable p1, p2
	string file = num2str(p1)+"to"+num2str(p2)

	string dataname = "root:mda_"+file
	setdatafolder(dataname)
	string filename = "ScatAng"+file
	Wave inputee = $"Fit_"+file+"_x0_unwiggled" // Wave vector from file
	Duplicate/O inputee, input
	Duplicate/O input, $filename
	Wave angles = $filename
	angles[] = 0
	Duplicate/O input, energy, position
	energy[] = 0
	variable i, h=4.14e-15, c=3e8, latticepar = 8.9e-10
	variable E1 = pnt2x(input,0), dE = pnt2x(input,1) - pnt2x(input,0)
	
	//linearize_to1(input,E1,energy[dimsize(energy,0)-1],1)
	for(i=0;i<dimsize(input,0);i+=1)
		energy[i] = E1 + i*dE
	endfor
	//Scattering Angles
		angles[] = asin(h*c/2/energy[p]*input[p]/latticepar)*180/pi //[1]
		
	//-----------------------------------------------------------------------------------------------------------------------------
	//Real part of refractive index
	
//	Uncommentize below if you want ir directly from peak position
		filename = "RefIndex_Re"+file
		Duplicate/O input , $filename
		wave Refractive = $filename
		Refractive[] = 1- 0.5*Sin(pi/180*angles[p])^2 + 0.5*(h*c)^2/(energy[p])^2/4/latticepar^2 //[2]
		//Refractive[] = sqrt( (h*c/2/energy[p]/latticepar)^2 + cos(angles[p])^2 ) //[1]- same as above
		filename = "RefIndex_Re"+file + "_unwiggled"
		Duplicate/O Refractive, $filename
		wave Refractive2 = $filename
		Refractive2[] = Refractive[p]
		//Error of real part
		filename = "RefIndex_Re"+file+"_sigma" //deviation
		wave PosSig = $"Fit_"+file+"_x0Sigma"
		Duplicate/O input , $filename
		wave Refsig = $filename
		Refsig[] = (h*c/2/latticepar/energy[p])^2*position[p]*PosSig[p]
	
		//Imaginary part of Refractive Index
		filename = "RefIndex_Im"+file
		Wave imagi = $"Fit_"+file+"_B" // width from file
		Duplicate/O imagi, $filename
		wave joe = $filename
		joe[] = (h*c/2/energy[p]/latticepar)^2*position[p]*sqrt(abs(imagi[p])) //[1]	
		//Error of imaginary part
		filename = "RefIndex_Im"+file+"_sigma" //Imaginary part deviation
		Wave Bsig = $"Fit_"+file+"_BSigma" // width from file
		Duplicate/O Bsig, $filename
		wave jochen = $filename
		jochen[] = sqrt( (Bsig[p]/2/imagi[p])^2 + (PosSig[p]/position[p])^2 )*joe[p]
	//-----------------------------------------------------------------------------------------------------------------------------
	//Real part of Scattering Amplitude
		filename = "ScatAmplitude_Re"+file
		Duplicate/O refractive, $filename
		wave f1 = $filename
		//  mass density [kg/m^3], atomic weight [u], avogadro const, electron radius [m]
		variable Ro = 8930, A = 63.55, NA = 6.022e-23, Re = 2.8179e-15 //values for Copper Cu (Z=29)
		variable AtomDen = Ro*NA/A //atomic density
		f1[] = 2*pi/AtomDen*Re*(energy[p]/h/c)^2*(1-Refractive[p]) //[6] equation 3.13, page 61	
		
		//Imaginary part of Scattering Amplitude
		filename = "ScatAmplitude_Im"+file
		Duplicate/O jochen, $filename
		wave f2 = $filename
		f2[] = 2*pi/AtomDen*Re*(energy[p]/h/c)^2*joe[p]//[6] equation 3.13, page 61
		
		//Absorption cross section
		filename = "AbsCrossSec"+num2str(p1)+"to"+num2str(p2)
		Duplicate/O jochen, $filename
		wave SigAbs = $filename
		SigAbs[] = 2*Re*h*c/energy[p]*f2[p] //[6] equation 3.28, page 64
	//-----------------------------------------------------------------------------------------------------------------------------
	//Real part of Dielectric permittivity
		filename = "Permittivity_Re"+file
		Duplicate/O refractive, $filename
		wave eps1 = $filename
		eps1[] = Refractive[p]^2 - joe[p]^2
		
		//Imaginary part of Dielectric permittivity
		filename = "Permittivity_Im"+file
		Duplicate/O jochen, $filename
		wave eps2 = $filename
		eps2[] = 2*Refractive[p]*joe[p]
	//-----------------------------------------------------------------------------------------------------------------------------
	//Reflectivity
		//	Real part-----------------------------------
		filename = "ReflectivityRe"+file
		Duplicate/O jochen, $filename
		wave refre = $filename
		refre[] = real( ( 1 - sqrt(cmplx(eps1,eps2)) )/( 1 + sqrt(cmplx(eps1,eps2)) ) )
		//	Imaginary part------------------------------------
		filename = "ReflectivityIm"+file
		Duplicate/O jochen, $filename
		wave refim = $filename
		refim[] = imag( (1 - sqrt(cmplx(eps1,eps2)))/(1 + sqrt(cmplx(eps1,eps2))) )
		//	Phase-----------------------------------
		filename = "PhaseShift"+file
		Duplicate/O refre, $filename
		wave phase = $filename
		phase[] = atan2(refim[p],refre[p]) //in radians
	//-----------------------------------------------------------------------------------------------------------------------------
	// Absorption
		filename = "LinearAbsorption"+file
		Duplicate/O joe, $filename
		wave Absorption = $filename
		Absorption[] = 2*energy[p]*joe[p]/c/h
	setdatafolder("root:")
	killwaves/z ReSigma, input, position
end

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This one does exactly the sam in tetragonal symmetry, but now we assume
// a magnetization in the z-direction
//	TEY according to Haverkord et al. [5]
function OpticalCondExtract_Magn()
	variable i, j, k
	variable E1 = 926.65, E2 = 961, dE = 0.15
	variable num = (E2-E1)/(dE)+1
	string filename
	Make/O/N=(num) omega
	omega[] = E1 + p*dE
	setdatafolder("root:OpticalConductivity")
	Make/O/N=(num,4) tey
	Make/O/N=(num) ang1, ang2, ang3
	variable h=4.14e-15, c = 3e8
		//LV polarisation
		setdatafolder("root:mda_"+num2str(960)+"to"+num2str(1300))
		//filename = "Fit_"+num2str(960)+"to"+num2str(1300)+"_abs"+"_norm" //LV polarisation
		//filename = "Xas"+num2str(960)+"to"+num2str(1300)
		filename = "XAS_"+num2str(960)+"to"+num2str(1300)+"nobackground"
		wave joe = $filename	
		
		tey[][0] = joe[x2pnt(joe,omega[p])]
		
		//---------------------------------------------------------------------------------------------------
		//LH polarisation
		setdatafolder("root:mda_"+num2str(1713)+"to"+num2str(2044))
		//filename = "Fit_"+num2str(1713)+"to"+num2str(2044)+"_abs"+"_norm" //LH polarisation
		//filename = "Xas"+num2str(1713)+"to"+num2str(2044)
		filename = "XAS_"+num2str(1713)+"to"+num2str(2044)+"nobackground"
		wave joe = $filename
		tey[][1] = joe[x2pnt(joe,omega[p])]
		
		filename = "ScatAng"+num2str(1713)+"to"+num2str(2044)
		wave angles = $filename
		ang1[] = angles[x2pnt(angles,omega[p])]
		ang1[] *= pi/180
		
		//---------------------------------------------------------------------------------------------------
		//LCP polarisation
		setdatafolder("root:mda_"+num2str(575)+"to"+num2str(815))
		//filename = "Fit_"+num2str(575)+"to"+num2str(815)+"_abs"+"_norm" //LCP polarisation
		//filename = "Xas"+num2str(575)+"to"+num2str(815)
		filename = "XAS_"+num2str(575)+"to"+num2str(815)+"nobackground"
		wave joe = $filename
		tey[][2] = joe[x2pnt(joe,omega[p])]
		
		filename = "ScatAng"+num2str(575)+"to"+num2str(815)
		wave angles = $filename
		ang2[] = angles[x2pnt(angles,omega[p])]
		ang2[] *= pi/180
		
		//---------------------------------------------------------------------------------------------------
		//RCP polarisation
		setdatafolder("root:mda_"+num2str(1324)+"to"+num2str(1689))
		//filename = "Fit_"+num2str(1324)+"to"+num2str(1689)+"_abs"+"_norm" //RCP polarisation
		//filename = "Xas"+num2str(1324)+"to"+num2str(1689)
		filename = "XAS_"+num2str(1324)+"to"+num2str(1689)+"nobackground"
		wave joe = $filename
		tey[][3] = joe[x2pnt(joe,omega[p])]
		
		filename = "ScatAng"+num2str(1324)+"to"+num2str(1689)
		wave angles = $filename
		ang3[] = angles[x2pnt(angles,omega[p])]
		ang3[] *= pi/180
		//---------------------------------------------------------------------------------------------------
		// Spin vector of Copper atoms
		Make/O/N=3 Spin
		Spin[] = {0,0,1}
		variable dummyVar = norm(Spin) 
		Spin[] = Spin[p]/dummyVar/2
		//---------------------------------------------------------------------------------------------------		
	
	setdatafolder("root:OpticalConductivity")			
	Duplicate/O ang1, sxx, szz, sxy, sxz
	Make/O/N=(4,4) M, W //conversion matrix
	M[][] = 0
	setscale/i x E1, E2, sxx, szz, sxy, sxz
	Make/c/N=(num,3,3) Sigma
	Sigma[][][] = 0
	tey[][] *= pi
	for(i=0;i<num-1;i+=1)	//0-LV,  1-LH,  2-LCP,  3-RCP
		M[0][0] = -1 //LV
		M[1][0] = -sin(ang1[i])^2 //LH
		M[1][3] = -cos(ang1[i])^2
		M[2][0] = -0.5*(1+sin(ang2[i])^2) //LCP
		M[2][1] = -2*Spin(2)*sin(ang2[i])
		M[2][2] = 2*Spin(1)*cos(ang2[i])
		M[2][3] = -0.5*cos(ang2[i])^2
		M[3][0] = -0.5*(1+sin(ang3[i])^2) //RCP
		M[3][1] = 2*Spin(2)*sin(ang3[i])
		M[3][2] = -2*Spin(1)*cos(ang3[i])
		M[3][3] = -0.5*cos(ang3[i])^2
		MatrixInverse  M //inverting the matrix equation
		wave W =$"M_inverse"
		sxx[i] = W[0][0]*tey[i][0] + W[0][1]*tey[i][1] + W[0][2]*tey[i][2] + W[0][3]*tey[i][3]
		sxy[i] = W[1][0]*tey[i][0] + W[1][1]*tey[i][1] + W[1][2]*tey[i][2] + W[1][3]*tey[i][3]
		sxz[i] = W[2][0]*tey[i][0] + W[2][1]*tey[i][1] + W[2][2]*tey[i][2] + W[2][3]*tey[i][3]
		szz[i] = W[3][0]*tey[i][0] + W[3][1]*tey[i][1] + W[3][2]*tey[i][2] + W[3][3]*tey[i][3]
	 endfor
	 
	KKtransform_extrapolated(sxx,0,0,0,10000)
	filename = NameOfWave(sxx)+"_KK"
	wave joe = $filename
	Sigma[][0][0] = joe[p] + cmplx(0,1)*sxx[p] //xx
	Sigma[][1][1] = Sigma[p][0][0] //yy = xx
	
	if(Spin[2] != 0)
		IKKtransform_extrapolated(sxy,0,0,0,10000)
		filename = NameOfWave(sxy)+"_KK"
		wave joe = $filename
		Sigma[][0][1] = 2*Spin[2]*sxy[p] + cmplx(0,1)*joe[p] //xy
		Sigma[][1][0] = -2*Spin[2]*sxy[p] - cmplx(0,1)*joe[p] //xy
	endif
	if(Spin[0] != 0 || Spin[1] != 0)
		IKKtransform_extrapolated(sxz,0,0,0,10000)
		filename = NameOfWave(sxz)+"_KK"
		wave joe = $filename 
		Sigma[][0][2] = -2*Spin[1]*sxz[p] + cmplx(0,1)*joe[p] //xz
		Sigma[][1][2] = 2*Spin[0]*sxz[p] + cmplx(0,1)*joe[p] //yz
		Sigma[][2][0] = 2*Spin[1]*sxz[p] + cmplx(0,1)*joe[p] //zx
		Sigma[][2][1] = -2*Spin[0]*sxz[p] + cmplx(0,1)*joe[p] //zy
	endif
	KKtransform_extrapolated(szz,0,0,0,10000)
	filename = NameOfWave(szz)+"_KK"
	wave joe = $filename
	Sigma[][2][2] = joe[p] + cmplx(0,1)*szz[p] //zz 
	
	CreateTensorElements(Sigma,"OptCond",E1,E2)
	DieTensFromOptic(Sigma,E1,E2,dE) //for now it creates the electric susceptibility
	
	filename = "OpticalCondFromTEY"
	Duplicate/c/O Sigma, $filename
	wave/c jochen = $filename
	jochen[][][] = Sigma[p][q][r]
	 killwaves/z  omega, tey, ang1, ang2, ang3, M, W, sxx, szz, sxy, sxz, W, Spin, Sigma
	 setdatafolder("root:")
end function
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//	This function recreates the Dielectric tensor from the dielectric function for given polarizations

function DielectricTensor()
	variable i, j, k, a, b
	setdatafolder("root:DielectricTensor")
	variable E1 = 926.65, E2 = 961, dE = 0.15
	variable num = (E2-E1)/(dE)+1
	variable/c dummyVar
	string filename
	Make/O/N=(num) omega
	omega[] = E1 + p*dE
	Make/c/O/N=(num,4) Epsilon
	Make/O/N=(num) ang1, ang2, ang3
	setscale/i x E1, E2, ang1, ang2, ang3
		Epsilon[][] = 0
		
		//LV polarisation
		setdatafolder("root:mda_"+num2str(960)+"to"+num2str(1300))
		filename = "RefIndex_Re"+num2str(960)+"to"+num2str(1300)+"_unwiggled"  
		wave joe = $filename
		filename = "RefIndex_Im"+num2str(960)+"to"+num2str(1300) 
		wave jochen = $filename
			
		Epsilon[][0] = ( joe[x2pnt(joe,omega[p])] +cmplx(0,1)* jochen[x2pnt(jochen,omega[p])] )^2
		
		//--------------------------------------------------------------------------------------------------- 
		//LH polarisation
		setdatafolder("root:mda_"+num2str(1713)+"to"+num2str(2044))
		filename = "ScatAng"+num2str(1713)+"to"+num2str(2044) 
		wave angles = $filename
		ang1[] = angles[x2pnt(angles,omega[p])]
		ang1[] *= pi/180
		
		filename = "RefIndex_Re"+num2str(1713)+"to"+num2str(2044)+"_unwiggled"  
		wave joe = $filename
		filename = "RefIndex_Im"+num2str(1713)+"to"+num2str(2044) 
		wave jochen = $filename
			
		Epsilon[][1] = ( joe[x2pnt(joe,omega[p])] +cmplx(0,1)* jochen[x2pnt(jochen,omega[p])] )^2 
		//---------------------------------------------------------------------------------------------------
		//LCP polarisation
		setdatafolder("root:mda_"+num2str(575)+"to"+num2str(815))
		filename = "ScatAng"+num2str(575)+"to"+num2str(815) 
		wave angles = $filename
		ang2[] = angles[x2pnt(angles,omega[p])]
		ang2[] *= pi/180
		
		filename = "RefIndex_Re"+num2str(575)+"to"+num2str(815)+"_unwiggled"  
		wave joe = $filename
		filename = "RefIndex_Im"+num2str(575)+"to"+num2str(815) 
		wave jochen = $filename
		
		Epsilon[][2] = ( joe[x2pnt(joe,omega[p])] +cmplx(0,1)* jochen[x2pnt(jochen,omega[p])] )^2
		//---------------------------------------------------------------------------------------------------
		//RCP polarisation
		setdatafolder("root:mda_"+num2str(1324)+"to"+num2str(1689))
		filename = "ScatAng"+num2str(1324)+"to"+num2str(1689) 
		wave angles = $filename
		ang3[] = angles[x2pnt(angles,omega[p])]
		ang3[] *= pi/180
		
		filename = "RefIndex_Re"+num2str(1324)+"to"+num2str(1689)+"_unwiggled" 
		wave joe = $filename
		filename = "RefIndex_Im"+num2str(1324)+"to"+num2str(1689)  
		wave jochen = $filename
	
		Epsilon[][3] = ( joe[x2pnt(joe,omega[p])] +cmplx(0,1)* jochen[x2pnt(jochen,omega[p])] )^2
		
	setdatafolder("root:DielectricTensor")			
	Make/c/O/N=(3,3) M //conversion matrix
	M[][] = 0
	Make/c/O/N=(num,3,3) EpsTens
	EpsTens[][][] = cmplx(0,0)
	for(i=0;i<num-1;i+=1)
//		M[0][0] = 1
//		M[1][0] = sin(ang1[i])^2
//		M[1][2] = cos(ang1[i])^2
//		M[2][0] = 0.5*( 1 + sin(ang2[i])^2)
//		M[2][1] = cmplx(0,1)*sin(ang2[i])
//		M[2][1] = cmplx(0,1)*cos(ang2[i])
//		M[2][3] = 0.5*cos(ang2[i])^2
//		M[3][0] = 0.5*( 1 + sin(ang3[i])^2)
//		M[3][1] = -cmplx(0,1)*sin(ang3[i])
//		M[3][1] = -cmplx(0,1)*cos(ang3[i])
//		M[3][3] = 0.5*cos(ang3[i])^2
//		MatrixInverse M
//		wave/c W = $"M_inverse"
//		EpsTens[i][0][0] = W[0][0]*Epsilon[i][0] + W[0][1]*Epsilon[i][1] + W[0][2]*Epsilon[i][2] + W[0][3]*Epsilon[i][3]  //xx
//		EpsTens[i][0][1] = W[1][0]*Epsilon[i][0] + W[1][1]*Epsilon[i][1] + W[1][2]*Epsilon[i][2] + W[1][3]*Epsilon[i][3]//xy
//		EpsTens[i][0][2] = W[2][0]*Epsilon[i][0] + W[2][1]*Epsilon[i][1] + W[2][2]*Epsilon[i][2] + W[2][3]*Epsilon[i][3] //xz
//		EpsTens[i][2][2] = W[3][0]*Epsilon[i][0] + W[3][1]*Epsilon[i][1] + W[3][2]*Epsilon[i][2] + W[3][3]*Epsilon[i][3] //zz
		M[0][0] = 1
		M[1][0] = sin(ang2[i])^2 
		M[1][2] = cos(ang2[i])^2
		M[2][0] = 0.5*( 1+sin(ang2[i])^2 )
		M[2][1] = cmplx(0,1)*sin(ang2[i])
		M[2][2] = 0.5*cos(ang2[i])^2
		MatrixInverse M
		wave/c W = $"M_inverse"
		EpsTens[i][0][0] = W[0][0]*Epsilon[i][0] + W[0][1]*Epsilon[i][1] + W[0][2]*Epsilon[i][2] //xx
		EpsTens[i][0][1] = W[1][0]*Epsilon[i][0] + W[1][1]*Epsilon[i][1] + W[1][2]*Epsilon[i][2] //xy
		EpsTens[i][2][2] = W[2][0]*Epsilon[i][0] + W[2][1]*Epsilon[i][1] + W[2][2]*Epsilon[i][2] //zz
	endfor
	EpsTens[][1][0] = -EpsTens[p][0][1] //yx
	//EpsTens[][2][0] = -EpsTens[p][0][2] //zx
	//EpsTens[][1][2] = EpsTens[p][0][2] //yz
	//EpsTens[][2][1] = -EpsTens[p][1][2] //zy
	EpsTens[][1][1] = EpsTens[p][0][0] //yy
	
	//Symmetric part only - assumption
//		EpsTens[][][] = 0
//		EpsTens[][0][0] = Epsilon[p][0]
//		EpsTens[][1][1] = EpsTens[p][0][0]
//		EpsTens[][2][2] = ( Epsilon[p][0]^2*Epsilon[p][1]^2*sin(ang1[p])^2 )/( Epsilon[p][0]^2 + Epsilon[p][1]^2*cos(ang1[p])^2 )
	//-----------------------------
	killwaves/z  omega, Epsilon, ang1, ang2, ang3,M,W, M_inverse 
	killwaves/z EpsdummyVar
	 
	CreateTensorElements(EpsTens,"Eps",E1,E2)
	OptTensFromDielectric(EpsTens,E1, E2, dE)
	
	 setdatafolder("root:")
end function

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//Creates the elements of a given Tensor with user name
function CreateTensorElements(tensor,filename,E1,E2)
	wave tensor
	variable E1, E2
	string filename
	string H = "xyz", name
	variable i, j, num = dimsize(tensor,0)
	make/O/N=(num) test
	for(i=0;i<3;i+=1)
		for(j=0;j<3;j+=1)
			name = filename+"Re_" + (H[i]+H[j])
			Duplicate/O test, $name
			wave joe = $name
			joe[] = real(tensor[p][i][j])
			setscale/i x, E1, E2, joe
			name = filename+"Im_" + (H[i]+H[j])
			Duplicate/O test, $name
			wave joe = $name
			joe[] = imag(tensor[p][i][j])
			setscale/i x, E1, E2, joe
		endfor
	endfor
	killwaves/z test
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	Extracts the optical conductivity from the Dielectric Tensor according to [4]: equation 3.7
function OptTensFromDielectric(Epsilon, E1, E2, dE)
	wave Epsilon
	variable E1, E2, dE
	variable h=4.14e-15, e0 = 8.85e-12
	variable num = dimsize(Epsilon,0), i, a, b, dummyVar, kmod
	setdatafolder("root:OpticalConductivity")
	Make/O/N=(num) angles, energy
	energy[] = E1 + dE*p	
	
	Duplicate/C/O Epsilon, Sigma
	for(i=0;i<num;i+=1)
		for(a=0;a<3;a+=1)
			for(b=0;b<3;b+=1)
				if(a==b)
					Sigma[i][a][b] = -( Epsilon[i][a][b] - 1 )*cmplx(0,1)*2*pi*e0*energy[i]/h
				else
					Sigma[i][a][b] = -Epsilon[i][a][b]*cmplx(0,1)*2*pi*e0*energy[i]/h
				endif
			endfor
		endfor
	endfor
	CreateTensorElements(Sigma,"OptCond",E1,E2)
	string filename = "OpticalCondFromDielectric"
	Duplicate/c/O Sigma, $filename
	wave/c joe = $filename
	joe[][][] = Sigma[p][q][r]
	killwaves/z k, angles, Eps, G, EpsdummyVar, Sigma, energy
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	Extracts the the Dielectric Tensor from optical conductivity according to [4]: equation 3.7
function DieTensFromOptic(Sigma, E1, E2, dE)
	wave Sigma
	variable E1, E2, dE
	variable h=4.14e-15, e0 = 8.85e-12, c = 3e8
	variable num = dimsize(Sigma,0), i, a, b
	setdatafolder("root:DielectricTensor")
	Make/O/N=(num) omega
	omega[] = (E1 + dE*p)/h*2*pi //frequancy
	
	Duplicate/C/O Sigma, Epsilon
	Epsilon[][][] = 0
	Make/N=(dimsize(Sigma,0),3) k
	wave angles = $"angles"
	k[][] = 0
	k[][1] = omega[p]/c*cos(angles[p])
	k[][2] = -omega[p]/c*sin(angles[p])
	variable k2, dummyVar
	make/o/c/N=(3,3) EpsdummyVar
	for(i=0;i<num;i+=1)
		for(a=0;a<3;a+=1)
			for(b=0;b<3;b+=1) 
				//k2 = omega[i]/c
				//dummyVar = 1/(omega[i]^2 - c^2*k2)
				if(a==b)
					EpsdummyVar[a][b] = 1 - Sigma[i][a][b]*cmplx(0,1)/omega[p]/e0 //effective dielectric tensor [4]
					//Epsilon[i][a][b] = 1 + dummyVar*(omega[p]^2 - c^2*k[i][a]*k[i][b])*cmplx(0,1)*Sigma[i][a][b]/e0/omega[p] 
				else
					EpsdummyVar[a][b] = -Sigma[i][a][b]*cmplx(0,1)/omega[p]/e0 //effective dielectric tensor [4]
					//Epsilon[i][a][b] = dummyVar*(omega[p]^2 - c^2*k[i][a]*k[i][b])*cmplx(0,1)*Sigma[i][a][b]/e0/omega[p]
				endif
			endfor
		endfor
		//MatrixInverse/O EpsdummyVar // with inversion (general solution) change minus sign in loop to '+'
		Epsilon[i][][] = EpsdummyVar[q][r]
	endfor
	CreateTensorElements(Epsilon,"Eps",E1,E2)
	string filename = "DielectricFromOptical"
	Duplicate/c/O Epsilon, $filename
	wave/c joe = $filename
	joe[][][] = Epsilon[p][q][r]
	killwaves/z k, omega, Eps, G, Epsilon, k, EpsdummyVar
	setdatafolder("root:OpticalConductivity")	
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	This function substracts the background of the continuum excitation using a shirley algorithm
function SubstractXASBackround(p1, p2, Eend, Em1, Em2,Em3)
	variable p1, p2, Eend, Em1, Em2, Em3
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	filename = "Fit_"+file+"_abs"
	wave joe = $filename
	Duplicate/O joe, dummy, energy
	variable E1 = pnt2x(joe,0), dE = pnt2x(joe,1) - pnt2x(joe,0)
	energy[] = E1 + dE*p
	Make/O/N=2 test1, test2
	Curvefit/Q=1 line kwCWave=test1 dummy[0,x2pnt(energy,932)] /X=energy
	
	dummy[] -= test1[0] + test1[1]*energy[p]
	//dummy[] -= joe[x2pnt(joe,930.5)]
	Duplicate/O dummy, dummy2
	
	variable y1 =dummy[x2pnt(dummy,932)]
	variable h = dummy[x2pnt(dummy,945)] - dummy[x2pnt(dummy,935)]
	variable h2 = dummy[x2pnt(dummy,975)] - dummy[x2pnt(dummy,955)]
	Duplicate/O joe, back, back2
	back2[] = y1 + 0.75*h/2 + 3/4*h/pi*atan((energy[p] - 938.05)) + 0.25*h/2 + 1/4*h/pi*atan((energy[p] - 941.05))
	back2[] += 0.75*h2/2 + 3/4*h2/pi*atan((energy[p] - 958.15)) + 0.25*h2/2 + 1.2/4*h2/pi*atan((energy[p] - 961))
	
	Duplicate/O dummy, dummy1
	dummy1[] = dummy[p] - back2[p]
	
	dummy[] = dummy2[dimsize(joe,0)-p-1]
	
	CreateShirleyBackground(dummy,back,dE,E1,Eend,Em1,Em2,Em3)
	
	dummy[] = dummy2[p] - back[p]
	filename = "XAS_"+file+"nobackground"
	Duplicate/O dummy, $filename
	wave joe = $filename
	joe[] = dummy[p]
	display joe, back, dummy2
	ModifyGraph lstyle(back)=3,lsize(back)=1.2,lsize(dummy2)=1.1;DelayUpdate
	ModifyGraph rgb(back)=(0,0,0)
	ModifyGraph rgb(dummy2)=(16384,28160,65280)
	
	setdatafolder("root:")
	killwaves/z test1, test2, energy, W_sigma, test, back2, dummy1, dummy
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

function RecreateRIXS(p1, p2)
	variable p1, p2

	variable h = 4.14e-15/2/pi, c = 3e8, eC = 1.602e-19
	string foldername = "root:OpticalConductivity"
	setdatafolder(foldername)
	string filename = "OpticalCondFromTEY"
	wave/c sigma = $filename
	variable dim2 = dimsize(sigma,0)
	variable E1 = 926.65, E2 = 961, dE = 0.15	
	Make/O/N=(dim2) omega
	omega[] = E1 + dE*p
	setscale/p x E1, dE, omega
	variable const = h^2*c/eC/2/pi
	
	Make/C/N=(3) ein, eout
	ein[] = 0;  eout[] = 0 //polarizatioon vecotrs for given angles
	
	QuantitiesFromFit(p1,p2)
	foldername = "root:mda_"+num2str(p1)+"to" + num2str(p2)
	setdatafolder(foldername)
	filename = "ScatAng"+num2str(p1)+"to" + num2str(p2)
	wave alfa = $filename //in degree's
	filename = "PhaseShift"+num2str(p1)+"to" + num2str(p2)
	wave ksi = $filename
	filename = "D3_"+num2str(p1)+"to" + num2str(p2)
	wave rixsD3 = $filename
	
	foldername = "root:RIXSrecreation"
	setdatafolder(foldername)
	
	variable dim1 = 200//dimsize(rixsD3,0)
	Make/O/N=(dim1, dim2) dummy
	dummy[][] = 0
	variable i, j, k
	variable dummyVar, a, da = 0.001
	variable/c tmp
	Make/O/N=(dim1) theta
	variable alfa0 = 15, alfa1 = 75
	//variable alfa0 = wavemin(alfa), alfa1 = wavemax(alfa)
	theta[] = ( alfa0 + (alfa1-alfa0)/(dim1-1)*p)*pi/180 //in radians
	//theta[] *= 2
	for(k=0;k<dim2;k+=1)
		for(i=0;i<dim1;i+=1)
			switchPolarization(ein, theta[i], p1) //incident polarization;  in dataanalysis.ipf
			switchPolarization(eout, pi - theta[i], 575) //scattered polarization
			//eout[0] = cos(pi/2)
			//eout[1] = sin(pi/2)*eout[1]
			//eout[2] = sin(pi/2)*eout[2]
			tmp   = conj(eout[0])*sigma[k][0][0]*ein[0] + conj(eout[0])*sigma[k][0][1]*ein[1] + conj(eout[0])*sigma[k][0][2]*ein[2]
			tmp += conj(eout[1])*sigma[k][1][0]*ein[0] + conj(eout[1])*sigma[k][1][1]*ein[1] + conj(eout[1])*sigma[k][1][2]*ein[2]
			tmp += conj(eout[2])*sigma[k][2][0]*ein[0] + conj(eout[2])*sigma[k][2][1]*ein[1] + conj(eout[2])*sigma[k][2][2]*ein[2]
			//tmp   = conj(eout[0])*sigma[k][0][0]*ein[0] + conj(eout[1])*sigma[k][1][1]*ein[1] + conj(eout[2])*sigma[k][2][2]*ein[2]
			//tmp   = conj(eout[0])*sigma[k][0][1]*ein[1] + conj(eout[0])*sigma[k][0][2]*ein[2] + conj(eout[1])*sigma[k][1][0]*ein[0]
			//tmp += conj(eout[1])*sigma[k][1][2]*ein[2] + conj(eout[2])*sigma[k][2][0]*ein[0] + conj(eout[2])*sigma[k][2][1]*ein[1] 
			dummy[i][k] = magsqr(tmp)*omega[k]^2
		endfor
	endfor
	filename = "RIXS_"+num2str(p1)+"to"+num2str(p2)
	Duplicate/O dummy, $filename
	wave RIXS = $filename
	RIXS[][] = dummy[p][q]
	//RIXS[][] = rixsD3[p][q] - dummy[p][q]
	setscale/p y E1, dE, RIXS
	setscale/i x alfa0, alfa1, RIXS
	
	//NewImageTool5("root:RIXSrecreation:"+filename)
	
	//RIXS recreation along bragg peak--------------------
		RIXSonpeak(p1,p2)
		foldername = "root:mda_"+num2str(p1)+"to" + num2str(p2)
		setdatafolder(foldername)
		filename = "RIXS_onpeak_"+num2str(p1)+"to" + num2str(p2)
		wave rixsD3 = $filename
		foldername = "root:RIXSrecreation"
		setdatafolder(foldername)
		filename = "RIXSalongpeak_"+num2str(p1)+"to" + num2str(p2)
		Duplicate/O rixsD3, $filename
		wave RIXSpeak = $filename
		RIXSpeak[] = 0
		Redimension/N=(dim2) RIXSpeak
		setscale/p x E1, dE, RIXSpeak
		for(k=0;k<dim2;k+=1)
				switchPolarization(ein, alfa[k]*pi/180, p1) //incident polarization;  in dataanalysis.ipf
				eout[0] = cos(ksi[k])
				eout[1] = sin(ksi[k])*sin(alfa[k]*pi/180)
				eout[2] = sin(ksi[k])*cos(alfa[k]*pi/180)
				tmp =   conj(eout[0])*sigma[k][0][0]*ein[0] + conj(eout[0])*sigma[k][0][1]*ein[1] + conj(eout[0])*sigma[k][0][2]*ein[2]
				tmp += conj(eout[1])*sigma[k][1][0]*ein[0] + conj(eout[1])*sigma[k][1][1]*ein[1] + conj(eout[1])*sigma[k][1][2]*ein[2]
				tmp += conj(eout[2])*sigma[k][2][0]*ein[0] + conj(eout[2])*sigma[k][2][1]*ein[1] + conj(eout[2])*sigma[k][2][2]*ein[2]
				RIXSpeak[k] = magsqr(tmp)/omega[k]^2
		endfor
		Duplicate/O RIXSpeak, dummy3
		Duplicate/O rixsD3, dummy2
		//variable maxx = wavemax(dummy2)
		//dummy2[] = dummy2[p]/maxx*wavemax(dummy3)
//		display dummy2, dummy3
//		ModifyGraph lsize=1.2,rgb(dummy3)=(0,0,0)
//		ModifyGraph nticks(left)=0
//		Legend/C/N=text0/J/A=RT/X=0.00/Y=0.0 "\\s(dummy2) directly from RIXS measurment \r\\s(dummy3) from conductivity tensor"
//		Label left "RIXS intensity (a.u.)";DelayUpdate
//		Label bottom "photon energy[eV]"
	
	killwaves/z ein, eout, dummy, dummy2, dummy3, omega, theta
	//setdatafolder("root:")
end function

function probability(a,fi) //probability of given polarization with a and fi
	variable a, fi
	return sqrt(1-a^2)*exp(-a*fi) + sqrt(1+a^2)*exp(a*fi)
end function 
//	Integrated RIXS over outoput polarzation
//	variable/c v1, v2, v3, int1, int2, int3, int4
//	int1 = Int2D(0,1,0,pi,1)
//	int2 = Int2D(0,1,0,pi,2)
//	int3 = Int2D(0,1,0,pi,3)
//	int4 = Int2D(0,1,0,pi,4)
//			v1 = sigma[k][0][0]*ein[0] + sigma[k][0][1]*ein[1] + sigma[k][0][2]*ein[2]
//			v2 = sigma[k][1][0]*ein[0] + sigma[k][1][1]*ein[1] + sigma[k][1][2]*ein[2]
//			v3 = sigma[k][2][0]*ein[0] + sigma[k][2][1]*ein[1] + sigma[k][2][2]*ein[2] 	
//			
//			dummyVar = 0
//			dummy[i][k]= 0
//			dummy[i][k] += int1*magsqr(v1)
//			dummy[i][k] += int2*sin(theta[i])^2*magsqr(v2)
//			dummy[i][k] += int2*cos(theta[i])^2*magsqr(v3)
//			dummy[i][k] -= int2*sin(theta[i])*cos(theta[i])*2*( real(v2)*real(v3) + imag(v2)*imag(v3) )
//			dummy[i][k] += 2*( real(v1)*real(v2) + imag(v1)*imag(v2) )*sin(theta[i])*int3 + 2*( real(v1)*imag(v2) - imag(v1)*real(v2) )*sin(theta[i])*int4
//			dummy[i][k] -= 2*( real(v1)*real(v3) + imag(v1)*imag(v3) )*cos(theta[i])*int3 + 2*( real(v1)*imag(v3) - imag(v1)*real(v3) )*cos(theta[i])*int4 

function fun1(a, fi, var)
	variable a, fi, var
	switch(var)
		case 1:
			return a^2
			break;
		case 2:
			return 1-a^2
			break;
		case 3:
			return a*sqrt(1-a^2)*cos(fi)
			break;
		case 4:
			return a*sqrt(1-a^2)*sin(fi)
			break;
	endswitch
end function

function Int2D(xmin, xmax, ymin, ymax,var)
	variable xmin, xmax, ymin, ymax, var
	variable xV, yV, dy = (ymax-ymin)/100, dx = 0.01//(xmax-xmin)/100
	variable int = 0
	for(xV=xmin;xV<=xmax;xV+=dx)
		for(yV=ymin;yV<=ymax;yV+=dy)
			int += fun1(xV,yV,var)*probability(xV,yV)*dy*dx
		endfor
	endfor
	return int;
end function
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
function BackRIXS(p1, p2) //RIXS background
	variable p1, p2
	string file = num2str(p1)+"to"+num2str(p2)
	setdatafolder("root:mda_"+file)
	string filename = "D3_"+file
	wave rixs = $filename
	Duplicate/O rixs, rixs_back
	variable dim1 = dimsize(rixs,0), dim2=dimsize(rixs,1), i
	variable x1 = pnt2x(rixs,0), dx = pnt2x(rixs,1) - pnt2x(rixs,0)
	Make/O/N=(dim1) lor, back, xdat
	Make/N=2 test1
	xdat[] = x1 + p*dx
	setscale/p x x1, dx, lor
	for(i=0;i<dim2;i+=1)
		lor[] = rixs[p][i]
		Curvefit/Q=1 line kwCWave=test1 lor /X=xdat
		back[] = test1[0] + test1[1]*xdat[p]
		lor[] -= back[p]
		rixs_back[][i] = lor[p]
	endfor
	killwaves/z, lor, back, xdat
end function

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
// Refrences for the formulas 
//Bibliography
//	[1] Determination of the Anomalous Scattering Factors in the Soft-X-ray Range using Diffraction
//from a Multilayer,  L.Seve
//	[2] Real-part EXAFS from multilayer Bragg reflections, U.Staub
//	[3] X-ray interactions, B.L. Henke
//	[4] Wavevector-dependent optical properties from wavevector-independent proper conductivity tensor, J. Kortus
//	[5] Theory of Resonant Inelastic X-Ray Scattering by Collective Magnetic Excitations, M.W. Haverkort
//	[6] Soft X-rays and extreme ultraviolet radiation, David Attwood




//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//	Evaluates the absorption cross section from the conductivity tensor extracted from the dielectricity tensor
// and evaluates the difference assossiated with magnon dispersion and higher order transitions
function RestIntensity(p1, p2)
	variable p1, p2
	
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	filename = "XAS_"+file+"nobackground"
	wave Xas = $filename
	filename = "ScatAng"+file
	wave alfa = $filename
	
	filename = "root:RestIntensity"
	setdatafolder(filename)
	filename ="RestXAS"+file
	Duplicate/O Xas, $filename
	wave joe = $filename
	filename ="XAS"+file
	Duplicate/O Xas, $filename
	wave jochen = $filename
	jochen[] = Xas[p]
	
	filename = "root:OpticalConductivity"
	setdatafolder(filename)
	filename = "OpticalCondFromDielectric"
	wave/c SigmaDiel = $filename
	
	variable num = dimsize(SigmaDiel,0)
	Make/O/N=(num) omega
	Variable E1 = 926.65, E2 = 961, dE = (E2 - E1)/(num-1)
	omega[] = E1 + p*dE
	Redimension/N=(num) joe
	setscale/i x ,E1, E2, joe
	Make/c/O/N=3 polar //wave polarization
	
	filename = "root:RestIntensity"
	setdatafolder(filename)
	variable i,j, k
	Make/c/O/N=(num) dummyVar
	dummyVar[] = 0
	for(k=0;k<num;k+=1)
		switchPolarization(polar, alfa[k]*pi/180, p1)
		for(i=0;i<3;i+=1)
			for(j=0;j<3;j+=1)
				dummyVar[k] += conj(polar[i])*SigmaDiel[k][i][j]*polar[j]
			endfor
		endfor
	endfor
	joe[] =  Xas[x2pnt(Xas,omega[p])] + imag(dummyVar[p]) 
	
	setdatafolder("root:")
	killwaves/z polar, ang, omega, dummyVar
	
end function





//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
// 	For dielectric tensor from refractive Index via change of basis

//	for(i=0;i<num-1;i+=1)	
//		M[0][0] = 1
//		M[1][1] = cmplx(0,1)*sin(ang1[i])/sqrt(2)
//		M[1][2] = cmplx(0,1)*cos(ang1[i])/sqrt(2)
//		M[2][1] = -cmplx(0,1)*sin(ang2[i])/sqrt(2) + cos(ang1[i]-ang2[i])*sin(ang1[i])/sqrt(2)
//		M[2][2] = -cmplx(0,1)*cos(ang2[i])/sqrt(2) + cos(ang1[i]-ang2[i])*cos(ang1[i])/sqrt(2)
//		MatrixInverse  M //inverting the matrix equation
//		wave/c W =$"M_inverse"
//		M[][] = W[p][q]
//		MatrixTranspose W
//		for(j=0;j<3;j+=1) // Matrix multiplication;  out_jk = M_ja * eps_ab * W_bk
//			for(k=0;k<3;k+=1) //Index i goes over the energy
//				for(a=0;a<3;a+=1)
//					for(b=0;b<3;b+=1)
//						EpsOut[i][j][k] += -omega[i]/4/pi*M[j][a]*Epsilon[i][a][b]*W[b][k]
//					endfor
//				endfor
//			endfor
//		endfor
//		
//	endfor
//	
//	 	for(j=0;j<3;j+=1) //Creating Output data with proper name
//	 		for(k=0;k<3;k+=1) 
//	 			//Real part
//	 			filename = "EpsRe_"+H[j]+H[k]
//				Duplicate/O ang1, $filename
//				wave joe = $filename
//				joe[] = real(EpsOut[p][j][k])
//				setscale/i x 926.65, 960.85, joe
//				
//				//Imaginary part
//				filename = "EpsIm_"+H[j]+H[k]
//				Duplicate/O ang1, $filename
//				wave joe = $filename
//				joe[] = imag(EpsOut[p][j][k]) 
//				setscale/i x 926.65, 960.85, joe
//			endfor
//		endfor



function DeleteNan(p1,p2)
	variable p1, p2
	string file = num2str(p1)+"to"+num2str(p2)
	string filename = "root:mda_"+file
	setdatafolder(filename)
	filename = "Fit_"+file+"_"
	variable pnt1 = 29
	variable pnt2 = 523
	wave joe = $filename+"y0"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"yoSigma"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"x0"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"x0Sigma"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"A"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"ASigma"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"B"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	wave joe = $filename+"BSigma"
	joe[pnt1] = joe[pnt1+1]
	joe[pnt2] = joe[pnt2-1]
	
end function