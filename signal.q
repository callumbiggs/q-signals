// External Function to be called to perform the FFT
// 
spectral:{[vec;fs]
	// Get the length of the series, ensure it is a power of 2
	// If not, we should abort operation and output a reason
	nr:count first vec;
	if[not 0f  =(2 xlog nr) mod 1;0N "The input vector is not a power of 1";:()];

	// Window the functions, using the Hanning window
	wndFactor:{[n] {[nr;x](sin((x*(acos -1))%(nr-1))) xexp 2} [n;] each til n};
	vec*:(wndFactor[nr];wndFactor[nr]);
	
	// Compute the fft and get the absolute value of the result
	fftRes:.signal.fftrad2[vec];
	magnitude:.signal.mag fftRes;

	// Scale an x axis
	xAxis:{[Ns;fs] (neg fs*0.5)+(fs%(Ns-1))*til Ns}[nr;fs];
	
	// Table these results
	([]freq:xAxis;real:fftRes 0;im:fftRes 1;magnitude:magnitude%nr)
	};


\d .signal
// Global constants
PI:acos -1; / pi ;
BR:2 sv reverse 2 vs til 256; / Reversed bits in bytes 00-FF
P2:1,prds[7#2],1; / Powers of 2 with outer guard
P256:prds 7#256; / Powers of 256

// Complex Functions
mult:{[vec1;vec2]
	// Performs the dot product of two vectors, 
	realOut:((vec1 0) * (vec2 0)) - ((vec1 1) * (vec2 1));
	imagOut:((vec1 1) * (vec2 0)) + ((vec1 0) * (vec2 1));
	(realOut;imagOut)};

division:{[vec1;vec2]
	// We can use the 
	denom:1%((vec1 1) xexp 2) + ((vec2 1) xexp 2);
	mul:mult[(vec1 0);(vec1 1);(vec2 0);neg (vec2 1)];
	realOut:(mul 0)*denom; 
	imagOut:(mul 1)*denom;
	(realOut;imagOut)};

conj:{[vec](vec 0;neg vec 1)};

mag:{[vec]
	sqrvec:(vec xexp 2);
	sqrt (sqrvec 0)+(sqrvec 1)};



// Bit Reversal
bitreversal:{[indices] 
	// Applies a bitwise reversal for a list of indices
    ct:ceiling 2 xlog last indices; / Number of significant bits (assumes sorted)
    bt:BR 256 vs indices; / Breakup into bytes and reverse bits
    dv:P2[8-ct mod 8]; / Divisor to shift leading bits down
    bt[0]:bt[0] div dv; / Shift leading bits
    sum bt*count[bt]#1,P256 div dv / Reassemble bytes back to integer
    };

// Radix 2 DIT 
fftrad2:{[vec]
	// This performs a FFT for samples of power 2 size
	// First, define some constants
	n:count vec 0;
	n2:n div 2;
	indexOrig:til n;

	// Twiddle Factors - precomputed the complex values over the discrete angles, using Euler formula
	angles:{[n;x]2*PI*x%n}[n;] til n2;
	.const.twiddle:(cos angles;neg sin angles);

	// Bit-reversed the vector and define it into a namespace so lambdas can access it 
	ind:bitreversal[indexOrig];
	.res.vec:`float$vec . (0 1;ind);

	// Precomputing the indices required to immplement each temporal phase of the DIT
		// Number of signals
		signalcount:`int$ {2 xexp x} 1+ til ( `int$(2 xlog n));
		// Number of points in each signal
		signalpoints:reverse signalcount;
		// Define an initial count, this can be used as a backbone to get the even and odd indices
		// by adjustment
		initial:{[n2;x]2*x xbar til n2}[n2;] peach (signalcount div 2);
		evens:{[n2;x]x + (n2)#til n2 div count distinct x}[n2;] peach initial;
		odds:evens (+)' (signalcount div 2);
		twiddleIndices:{[n;n2;x](n2)#(0.5*x) * til (n div (x))}[n;n2;] peach signalpoints;


	// Butterfly Implementation
	bflyComp:{[bfInd] 
		tmp:mult[.res.vec[;bfInd 1];.const.twiddle[;bfInd 2]];
		.[`.res.vec;(0 1;bfInd 1);:;.res.vec[;bfInd 0]-tmp];
		.[`.res.vec;(0 1;bfInd 0);+;tmp];
		.res.vec};

	bflyComp each (flip `int$(evens;odds;twiddleIndices));
	.res.vec
	};


movAvg:{[list;N]
	$[0=N mod 2;
		(floor N%2) rotate 2 mavg (N mavg list);
		(floor N%2) rotate N mavg list]};

movDev:{[list;N]
	$[0=N mod 2;
		(floor N%2) rotate 2 mdev (N mdev list);
		(floor N%2) rotate N mdev list]};



/------ The following are not used in the whitepaper but are functional. 
// A more fully fledged application would use both a Bluestein and Radix-2 Transformation to remove the requirement that 
// processed signals must be in powers of 2. 
// This will be developed further in the future

// This generic FFT will sort out if the a FFT requires a Bluestein or Radix 2 Transformation
fft:{[vec]
	nr:count vec 1;
	// Distributing to either radix-2 or Bluesteins Z-chrip
	$[0f  =(2 xlog nr) mod 1;res:.fft.fftrad2[vec];
		res:.fft.fftBluestein[vec]];
	res
	};

ifft:{[vec]
	n:count vec 1;
	inverseScale:1%n;
	vecR:((vec 1;vec 0));
	invr:fft[vecR];
	inverseScale * ((invr 1);(invr 0))
	};



fftBluestein:{[vec]
	// This method will evaluate the FFT by a convolution of two transforms of length
	// that is the first power of 2 greater than 2n-1
	// This allows for radix-2 fft to be used, for the two forward transformns, and the 
	// final reverse transform of the result

	// First we will define some utility definitions 
	PI::acos -1;
	n:count vec 0;
	indexOrig:til n;
	m:4 * 2 sv (1b,((count 2 vs n)-1)#0b);

	// Precompute the Trig values
	angles:{[n;x]PI*x%n}[n;] ((indexOrig * indexOrig) mod (2*n));
	blueTwiddle:(cos angles;neg sin angles);

	// Define temporary vectors, which are padded with 0's to make a length of m
	// tempA has 0's prepended
	// tempB is defined mathematically with negative indices, the result is 0's in the middle
	tempA: flip (flip mult[vec;blueTwiddle]),(flip ((m-n)#0f;(m-n)#0f));
	tempB:(((blueTwiddle 0),((m-((2*n)-1))#0f),(1_(blueTwiddle 0)));((blueTwiddle 1),((m-((2*n)-1))#0f),(1_(blueTwiddle 1))));

	// Using the circular convolution method. Compute the FFT of tempA and tempB, then scale the inverse transformation
	ft1:fftrad2[tempA];
	ft2:fftrad2[tempB];
	tempProd:mult[ft1;ft2];

	// Clear variables
	delete tempA, ft1, ft2 from `.;


	// Inverse the results, then scale and take the first n results
	ires:fftrad2[(tempProd 1;tempProd 0)];
	ires: (ires 1;ires 0);
	convolved:(1%m) * ((n#(ires 0));(n#(ires 1)));
	
	mult[convolved; blueTwiddle]
	};

\d .