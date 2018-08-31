# Signal Processing and Q
## Acknowledgements
I would like to acknowledge the Machine Learning Repository that is hosted by the University of California, Irvine for providing the real-world datasets that I used for research and as example sets used within this paper. This repository has an enormous number of freely-available datasets that can be accessed from http://archive.ics.uci.edu/ml/datasets.html.

<abbr title="Tooltip text">Text</abbr>

## Signal Processing and Q
Signal Processing is the analysis, interpretation and manipulation of signals to reveal important information. Signals of interest include human speech, seismic waves, images of faces or handwriting, brainwaves, radar, traffic counts and many others. This processing reveals information in a signal that can be obscured by non-useful information, commonly called ‘noise’. This noise can be due to the stochastic  nature of signals and/or interference from other signals. Traditionally, signal processing this has been performed by either a physical analogue process, through the use of a series of electronic circuits, or through the use of dedicated hardware solutions, such as SIP (Systems In Package) or SoC (Systems on a Chip).

The use of Internet of Things  (IoT) devices to capture signals is driving the trend towards software based signal processing solutions. Software based solutions are not only cheaper and more widely accessible than their hardware alternatives, but their highly configurable nature is better suited to the modular nature of IoT sensor setups. As a result there is an increased availability of cheaply processed signal data, enabling more data driven decision making, particularly in the manufacturing sector [1]. 

Currently, popular software implementations of digital signal processing techniques can found in the open source-source libraries of Python (e.g., Sci Py, NumPy, PyPlot) and C++ (e.g., SlgPack, Aquila, Liquid-DSP ). The convenience of these libraries is offset by the lack of a robust, high throughput data capture system, such as kdb+ tick. Whilst it is possible to integrate a q process with Python or C++ , and utilize these signal processing routines on a q dataset, it is entirely possible to implement them natively within q.

This whitepaper will explore how statistical signal processing operations (those which assume that signals are stochastic), can be implemented natively within q to remove noise, extract useful information, and quickly identify anomalies. This integration allows for q/kdb+ to be used as a single platform for the capture, processing, analysis and storage of large volumes of sensor data.

In order to run the q code given in this paper, it is necessary to place all the functions and constants within a namespace, which is most easily achieved by adding a directory declaration at the begging and end of the code. The technical specifications for the machine used to perform the computations described in this paper are as follows
•	CPU: Intel® Core™i7 Quad Core Processor i7-7700 (3.6GHz) 8MB Cache
•	RAM: 16GB Corsair 2133MHz SODIMM DDR4 (1 x 16GB)
•	OS: Windows 10 64-bit


## Basic Steps for Signal Processing
For the purpose of this paper, Signal Processing will be composed of the follow steps:
1.	Data Capture
2.	Spectral Analysis
3.	Smoothing
4.	Anomaly Detection

This paper will focus on the last three steps listed above, as the first topic of Data capture has already been covered extensively in previous Kx whitepapers, including (Kdb+tick profiling for throughput optimization, Data recovery for kdb+tick and Query Routing: a kdb+ framework for a scalable load-balanced). We will investigate a sensor dataset (see Figure 1) collected by the University of California, Irvine, which contains the total power load for a single house over 4 years, with a resolution of 1 minute.


![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure1.png "Figure 1: The grid load data of a single household, in 1 minute increments, from January 2007 to September 2010")

## Spectral Analysis
A fundamental property of any wave (and hence signal) is that of superpositioning , which is the combining of different signals (of the same or different frequency), to create a new wave with different amplitude and/or frequency. This means that any complicated signal is really just the result of the combination of more fundamental, simple signals. Spectral analysis is a broad term for the family of transformations that ‘decompose’  a signal from the time domain to the frequency domain ,  revealing these fundamental components that make up a more complex signal.

This results in a series of complex vectors associated to a frequency 'bin', whose absolute value is representative of the relative strength of that frequency within the original signal. This is a powerful tool that allows insight to be gained on the periodic nature of the signal, to troubleshoot the sensor in question, and to guide the subsequent transformation of the signal to clean it up.

In general, some of the basic rules to keep in mind when using the frequency decomposition of a signal are:
* Sharp distinct lines indicate the presence of a strong periodic nature;
* Wide peaks indicate there is periodic nature, and possibly some spectral leakage (which won't be fully discussed in this paper);
* An overall trend in the frequency distribution indicates aperiodic (and hence of an infinite period) behaviour present.

Consider the following artificial signal in Figure 2, which is a combination of a 10Hz and 20Hz signal and a constant background Gaussian noise:

```pi:acos -1;
x:(2*pi*(1%2048)) * til 2048;
y:(5 *  sin (20*x))+(5 *  sin (3 * x)) + `float$({first 1?15.0} each x);
```
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure2.png "Figure 2: A simulated noise signal")
This time-series does show a clear periodic nature, but the details are hard to discern from visual inspection. The following graph (Figure 3) shows the frequency distribution of this signal, produced through the application of an un-windowed Fourier Transform, a technique that will be covered in detailed in the chapter on the Fast Fourier Transform. 

![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure3.png "Figure 3: The frequency distribution of the simulated signal show in Figure 2")

The existence of distinct, sharp peaks and a low level random ‘static’ in Figure 3 demonstrate the existence of both a period trends and some background noise within the signal respectively. But why are there 3 different frequencies (i.e., 3, 20 and 50Hz) present? This is because the original signal contained a low-level 50Hz signal, one that you might typically observe if a sensor wasn't correctly shielded from 50Hz main power. This “leakage” would be influencing any decisions made from this signal, and may not have been immediately detected without spectral analysis.

For the purpose of this paper, spectral analysis is used to provide information to make a better decision on how to process the signal. In the case above, an analyst may decide to implement a high pass filter to reduce the 50Hz noise, and then implement a moving mean (or Finite Impulse) filter to reduce the noise.

In general, this spectral analysis can provide essential information to the analyst on the periodic nature of the signal, and hence appropriate techniques for detecting anomalies and smoothing the signal . 

In the following chapters, an industry standard technique for spectral analysis, the radix-2 Fast Fourier Transform is implemented natively within q. This implementation is reliant on the introduction of a framework for complex numbers within q.


### Complex Numbers and Q
All spectral analysis methods are performed in the complex number space. Whilst complex numbers are not natively supported in q, it is relatively straightforward to create a valid complex vector space by representing any complex number as a list of reals and imaginary parts. 

![equation](http://latex.codecogs.com/gif.latex?z%3D5&plus;3i)

```
z:(5;3)
```
This means that a complex vector can be represented like so,

![equation](http://latex.codecogs.com/gif.latex?z%3D%28%285&plus;3i%29%2C%287-9i%29%2C%2821&plus;2i%29%29)

```
\d .signal
mult:{[vec1;vec2]
	// Performs the dot product of two vectors, 
	realOut:((vec1 0) * (vec2 0)) - ((vec1 1) * (vec2 1));
	imagOut:((vec1 1) * (vec2 0)) + ((vec1 0) * (vec2 1));
	(realOut;imagOut)};

division:{[vec1;vec2]

	denom:1%((vec1 1) xexp 2) + ((vec2 1) xexp 2);
	mul:mult[(vec1 0);(vec1 1);(vec2 0);neg (vec2 1)];
	realOut:(mul 0)*denom; 
	imagOut:(mul 1)*denom;
	(realOut;imagOut)};

conj:{[vec](vec 0;neg vec 1)};

mag:{[vec]
	sqrvec:(vec xexp 2);
	sqrt (sqrvec 0)+(sqrvec 1)};
\d .

q).signal.mult[(5 -3);(9 2)]
51 -17
q).signal.mult[(5 2 1;-3 -8 10);(9 8 -4;2 3 6)]
51  40  -64   / Reals
-17 -58 -34   / Imaginary

```

#### Fast Fourier Transform
The Fourier Transform is a standard transformation to decompose a real or complex signal into its frequency distribution. The Fast Fourier Transform is actually a family of algorithms that allow for the Fourier transform to be computed in an efficient manner by utilizing symmetries within the computation and removing redundant calculations. This reduces the complexity of the algorithm from O(n2) to around O(n log(n)), which scales impressively for large samples.

Traditionally the Fast Fourier Transform and its family of transforms (such as Bluestein Z-Chirp or Raders), are often packaged up in libraries for languages such as Python or Matlab, or as external libraries such as FFTW (Fasted Fourier Transform in the West). However, with the creation of a Complex Number structure in q, this family of transformations can be written natively in q and therefore benefit from the speed of vector operations that q has over other languages.

The following shows how the Radix-2 FFT (specifically a Decimation-In-Time, Bit-Reversed Input Radix-2 algorithm) can be implemented in q, this function takes equi length real and imaginary components of a signal (a real signal of length N would have an N length, 0 valued imaginary component) and produces a complex vector of length N, representing the frequency decomposition of the signal.


``` 
\d .signals
// Global Definitions
PI:acos -1; / pi ;
BR:2 sv reverse 2 vs til 256; / Reversed bits in bytes 00-FF
P2:1,prds[7#2],1; / Powers of 2 with outer guard
P256:prds 7#256; / Powers of 256

bitreversal:{[indices] 
    // Applies a bitwise reversal for a list of indices
    ct:ceiling 2 xlog last indices; / Number of significant bits (assumes sorted)
    bt:BR 256 vs indices; / Breakup into bytes and reverse bits
    dv:P2[8-ct mod 8]; / Divisor to shift leading bits down
    bt[0]:bt[0] div dv; / Shift leading bits
    sum bt*count[bt]#1,P256 div dv / Reassemble bytes back to integer
    };

fftrad2:{[vec]
    // This performs a FFT for samples of power 2 size
    // First, define some constants
    n:count vec 0;
    n2:n div 2;
    indexOrig:til n;

    // Twiddle Factors - precomputed the complex values over the 
    //discrete angles, using Euler formula
    angles:{[n;x]2*PI*x%n}[n;] til n2;
    const.twiddle:(cos angles;neg sin angles);

     // Bit-reversed the vector and define it into a namespace so lambdas can access it 
    ind:bitreversal[indexOrig];
    .res.vec:`float$vec . (0 1;ind);

    // Precomputing the indices required to immplement each temporal phase of the DIT
	// Number of signals
	signalcount:`int$ {2 xexp x} 1+ til ( `int$(2 xlog n));
	// Number of points in each signal
	signalpoints:reverse signalcount;
	// Define an initial count, this can be used as a 
	// backbone to get the even and odd indices by adjustment
	initial:{[n2;x]2*x xbar til n2}[n2;] peach (signalcount div 2);
	evens:{[n2;x]x + (n2)#til n2 div count distinct x}[n2;] peach initial;
	odds:evens (+)' (signalcount div 2);
	twiddleIndices:{[n;n2;x](n2)#(0.5*x) * til (n div (x))}[n;n2;] peach signalpoints;


    // Butterfly Implementation
	bflyComp:{[bfInd] 
        tmp:mult[.res.vec[;bfInd 1];.const.twiddle[;bfInd 2]];
        .[`.res.vec;(0 1;bfInd 1);:;.res.vec[;bfInd 0]-tmp];
        .[`.res.vec;(0 1;bfInd 0);+;tmp]};

    bflyComp each (flip `int$(evens;odds;twiddleIndices));
    res.vec};
\d .

```
#### Sampling Frequency and Length
The number of samples Ns used and the frequency at which they occur fs (determined by the constant time period between samples), are the important factors in the spectral analysis as they determine the spacing between detected frequencies and the maximum frequency that can be observed. In short, these factors scale the x-axis of any transformation.

The number of frequency bins (i.e., x-axis ticks) that a given transformation will produce is equal to the number of samples used, and as a result of the complex math used, will range from negative to positive values. The Nyquist Theorem states that for a Fourier Transformation, the absolute maximum frequency that can be detected is half that of the sampling frequency f_s. Given this information, the x-axis of a frequency distribution can be scaled by using the number and frequency of samples used to produce it.
xAxis:{[Ns;fs] (neg fs*0.5)+(fs%(Ns-1))*til Ns}


```    
xAxis:{[Ns;fs] (neg fs*0.5)+(fs%(Ns-1))*til Ns}
```

#### Windowing the Signal
There is an important assumption that is carried forward from the continuous Fourier transform, to the discrete version, that the data set being operated on is infinitely large. Obviously this cannot be true for a discrete transform, but it is simulated by operating with a circular topology, i.e., by assuming that the end of the sample connects to the start. This means that there can be discontinuities in the data, if a non-integer number of wavelengths are captured by the sampling data. The resulting effect is called ‘spectral leakage’, where results in the frequency domain are diffused over adjacent frequency bins.

![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure4.png "Figure 4: A simple sinusoid signal, with a red overlay for the signal that is being sampled.")

![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure5.png "Figure 5: The signal that would be used by a Fourier transform, from the sample shown in Figure 4")

The solution to this problem is to window the signal which adjusts the overall amplitude of the signal to limit the discontinuities. For this paper, the commonly used Hanning Window (which finds a middle ground between many different windowing functions),  shown in Figure 6, will be applied. 

![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure6.png "Figure 6: The shape of a Hann Window, which is applied to an input signal to 'tamper down' the signal so that there is greater continuity under the transform. Distributed under CC Public Domain")

#### The Fourier Transformation for a Real-Valued Signal
In general, the final step to getting the frequency distribution from the transformed vector series is to get the normalized magnitude of the complex result. The normalised magnitude is computed by taking the absolute value of the real and imaginary components, divided by the number of samples. However, for a real-valued input, where there is symmetry about the central axis, the dataset can be halved by only using the first half of the result.

### Putting it all together
A single function will be defined to perform a spectral analysis (through a Fourier Transform) on a complex vector. This function should window the vector, perform a FFT, and the scale the results appropriately.  
```
spectral:{[vec;fs]
    // Get the length of the series, ensure it is a power of 2
    // If not, we should abort operation and output a reason
    nr:count first vec;
    if[not 0f  =(2 xlog nr) mod 1;0N! "The input vector is not a power of 2";`error];
    
    // Window the functions, using the Hanning window
    wndFactor:{[nr;x](sin((x*PI)%(nr-1))) xexp 2} [nr;] each til nr;
    vec*:(wndFactor;wndFactor);

    // Compute the fft and get the absolute value of the result
    fftRes:.fft.fftrad2[vec];
    mag:.cmplx.mag fftRes;

    // Scale an x axis
    xAxis:{[Ns;fs] (neg fs*0.5)+(fs%(Ns-1))*til Ns}[nr;fs];
    
    // Table these results
    ([]freq:xAxis;magnitude:mag%nr)
    };

```

Let's apply this to the load data, loading it from the txt files and then filling in null values to ensure a constant sampling frequency.
```
q) household:select Time:Date+Time,Power:fills Global_active_power from 
    ("DTF      "; enlist ";") 0:`:household_power_consumption.txt 
```
Let's sample at 2-hour intervals, ensuring the sample size is a power of 2 and perform the spectral analysis
```
q)household2hour:select from household where 0=(i mod 120)
q)n:count household2hour
17294
q)cn:`int$2 xexp (count 2 vs n)-1
16384
q)spec:spectral[(cn sublist household2hour`Power;cn#0f);1%(120*60)]
q)\t spectral[(cn sublist household2hour`Power;cn#0f);1%(120*60)]
6

```
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure7.png "The frequency distribution of the household power load, excluding the first result to demonstrate the Hermitian symmetry")
The above FFT was formed on 16,384 points, and was computed in 6ms . In this case, the frequencies are so small that they are difficult to make inferences from, it's hard to comprehend  what a frequency of 1e-6 Hz corresponds to. A nice solution is to remove the symmetry present (taking only the first half), then convert from Frequency to Period.
```
q)specReal:select period:neg (1%freq)%3600, magnitude from ((cn div 2) sublist spec)
```
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure 8.png "Figure 8: The frequency distribution of the household grid data with the symmetry removed")

It can see from the above plot that the first value is nearly an order of magnitude larger than any other, this suggests that there is a strong, high frequency signal that hasn’t been captured. At the very least we know for certain that there is a 50Hz signal due to mains power, the existence of others isn’t a shock. After the initial large value, there are distinct 5, 6, 12 and 24 hour cycles, which could be explained as part of the ebb and flow of day-to-day life, coming home from school and work. This result, with strong sharp lines, definitively shows the existence of periodic trends, along with a constant background noise. This knowledge will be used to further process this signal and obtain more useful, long-term information.

## Comparing the FFT
In the following chapter, a comparison between two alternative implementations of FFT will be performed to benchmark the native q based algorithm. First a comparison will be made to a similar algorithm written natively in Python3.6, then a comparison will be made to calling into a well refined, and highly optimized C based library, via the q fusion with python – embedPy.

In order to get an accurate comparison, each method should be tested on a sufficiently large dataset, the largest 2N subset of the previously used grid data is suitable for this purpose. 

```
//Get the largest possible 2N dataset size
q)n:count household;
q)nr:`int$2 xexp (count 2 vs n)-1;
comparison:select real:Power, imag:nr#0f from (nr#household)
save `comparison.csv
q)\t .fft.fftrad2[(comparison`real;comparison`imag)]
1020
```
So our q algorithm can determine the Fourier Transform of a million point dataset in about a second, a nice whole number to gauge the comparative performance. 

In python, the native complex math definitions field can be used to build a similar radix-2 FFT Decimation in Time routine as was implemented above, from a publically available GitHub https://github.com/peterhinch/micropython-fft/blob/master/algorithms.py

```python
import pandas as pd
import numpy as np
import time
import math
from cmath import exp, pi

data = pd.read_csv('c:/q/w64/comparison.csv')
fftinput = data.real.values + 1j * data.imag.values
def ffft(nums, forward=True, scale=False):
    n = len(nums)
    m = int(math.log2(n))
    #n= 2**m #Calculate the number of points
    #Do the bit reversal
    i2 = n >> 1
    j = 0
    for i in range(n-1):
        if i<j: nums[i], nums[j] =  nums[j], nums[i]
        k = i2
        while (k <= j):
            j -= k
            k >>= 1
        j+=k
    #Compute the FFT
    c = 0j-1
    l2 = 1
    for l in range(m):
        l1 = l2
        l2 <<= 1
        u = 0j+1
        for j in range(l1):
            for i in range(j, n, l2):
                i1 = i+l1
                t1 = u*nums[i1]
                nums[i1] = nums[i] - t1
                nums[i] += t1
            u *= c
        ci = math.sqrt((1.0 - c.real) / 2.0) # Generate complex roots of unity
        if forward: ci=-ci                   # for forward transform
        cr = math.sqrt((1.0 + c.real) / 2.0) # phi = -pi/2 -pi/4 -pi/8...
        c = cr + ci*1j#complex(cr,ci)
    # Scaling for forward transform
    if (scale and forward):
        for i in range(n):
            nums[i] /= n
    return nums  

t0 = time.time()
ffft(fftinput)
t1 = time.time()
print(t1 - t0)
8.336997509002686
```

And finally using embedPy to call upon the numpy C based library to perform our FFT, which will likely not use a radix-2 algorithm, but a mixed radix implementation that allows for multi-threading (which the in place Decimation-In-Time algorithm fundamentally cannot do). 

```
p)import numpy as np;
p)import pandas as pd;
p)import time;

p)data = pd.read_csv('comparison.csv');


// Convert to a complex pair

p)fftinput = data.real.values + 1j * data.imag.values;
p)t0 = time.time();
p)np.abs(np.fft.fft (fftinput));
p)t1 = time.time();
p)print (t1 - t0)
0.04644131660461426
```
And so the end results are somewhat expected. When implementing a simplistic algorithm, the performance of q is nearly a magnitude faster than that seen by a natively implemented python algorithm. However, it can be seen that calling into a well defined C library from q (if the environment is already setup to allow for embedPy to be used), yields a significant improvement in computation time, which is to be expected. It would be unreasonable to expect an algorithm developed initially to be done by hand in 1805 to outperform a modern, highly optimized library of functions. 
## Smoothing
From the spectral analysis above, it can see that the most dominant signals have a period of at most 24 hours. A simple method to see the longer-term trends and remove the daily noise, is to build a filter that is targeted at those frequencies, a low pass filter. A simple low pass filter can be constructed as a windowed moving average, which for q removes the need for any conclusions due to the inbuilt function mavg. All that is needed is to ensure the moving average is centred on the given result, which can be achieved by constructing a function that uses rotate to shift the input lists to appropriately centre the result. The 

```
\d .signal
movAvg:{[list;N]
    $[0=N mod 2;
        (floor N%2) rotate 2 mavg (N mavg list);
        (floor N%2) rotate N mavg list]};
\d .
```

Let's implement a 24-hour moving average on the data
```
q)householdMA:select Time, Power: .signal.movAvg[Power;24*60] from household
q)\t householdMA:select Time, Power: .signal.movAvg[Power;24*60] from household
109
```
The following figures (Figure 9 and Figure 10) show the result of the application of a 1,440 point moving average on over 2 million data points (achieved in just over 100ms).

The results are as follows
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure9.png "Figure 9: The complete smoothed grid data")
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure10.png "Figure 10: The smoothed dataset zoomed in to better highlight the increased clarity in the longer-term trends ")


Comparison of the smoothed data (Figure 9) and the raw data (Figure 10) demonstrates the clarity achieved by applying a simple moving average smoothing. The use of a carefully selected moving average (24 hours in this case), has removed a significant amount of noise present in the data shown the existence of a year long trend in power consumption. 

Application of a FFT to this smoothed data set can reveal trends hidden within this smoothed dataset, which would have been difficult to obtain from the raw data. For this, a lower sampling rate will be used to assist this process, sampling at 6 hour intervals instead of 2. 

![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure11.png "Fourier Transform")

The figure above shows the scaled results of this second FFT. It can be seen that whilst there is still a good bit of high frequency noise, there is a definite 3.5 day trend, and a smaller, but distinct 7 day trend in the power usage. In the context of household power consumption, this result is not unexpected.

## Anomaly Detection
Anomaly detection can be implemented independently of any other signal processing, serving as an early warning indicator of deviations from expected behaviour. Assuming the sensors are behaving stochastically, outliers can be identified by the level of deviation they are experiencing from the expected behaviour, i.e., using the mean and standard deviation. Ideally, these deviations should be measured from the local behaviour in the signal, which allows for periodic behaviour to be accounted for. This can be achieved by creating a windowed moving average and standard deviation, which conveniently allows the reuse the moving average filtering code.

A moving standard deviation and moving average function can be written, taking advantage of the pre-existing and optimized mdev and mavg functions.

```
\d .signal
movDev:{[list;N]
    $[0=N mod 2;
	(floor N%2) rotate 2 mavg (N mavg list);
	(floor N%2) rotate N mavg list]};
\d .
```
We can determine if an anomaly is present if the actual reading is more than a factor (sigma) * standard deviation, away from the moving mean. To implement this effectively we can use a vector conditional, let's use a 5-point moving window and a sigma level of 3.

```
q)sigma:3;
q)outliers:update outlier:?[(Power > movAvg[Power;5] + (fills sigma * movDev[Power;5]))|
    (Power < movAvg[Power;5] - (fills sigma *movDev[Power;5]));Power;0n] from household;
q)\t update outlier:?[(Power > movAvg[Power;5] + (fills sigma * movDev[Power;5]))|
    (Power < movAvg[Power;5] - (fills sigma *movDev[Power;5]));Power;0n] from household;
380
```
It takes approximately 380ms to perform this detection routine on over 2 million data points, and the results can be seen below.
![alt text](https://github.com/callumjbiggs/q-signals/blob/master/images/figure12.png "Figure 12: The results of an anomaly detection on the original grid data")


## Conclusion
In this paper, we have shown the ease at which simple, fast signal processing can be implemented natively within q. All that was required was a basic framework to support complex number operations to be constructed to allow for complex field mathematics within q. Once the framework was constructed, the signal processing algorithms themselves can easily be implemented in such a manner that takes advantage of the speed at which vector operations can be implemented in q.

This complex framework was essential to the successful implementation of a Fast Fourier Transform (FFT), which is an essential operation in the processing of many signals. This FFT was roughly 10 times faster than an algorithm working natively in Python, and only 10 times slower than calling to a highly optimized C library of FFT algorithms. Overall, this demonstrates the ease at which operations that required complex numbers, can be integrated natively into q. 

The collection of Signal Processing Functions that were shown in this paper were able to rapidly process the grid load dataset obtained from the University of California, Irvine. A FFT was implemented on a sample of 16,000 data points in around 6ms and guided from that spectral analysis, a simple low pass filter was implemented directly in the time domain, on over 2 million data points, in around 100ms. The resulting data had a clear reduction in noise and revealed the longer-term trends present in the signal. It has also been shown that a simple and robust anomaly detection routine can be implemented through q, allowing anomalies to be detected in large datasets very rapidly or alternatively, real-time anomaly detection to be implemented with other data capture and processing routines.

Clearly, native application of Signal Processing, which historically have been the realm of libraries accessed through Python or C++, can be natively integrated into q/kdb+ data systems. This would allow for data capture, real-time signal processing, data storage and recall from the same integrated platform.

## References
McElheran, E. B. (2016). Data in Action: Data-Driven Decision Making in U.S. Manufacturing. Center for Economic Studies, U.S. Census Bureau.

