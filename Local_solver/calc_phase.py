"""
    Local Solver

    A Python program to extract phases of the reflection of interest from shotgather window 
    Base on Xinfa Zhu SEG abstract with title "COMPARISON OF METHODS FOR MODELING PHASE VARIATION WITH ANGLE' FROM 2011"
    and Bram Willemsen and Alison Malcolm "An efficient coupled acoustic-elastic local solver applied to phase inversion"
   
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Bram Willemsen, Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""
import numpy as np
import warnings
from scipy.signal import hilbert

__all__ = ['calc_phase',
           'envelope_trace_scipy',
           'envelope_trace_self',
           'calc_phase_trace',
           'correlation_coefficient',
           'rotate_phase'
           ]

def calc_phase(data, ret_corr_coeffs=False): 
    """
    Calculates the phase from the window shot-gather.

    Parameters:

    * data: 2d-arrays
        Arrays with traces in the window.
 
    Returns:

    * best_phase: 1d-arrays
        Phase extracted for each trace.

    """
    sh = data.shape
    if len(sh) == 1: #single trace
        if ret_corr_coeffs:
            corr_coeffs = np.zeros(360)
        else:
            corr_coeffs = None
            
        best_phase = calc_phase_trace(data, corr_coeffs)
        
        if ret_corr_coeffs:
            return best_phase, corr_coeffs
        else:
            return best_phase
        
    if len(sh) == 2: 
    
        [_,nr] = sh
        best_phases = np.zeros(nr)
        
        if ret_corr_coeffs:
            corr_coeffs = np.zeros((nr,360))
        else:
            corr_coeffs = None        
        
        corr_coeffs_iter = None
        for ir in range(nr):
            if ret_corr_coeffs:
                corr_coeffs_iter = corr_coeffs[ir,:]
                
            best_phases[ir] = calc_phase_trace(data[:,ir], corr_coeffs_iter)

        if ret_corr_coeffs:
            return best_phases, corr_coeffs
        else:
            return best_phases
            
    else:
        raise Exception('Wrong input?')

def envelope_trace_scipy(trace):
    """
    Calculates the envelope phase from each trace. 
    (-imag part is hilbert. real part is original signal.
    The hilbert is has zero freq zero (zero mean) in this implementation 
    So abs gives amplitude complex valued analytic function)

    Parameters:

    * trace: 1d-arrays
        Array with each trace.
 
    Returns:

    * np.abs(hilbert(trace))
        absolute value of the Hilbert trace 

    """    
    return np.abs(hilbert(trace))

def envelope_trace_self(trace):
    """
    Calculates the envelope from each trace. 
    #(alternative to getting directly from the 'hilbert' function. I put zero frequency to zero as well. 
    #This gives an envelope that is much better for certain functions like the box-car. 
    #Without zeroing the mean of the hilbert, the envelope will deviate much from zero far away from the box due to the mean shift

    Parameters:

    * trace: 1d-arrays
        Array with each trace.
 
    Returns:

    * np.sqrt(trace**2 + hilbert_of_trace**2) 

    """ 
    
    hilbert_of_trace = rotate_phase(trace,90)
    hilbert_of_trace -= np.mean(hilbert_of_trace)
    return np.sqrt(trace**2 + hilbert_of_trace**2) 
    
def calc_phase_trace(trace, corr_coeffs = None):

    """
    Calculates phase from each trace. 
    #WARNING SECTION
    warnings.warn('Be careful when the mean is not zero. 
    Box-car example works correctly, but not sure if this will generally be the case.')
    
    Parameters:

    * trace: 1d-arrays
        Array with each trace.
 
    Returns:

    * (360 - best_match_shift)%360
    
    """ 
    
    
    trace_zero_mean = trace 
    
    #Prepare correlation coefficient array.
    if type(corr_coeffs)==np.ndarray:
        if corr_coeffs.size == 360:
            corr_coeffs[:] = 0 #initialize
        else:
            raise Exception('wrong input')
    else:
        corr_coeffs = np.zeros(360)
    
    for shift in range(360): #Doing 1 degree at a time is somewhat arbitrary
        rotated_trace = rotate_phase(trace_zero_mean, shift)
        
        # The hilbert transform is made zero mean. 
        # This is the default in scipy, and also seems to result in better envelopes in the boxcar example where 
        # The mean is significantly nonzero
        
        envelope_of_rotated_trace = envelope_trace_scipy(rotated_trace) 
        
        corr_coeffs[shift] = correlation_coefficient(rotated_trace, envelope_of_rotated_trace) 

       
    best_match_shift = np.argmax(corr_coeffs) 
    return (360 - best_match_shift)%360 
    
        
def correlation_coefficient(a, b):
    """
    Define correlation coefficient following Fomel and van der Baan, 2010
    Parameters:

    * a and b: two ndarrays of the same length
 
    Returns:

    *  corr_coeff
    
    """ 

    corr_coeff = a.dot(b) / np.sqrt( a.dot(a) * b.dot(b))  
    return corr_coeff

def rotate_phase(trace,shift_deg):
    """
    Rotated_trace
    Parameters:
    * trace: one ndarrays
    * shift_deg: shift rotation range 
 
    Returns:

    *  rotated_trace
    
    """      
    shift_rad = shift_deg/180.*np.pi
    
    #Do FFT
    trace_fft = np.fft.fft(trace)
    freq = np.fft.fftfreq(trace_fft.size)
    freq_shifted = np.fft.fftshift(freq)
    zero_index = np.where(freq_shifted == 0)[0][0]
    trace_fft_shifted = np.fft.fftshift(trace_fft)
    
    neg_shift = np.exp(-1j*shift_rad)
    pos_shift = np.exp( 1j*shift_rad)
    
    #Shift phase by 'shift'
    rotated_trace_fft_shifted = trace_fft_shifted 
    rotated_trace_fft_shifted[0:zero_index]    *= neg_shift
    rotated_trace_fft_shifted[(zero_index+1):] *= pos_shift 
    
    
    #do IFFT
    rotated_trace_fft = np.fft.ifftshift(rotated_trace_fft_shifted)
    rotated_trace = np.fft.ifft(rotated_trace_fft)
    return np.real(rotated_trace) 
    
