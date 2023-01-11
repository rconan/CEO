import numpy as np

class IIRfilter:
    """
    Class to simulate an IIR type of filter, i.e. 
    a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M] - a[1]*y[n-1] - ... - a[N]*y[n-N]
    where a[] is the list of denominator coefficients, and b[] is the list of numerator coefficients.
    
    Parameters:
    -----------
    num_coeff : list of floats or numpy array
        Numerator coefficients for the filter
    den_coeff : list of floats or numpy array
        Denominator coefficients for the filter
    """
    
    def __init__(self,num_coeff, den_coeff):
        self.num_coeff = num_coeff
        self.den_coeff = den_coeff
        self.coeffmatrix = np.array(np.concatenate((num_coeff, [0], -den_coeff[1:])))
        self._outidx = len(self.num_coeff)
        self._state_initialized = False
        
    def initState(self,input_size=(1,)):
        """
        Allocate the state vector given the size of the input.
        
        Parameters:
        -----------
        input_size : tuple containing an integer
            Size of input vector. Default: 1
        """
        if type(input_size) != tuple:
            raise Exception("input_size must be a tuple.")
        state_size = ((len(self.num_coeff) + len(self.den_coeff)),) + input_size
        self._state = np.zeros(shape=state_size)
        self._state_initialized = True
        
    def reset(self):
        """
        Reset the state vector of the IIR filter.
        """
        if self._state_initialized == True:
            self._state *= 0.0
        
    def apply_filter(self, _input_):
        """
        Compute the output of the filter given a provided input at a particular instant.
        
        Parameters:
        ----------
        _input_ : numpy array
            The input to the IIR filter at a given instant.
        """
        self._state = np.roll(self._state,1,axis=0)
        self._state[0,:] = _input_
        _output_ = self.coeffmatrix @ self._state
        self._state[self._outidx,:] = _output_
        return np.squeeze(_output_)
    
    def modifyState(self, modeidx, values):
        """
        Modify the state of the state vector.
        
        Parameters:
        -----------
        modeidx : list of integers
            Index to elements of the output vector to be modified.
        values : numpy array (same size as modeidx)
            Values of elements of the output vector to be modified.
        """
        self._state[self._outidx:,modeidx] += values
