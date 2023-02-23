import numpy as np
import copy

class MirrorPositioner:
    """
    Encapsulates the state of a Mirror (M1 or M2) Positioner state (Rxyz, Txyz), taking into account of:
        1. Initial state of the Positioner;
        2. Input command (with the possibility to simulate a given temporal response);
        3. Disturbance to the Positioner (expressed as Rxyz, Txyz disturbances).
    
    Parameters:
    ----------
        mirror : gmtMirrors object
            M1 or M2 mirror of the GMT class
        
        temporal_response : string
            Type of temporal response of the M2 Positioner: "ideal", "lpf butter", "lpf bessel". Default: "ideal"
        
        filter_order : integer
            Filter order when selecting a LPF as temporal response. Default: 1
        
        cutoff_freq : float
            Cut-off frequency in Hz of low-pass filter When selecting a LPF as temporal_response.
        
        Ts : float
            Sampling time (needed for realizing the digital temporal response filter).
    """
    
    def __init__(self, mirror, temporal_response='ideal', filter_order=1, cutoff_freq=None, Ts=None):
        
        if not 'gmtMirrors' in str(type(mirror)):
            raise TypeError('"mirror" must be of ceo.gmtMirrors derived class.')
        self._mirror = mirror
        
        #--- setup the buffer to hold the initial state of the positioner.
        self.init_state = self.__newStateVector()
        
        #--- setup the buffer to "sample and hold" the command vector to the positioner.
        self.comm_state = self.__newStateVector()
        
        #--- setup the buffer to hold the time-evolving state vector (when simulating a temporal response)
        self.filt_state = self.__newStateVector()

        #--- setup the buffer to simulate a disturbance vector to the positioner.
        self.disturb_state = self.__newStateVector()
        
        #-- setup the temporal response model for the positioner.
        self.set_temporal_response(temporal_response, filter_order=filter_order, cutoff_freq=cutoff_freq, Ts=Ts)
    
    
    def set_temporal_response(self, temporal_response, filter_order=1, cutoff_freq=None, Ts=None):
        """
        Sets the temporal response model of the positioner.
        
        Parameters:
        -----------
            temporal_response : string
                Type of temporal response of the M2 Positioner: "ideal", "lpf butter", "lpf bessel". Default: "ideal"

            filter_order : integer
                Filter order when selecting a LPF as temporal response. Default: 1

            cutoff_freq : float
                Cut-off frequency in Hz of low-pass filter When selecting a LPF as temporal_response.

            Ts : float
                Sampling time (needed for realizing the digital temporal response filter).
        """
        if temporal_response not in ['ideal', 'lpf butter', 'lpf bessel']:
            raise ValueError("'temporal_response' must be either ['ideal','lpf butter', 'lpf bessel'].")
        self.temporal_response = temporal_response
        
        if temporal_response != 'ideal':
            from scipy import signal
            from ceo import IIRfilter
            
            if cutoff_freq is None or Ts is None:
                raise ValueError("'cutoff_freq [Hz]' and Ts [s] need to be specified.")
            
            lpf_w = cutoff_freq / (1/(2*Ts)) # Normalize the frequency
            
            if temporal_response == 'lpf butter':
                lpf_b, lpf_a = signal.butter(filter_order, lpf_w, 'low', analog=False)
            elif temporal_response == 'lpf bessel':
                lpf_b, lpf_a = signal.bessel(filter_order, lpf_w, 'low', analog=False)

            self._filter = IIRfilter(lpf_b, lpf_a)
            self._filter.initState(input_size=(6*7,)) 


    def update_temporal_response(self):
        """
        Update state according to the temporal response for the stored command.
        """
        if self.temporal_response != 'ideal':
            Txyz = self.comm_state['Txyz']
            Rxyz = self.comm_state['Rxyz']
            _cin_ = np.concatenate((Txyz, Rxyz), axis=1).flatten()
            _cout_ = self._filter.apply_filter(_cin_)
            Txyz, Rxyz = np.split(np.reshape(_cout_, [7,6]), [3], axis=1)
            self.filt_state['Txyz'] = Txyz
            self.filt_state['Rxyz'] = Rxyz
        else:
            self.filt_state = self.comm_state
            
    
    def __newStateVector(self, Txyz=np.zeros((7,3)), Rxyz=np.zeros((7,3))):
        """
        Create a state vector to be used internally.
        """
        return mirror_state_vector( dict(Txyz=Txyz, Rxyz=Rxyz) )
    
    
    def initState(self, Txyz=None, Rxyz=None):
        """
        Initialize the state of the positioner with the provided Txyz and Rxyz values.
        """
        if Txyz is not None:
            if not np.array_equal( Txyz.shape, (7,3)):
                raise ValueError("Txyz must be a (7,3) array.")
            self.init_state['Txyz'] = Txyz
            
        if Rxyz is not None:
            if not np.array_equal( Rxyz.shape, (7,3)):
                raise ValueError("Rxyz must be a (7,3) array.")            
            self.init_state['Rxyz'] = Rxyz
    

    def command(self, Txyz=None, Rxyz=None):
        """
        Sample and Hold the absolute command to the positioner with the provided Txyz and Rxyz.
        """
        if Txyz is not None:
            if not np.array_equal( Txyz.shape, (7,3)):
                raise ValueError("Txyz must be a (7,3) array.")
            self.comm_state['Txyz'] = Txyz

        if Rxyz is not None:
            if not np.array_equal( Rxyz.shape, (7,3)):
                raise ValueError("Rxyz must be a (7,3) array.")            
            self.comm_state['Rxyz'] = Rxyz

    
    def update(self, Txyz=None, Rxyz=None):
        """
        Update the state of the Mirror Positioner using as disturbance the provided Txyz and Rxyz.
        The updated state will be equal to: init_state + command state (with temporal dynamics) + disturbance.
        """
        
        #--- Update disturbance state vector if present
        if Txyz is not None:
            if not np.array_equal( Txyz.shape, (7,3)):
                raise ValueError("Txyz must be a (7,3) array.")
            self.disturb_state['Txyz'] = Txyz
        else:
            self.disturb_state['Txyz'] *= 0  # disturbance buffer shall not have memory

        if Rxyz is not None:
            if not np.array_equal( Rxyz.shape, (7,3)):
                raise ValueError("Rxyz must be a (7,3) array.")            
            self.disturb_state['Rxyz'] = Rxyz
        else:
            self.disturb_state['Rxyz'] *= 0  # disturbance buffer shall not have memory
        
        #--- Update temporal response of positioner to stored command
        self.update_temporal_response()
        
        #--- Update the state of the positioner
        updated_state = self.init_state + self.filt_state + self.disturb_state       
        self._mirror ^= updated_state.state

    
    def reset(self):
        self._mirror.motion_CS.reset()
        self.comm_state *= 0
        self.filt_state *= 0
        self.init_state *= 0
        self.disturb_state *= 0
        try:
            self._filter.reset()
        except:
            pass
        
        
class mirror_state_vector:
    """
    Class that represents a state vector.
    Parameters
    ----------
    mirror_state : dict
        Dictionary that defines the mirror state vector. It should have the following structure:
        mirror_state = {'dof1': np.array(<shape1>), 'dof2': np.array(<shape2>), ...}
    """
    def __init__(self, mirror_state):
        if type(mirror_state) is not dict:
            raise TypeError('state_template must be a dictionary.')
        self.state = copy.deepcopy(mirror_state)

    #================= State vector arithmetics ====================
    def __add__(self, other_state):
        sum_state = copy.deepcopy(self.state)
        for dof in sum_state:
            sum_state[dof] += other_state[dof]
        return mirror_state_vector(sum_state)

    def __sub__(self, other_state):
        sum_state = copy.deepcopy(self.state)
        for dof in sum_state:
            sum_state[dof] -= other_state[dof]
        return mirror_state_vector(sum_state)
    
    def __mul__(self, constant):
        sum_state = copy.deepcopy(self.state)
        for dof in sum_state:
            sum_state[dof] *= constant
        return mirror_state_vector(sum_state)    
    
    def __rmul__(self, constant):
        self.__mul__(constant)

    def __str__(self):
        return str(self.state)

    def __getitem__(self, key):
        return self.state[key]

    def __setitem__(self, key, value):
        self.state[key][:] = value
