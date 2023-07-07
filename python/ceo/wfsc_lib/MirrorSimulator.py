import numpy as np
import copy

class MirrorSimulator:
    """
    Encapsulates the state of a Mirror (M1 or M2) state (Rxyz, Txyz, modes), taking into account of:
        1. Initial state of the Mirror;
        2. Input command (with the possibility to simulate a given temporal response);
        3. Disturbance to the mirror (expressed as Rxyz, Txyz, and modes disturbances).
    
    Parameters:
    ----------
        mirror : gmtMirrors object
            M1 or M2 mirror of the GMT class
        
        segment_positioner : bool
            if True, include Txyz and Rxyz in the mirror state.
        
        segment_modes : bool
            if True, include modes in the mirror state.
        
        pos_temporal_response : string
            Type of temporal response of the mirror positioner: "ideal", "lpf butter", "lpf bessel". Default: "ideal"
        
        pos_filter_order : integer
            Filter order when selecting a LPF as temporal response of the mirror positioner. Default: 1
        
        pos_cutoff_freq : float
            Cut-off frequency in Hz of the LPF simulating the temporal response of the mirror positioner.

        modes_temporal_response : string
            Type of temporal response of the mirror modes: "ideal", "lpf butter", "lpf bessel". Default: "ideal"
        
        modes_filter_order : integer
            Filter order when selecting a LPF as temporal response of the mirror modes. Default: 1
        
        modes_cutoff_freq : float
            Cut-off frequency in Hz of the LPF simulating the temporal response of the mirror modes.
            
        Ts : float
            Sampling time (needed for realizing the digital temporal response filter).
    """
    def __init__(self, mirror, segment_positioner=True, segment_modes=True, Ts=None,
                 pos_temporal_response='ideal', pos_filter_order=1, pos_cutoff_freq=None,
                 modes_temporal_response='ideal', modes_filter_order=1, modes_cutoff_freq=None):
        
        if not 'gmtMirrors' in str(type(mirror)):
            raise TypeError('"mirror" must be of ceo.gmtMirrors derived class.')
        self._mirror = mirror
        
        #--- define state template and set temporal response model
        state_template = {}
        self.__simul_positioner = segment_positioner
        self.__simul_modes = segment_modes
        
        if segment_positioner==True:
            state_template['Txyz'] = np.zeros((7,3))
            state_template['Rxyz'] = np.zeros((7,3))
            self.set_positioner_response(pos_temporal_response, filter_order=pos_filter_order, 
                                         cutoff_freq=pos_cutoff_freq, Ts=Ts)
            
        if segment_modes==True:
            state_template['modes'] = np.zeros((7,mirror.modes.n_mode))
            self.set_modes_response(modes_temporal_response, filter_order=modes_filter_order, 
                                       cutoff_freq=modes_cutoff_freq, Ts=Ts)
        
        #--- setup the buffer to hold the initial state of the positioner.
        self.init_state = mirror_state_vector(state_template)
        
        #--- setup the buffer to "sample and hold" the command vector to the positioner.
        self.comm_state = mirror_state_vector(state_template)
        
        #--- setup the buffer to hold the time-evolving state vector (when simulating a temporal response)
        self.filt_state = mirror_state_vector(state_template)

        #--- setup the buffer to simulate a disturbance vector to the positioner.
        self.disturb_state = mirror_state_vector(state_template)
    
    
    def set_positioner_response(self, temporal_response, filter_order=1, cutoff_freq=None, Ts=None):
        """
        Sets the temporal response model of the positioner.
        
        Parameters:
        -----------
            temporal_response : string
                Type of temporal response of the positioner: "ideal", "lpf butter", "lpf bessel".

            filter_order : integer
                Filter order when selecting a LPF as temporal response. Default: 1

            cutoff_freq : float
                Cut-off frequency in Hz of low-pass filter When selecting a LPF as temporal_response.

            Ts : float
                Sampling time (needed for realizing the digital temporal response filter).
        """
        if temporal_response != 'ideal':
            self._pos_filter = get_temporal_filter(temporal_response, 7*6, filter_order, cutoff_freq, Ts)
            self._pos_response_params = dict(temporal_response=temporal_response, filter_order=filter_order,
                                         cutoff_freq=cutoff_freq, Ts=Ts)
        else:
            self._pos_response_params = dict(temporal_response=temporal_response)        
    
    
    def set_modes_response(self, temporal_response, filter_order=1, cutoff_freq=None, Ts=None):
        """
        Sets the temporal response model of the mirror modes.
        
        Parameters:
        -----------
            temporal_response : string
                Type of temporal response of the mirror modes: "ideal", "lpf butter", "lpf bessel".

            filter_order : integer
                Filter order when selecting a LPF as temporal response. Default: 1

            cutoff_freq : float
                Cut-off frequency in Hz of low-pass filter When selecting a LPF as temporal_response.

            Ts : float
                Sampling time (needed for realizing the digital temporal response filter).
        """
        if temporal_response != 'ideal':        
            self._modes_filter = get_temporal_filter(temporal_response, 7*self._mirror.modes.n_mode, 
                                        filter_order, cutoff_freq, Ts)
            self._modes_response_params = dict(temporal_response=temporal_response, filter_order=filter_order,
                                         cutoff_freq=cutoff_freq, Ts=Ts)
        else:
            self._modes_response_params = dict(temporal_response=temporal_response)
    
    
    def update_temporal_response(self):
        """
        Update state according to the temporal response for the stored command.
        """
        if self.__simul_positioner:
            if self._pos_response_params['temporal_response'] != 'ideal':
                Txyz = self.comm_state['Txyz']
                Rxyz = self.comm_state['Rxyz']
                _cin_ = np.concatenate((Txyz, Rxyz), axis=1).flatten()
                _cout_ = self._pos_filter.apply_filter(_cin_)
                Txyz, Rxyz = np.split(np.reshape(_cout_, [7,-1]), [3], axis=1)
                self.filt_state['Txyz'] = Txyz
                self.filt_state['Rxyz'] = Rxyz
            else:
                self.filt_state['Txyz'] = self.comm_state['Txyz']
                self.filt_state['Rxyz'] = self.comm_state['Rxyz']
        
        if self.__simul_modes:
            if self._modes_response_params['temporal_response'] != 'ideal':
                _cin_ = self.comm_state['modes'].flatten()
                _cout_ = self._modes_filter.apply_filter(_cin_)
                modes = np.reshape(_cout_, [7,-1])
                self.filt_state['modes'] = modes
            else:
                self.filt_state['modes'] = self.comm_state['modes']
    
    
    def initState(self, **state):
        """
        Initialize the state of the mirror.
        """
        for key, val in state.items():
            self.init_state[key] = val
    
    
    def command(self, **state):
        """
        Sample and Hold the absolute command to the positioner with the provided Txyz and Rxyz.
        """
        for key, val in state.items():
            self.comm_state[key] = val
    
    
    def update(self, **state):
        """
        Update the state of the Mirror Positioner using as disturbance the provided state.
        The updated state will be equal to: init_state + command state (with temporal dynamics) + disturbance.
        """
        
        #--- Update disturbance state vector if present
        self.disturb_state *= 0  # disturbance buffer shall not have memory
        for key, val in state.items():
            self.disturb_state[key] = val
        
        #--- Update temporal response of positioner to stored command
        self.update_temporal_response()
        
        #--- Update the state of the positioner
        updated_state = self.init_state + self.filt_state + self.disturb_state       
        self._mirror ^= updated_state.state

    
    def reset(self):
        if self.__simul_positioner:
            self._mirror.motion_CS.reset()
            if self._pos_response_params['temporal_response'] != 'ideal':
                self._pos_filter.reset()
        
        if self.__simul_modes:
            self._mirror.modes.reset()
            if self._modes_response_params['temporal_response'] != 'ideal':
                self._modes_filter.reset()
        
        self.comm_state *= 0
        self.filt_state *= 0
        self.init_state *= 0
        self.disturb_state *= 0

        
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


def get_temporal_filter(temporal_response, state_size, filter_order, cutoff_freq, Ts):
    """
    Sets the temporal response model of the mirror.
    """
    from scipy import signal
    from ceo.wfsc_lib import IIRfilter
    
    if temporal_response not in ['lpf butter', 'lpf bessel']:
        raise ValueError("'temporal_response' must be either ['ideal','lpf butter', 'lpf bessel'].")

    if cutoff_freq is None or Ts is None:
        raise ValueError("'cutoff_freq [Hz]' and Ts [s] need to be specified.")

    lpf_w = cutoff_freq / (1/(2*Ts)) # Normalize the frequency

    if temporal_response == 'lpf butter':
        lpf_b, lpf_a = signal.butter(filter_order, lpf_w, 'low', analog=False)
    elif temporal_response == 'lpf bessel':
        lpf_b, lpf_a = signal.bessel(filter_order, lpf_w, 'low', analog=False)

    _filter_ = IIRfilter(lpf_b, lpf_a)
    _filter_.initState(input_size=(state_size,)) 

    return _filter_