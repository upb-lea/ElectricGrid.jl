import math
import numpy as np
import cupy as cp
from numpy import any, asarray, concatenate, cos, delete, \
    empty, exp, eye, isinf, ones, pad, sin, zeros, squeeze
from numpy.random import rand, randn
from numpy.linalg import solve, eigvals, matrix_rank
from numpy.linalg.linalg import LinAlgError
import scipy as sp
from scipy.signal import cont2discrete
from scipy.signal import StateSpace as signalStateSpace
from warnings import warn
#from .lti import LTI, common_timebase, isdtime, _process_frequency_response
#from . import config
from copy import copy, deepcopy



__all__ = ['StateSpace', 'ss']


def _ssmatrix(data, axis=1, bit32=False):
    # Convert the data into an array or matrix, as configured
    # If data is passed as a string, use (deprecated?) matrix constructor

    if bit32:
        arr = np.array(data, dtype=np.float32)
    else:
        arr = np.array(data, dtype=float)

    ndim = arr.ndim
    shape = arr.shape
    # Change the shape of the array into a 2D array
    if (ndim > 2):
        raise ValueError("state-space matrix must be 2-dimensional")
    elif (ndim == 2 and shape == (1, 0)) or \
         (ndim == 1 and shape == (0, )):
        # Passed an empty matrix or empty vector; change shape to (0, 0)
        shape = (0, 0)
    elif ndim == 1:
        # Passed a row or column vector
        shape = (1, shape[0]) if axis == 1 else (shape[0], 1)
    elif ndim == 0:
        # Passed a constant; turn into a matrix
        shape = (1, 1)
    #  Create the actual object used to store the result
    return arr.reshape(shape)



class StateSpace():
    # Allow ndarray * StateSpace to give StateSpace._rmul_() priority
    __array_priority__ = 11     # override ndarray and matrix types
    def __init__(self, *args, **kwargs):

        if 'bit32' in kwargs:
            self.bit32 = kwargs['bit32']
        if 'use_cuda' in kwargs:
            self.use_cuda = kwargs['use_cuda']
        if 'offline_expm' in kwargs:
            self.offline_expm = kwargs['offline_expm']

        # first get A, B, C, D matrices
        if len(args) == 4:
            # The user provided A, B, C, and D matrices.
            (A, B, C, D) = args
        elif len(args) == 5:
            # Discrete time system
            (A, B, C, D, _) = args
        elif len(args) == 1:
            # Use the copy constructor.
            if not isinstance(args[0], StateSpace):
                raise TypeError(
                    "The one-argument constructor can only take in a "
                    "StateSpace object. Received %s." % type(args[0]))
            A = args[0].A
            B = args[0].B
            C = args[0].C
            D = args[0].D
        else:
            raise ValueError(
                "Expected 1, 4, or 5 arguments; received %i." % len(args))
        # Convert all matrices to standard form
        A = _ssmatrix(A, bit32=self.bit32)
        # if B is a 1D array, turn it into a column vector if it fits
        if np.asarray(B).ndim == 1 and len(B) == A.shape[0]:
            B = _ssmatrix(B, axis=0, bit32=self.bit32)
        else:
            B = _ssmatrix(B, bit32=self.bit32)
        if np.asarray(C).ndim == 1 and len(C) == A.shape[0]:
            C = _ssmatrix(C, axis=1, bit32=self.bit32)
        else:
            C = _ssmatrix(C, axis=0, bit32=self.bit32)    # if this doesn't work, error below
        if np.isscalar(D) and D == 0 and B.shape[1] > 0 and C.shape[0] > 0:
            # If D is a scalar zero, broadcast it to the proper size
            D = np.zeros((C.shape[0], B.shape[1]))
        D = _ssmatrix(D, bit32=self.bit32)

        self.ninputs = D.shape[1]
        self.noutputs = D.shape[0]
        self.A = A
        self.B = B
        self.C = C
        self.D = D

        # now set dt
        if len(args) == 4:
            if 'dt' in kwargs:
                dt = kwargs['dt']
            elif self._isstatic():
                dt = None
            else:
                dt = 0
        elif len(args) == 5:
            dt = args[4]
            if 'dt' in kwargs:
                warn("received multiple dt arguments, "
                     "using positional arg dt = %s" % dt)
        elif len(args) == 1:
            try:
                dt = args[0].dt
            except AttributeError:
                if self._isstatic():
                    dt = None
                else:
                    dt = 0
        self.dt = dt
        self.nstates = A.shape[1]

        if self.offline_expm and dt is not None and dt == 0:
            test_dt = 1e-4
            M = np.block([[A * test_dt, B * test_dt, np.zeros((self.nstates, self.ninputs))],
                          [np.zeros((self.ninputs, self.nstates + self.ninputs)),
                           np.identity(self.ninputs)],
                          [np.zeros((self.ninputs, self.nstates + 2 * self.ninputs))]])
            expM = sp.linalg.expm(M)
            self.Ad = expM[:self.nstates, :self.nstates]
            self.Bd1 = expM[:self.nstates, self.nstates + self.ninputs:]
            self.Bd0 = expM[:self.nstates, self.nstates:self.nstates + self.ninputs] - self.Bd1

        if 0 == self.nstates:
            # static gain
            # matrix's default "empty" shape is 1x0
            A.shape = (0, 0)
            B.shape = (0, self.ninputs)
            C.shape = (self.noutputs, 0)

        # Check that the matrix sizes are consistent.
        if self.nstates != A.shape[0]:
            raise ValueError("A must be square.")
        if self.nstates != B.shape[0]:
            raise ValueError("A and B must have the same number of rows.")
        if self.nstates != C.shape[1]:
            raise ValueError("A and C must have the same number of columns.")
        if self.ninputs != B.shape[1]:
            raise ValueError("B and D must have the same number of columns.")
        if self.noutputs != C.shape[0]:
            raise ValueError("C and D must have the same number of rows.")


    ninputs = 0
    noutputs = 0
    nstates = 0
    A = []
    B = []
    C = []
    D = []
    bit32 = False
    use_cuda = False
    offline_expm = False
    Ad = []
    Bd1 = []
    Bd0 = []

    def issiso(self):
        '''Check to see if a system is single input, single output'''
        return self.ninputs == 1 and self.noutputs == 1

    def isctime(self, strict=False):
        if self.dt is None:
            return True if not strict else False
        return self.dt == 0

    def isdtime(self, strict=False):
        if self.dt == None:
            return True if not strict else False
        return self.dt > 0

    def __str__(self):
        """Return string representation of the state space system."""
        string = "\n".join([
            "{} = {}\n".format(Mvar,
                               "\n    ".join(str(M).splitlines()))
            for Mvar, M in zip(["A", "B", "C", "D"],
                               [self.A, self.B, self.C, self.D])])
        if self.isdtime(strict=True):
            string += f"\ndt = {self.dt}\n"
        return string

    def __getitem__(self, indices):
        """Array style access"""
        if len(indices) != 2:
            raise IOError('must provide indices of length 2 for state space')
        i = indices[0]
        j = indices[1]
        return StateSpace(self.A, self.B[:, j], self.C[i, :],
                          self.D[i, j], self.dt)

    def _isstatic(self):
        """True if and only if the system has no dynamics, that is,
        if A and B are zero. """
        return not np.any(self.A) and not np.any(self.B)


def ss(*args, **kwargs):
    if len(args) == 4 or len(args) == 5:
        return StateSpace(*args, **kwargs)
    elif len(args) == 1:
        sys = args[0]
        if isinstance(sys, StateSpace):
            return deepcopy(sys)
        else:
            raise TypeError("ss(sys): sys must be a StateSpace or "
                            "TransferFunction object.  It is %s." % type(sys))
    else:
        raise ValueError(
            "Needs 1, 4, or 5 arguments; received %i." % len(args))


# Helper function for checking array-like parameters
def _check_convert_array(in_obj, legal_shapes, err_msg_start, squeeze=False,
                         transpose=False):
    # convert nearly everything to an array.
    out_array = np.asarray(in_obj)
    if (transpose):
        out_array = np.transpose(out_array)

    # Test element data type, elements must be numbers
    legal_kinds = set(("i", "f", "c"))  # integer, float, complex
    if out_array.dtype.kind not in legal_kinds:
        err_msg = "Wrong element data type: '{d}'. Array elements " \
                  "must be numbers.".format(d=str(out_array.dtype))
        raise TypeError(err_msg_start + err_msg)

    # If array is zero dimensional (in_obj is scalar):
    # create array with legal shape filled with the original value.
    if out_array.ndim == 0:
        for s_legal in legal_shapes:
            # search for shape that does not contain the special symbol any.
            if "any" in s_legal:
                continue
            the_val = out_array[()]
            out_array = np.empty(s_legal, 'd')
            out_array.fill(the_val)
            break

    # Test shape
    def shape_matches(s_legal, s_actual):
        """Test if two shape tuples match"""
        # Array must have required number of dimensions
        if len(s_legal) != len(s_actual):
            return False
        # All dimensions must contain required number of elements. Joker: "all"
        for n_legal, n_actual in zip(s_legal, s_actual):
            if n_legal == "any":
                continue
            if n_legal != n_actual:
                return False
        return True

    # Iterate over legal shapes, and see if any matches out_array's shape.
    for s_legal in legal_shapes:
        if shape_matches(s_legal, out_array.shape):
            break
    else:
        legal_shape_str = " or ".join([str(s) for s in legal_shapes])
        err_msg = "Wrong shape (rows, columns): {a}. Expected: {e}." \
                  .format(e=legal_shape_str, a=str(out_array.shape))
        raise ValueError(err_msg_start + err_msg)

    # Convert shape
    if squeeze:
        out_array = np.squeeze(out_array)
        # We don't want zero dimensional arrays
        if out_array.shape == tuple():
            out_array = out_array.reshape((1,))

    return out_array


def _process_labels(labels, signal, length):
    if labels is None or len(labels) == 0:
        return None

    # See if we got passed a dictionary (from iosys)
    if isinstance(labels, dict):
        # Form inverse dictionary
        ivd = {v: k for k, v in labels.items()}

        try:
            # Turn into a list
            labels = [ivd[n] for n in range(len(labels))]
        except KeyError:
            raise ValueError("Name dictionary for %s is incomplete" % signal)

    # Convert labels to a list
    labels = list(labels)

    # Make sure the signal list is the right length and type
    if len(labels) != length:
        raise ValueError("List of %s labels is the wrong length" % signal)
    elif not all([isinstance(label, str) for label in labels]):
        raise ValueError("List of %s labels must all  be strings" % signal)

    return labels

def _process_time_response(
        tout, yout, issiso=False, transpose=None, squeeze=None):
    if squeeze is None:
        squeeze = None

    # Figure out whether and how to squeeze output data
    if squeeze is True:         # squeeze all dimensions
        yout = np.squeeze(yout)
    elif squeeze is False:      # squeeze no dimensions
        pass
    elif squeeze is None:       # squeeze signals if SISO
        if issiso:
            if yout.ndim == 3:
                yout = yout[0][0]       # remove input and output
            else:
                yout = yout[0]          # remove input
    else:
        raise ValueError("Unknown squeeze value")

    # See if we need to transpose the data back into MATLAB form
    if transpose:
        # Transpose time vector in case we are using np.matrix
        tout = np.transpose(tout)

        # For signals, put the last index (time) into the first slot
        yout = np.transpose(yout, np.roll(range(yout.ndim), 1))

    # Return time, output, and (optionally) state
    return tout, yout



class TimeResponseData():
    def __init__(
            self, time, outputs, states=None, inputs=None, issiso=None,
            output_labels=None, state_labels=None, input_labels=None,
            transpose=False, return_x=False, squeeze=None, multi_trace=False
    ):
        # Time vector
        self.t = np.atleast_1d(time)
        if self.t.ndim != 1:
            raise ValueError("Time vector must be 1D array")

        #
        # Output vector (and number of traces)
        #
        self.y = np.array(outputs)

        if self.y.ndim == 3:
            multi_trace = True
            self.noutputs = self.y.shape[0]
            self.ntraces = self.y.shape[1]

        elif multi_trace and self.y.ndim == 2:
            self.noutputs = 1
            self.ntraces = self.y.shape[0]

        elif not multi_trace and self.y.ndim == 2:
            self.noutputs = self.y.shape[0]
            self.ntraces = 0

        elif not multi_trace and self.y.ndim == 1:
            self.noutputs = 1
            self.ntraces = 0

            # Reshape the data to be 2D for consistency
            self.y = self.y.reshape(self.noutputs, -1)

        else:
            raise ValueError("Output vector is the wrong shape")

        # Check and store labels, if present
        self.output_labels = _process_labels(
            output_labels, "output", self.noutputs)

        # Make sure time dimension of output is the right length
        if self.t.shape[-1] != self.y.shape[-1]:
            raise ValueError("Output vector does not match time vector")

        #
        # State vector (optional)
        #
        # If present, the shape of the state vector should be consistent
        # with the multi-trace nature of the data.
        #
        if states is None:
            self.x = None
            self.nstates = 0
        else:
            self.x = np.array(states)
            self.nstates = self.x.shape[0]

            # Make sure the shape is OK
            if multi_trace and \
               (self.x.ndim != 3 or self.x.shape[1] != self.ntraces) or \
               not multi_trace and self.x.ndim != 2 :
                raise ValueError("State vector is the wrong shape")

            # Make sure time dimension of state is the right length
            if self.t.shape[-1] != self.x.shape[-1]:
                raise ValueError("State vector does not match time vector")

        # Check and store labels, if present
        self.state_labels = _process_labels(
            state_labels, "state", self.nstates)

        #
        # Input vector (optional)
        #
        # If present, the shape and dimensions of the input vector should be
        # consistent with the trace count computed above.
        #
        if inputs is None:
            self.u = None
            self.ninputs = 0

        else:
            self.u = np.array(inputs)

            # Make sure the shape is OK and figure out the nuumber of inputs
            if multi_trace and self.u.ndim == 3 and \
               self.u.shape[1] == self.ntraces:
                self.ninputs = self.u.shape[0]

            elif multi_trace and self.u.ndim == 2 and \
                 self.u.shape[0] == self.ntraces:
                self.ninputs = 1

            elif not multi_trace and self.u.ndim == 2 and \
                 self.ntraces == 0:
                self.ninputs = self.u.shape[0]

            elif not multi_trace and self.u.ndim == 1:
                self.ninputs = 1

                # Reshape the data to be 2D for consistency
                self.u = self.u.reshape(self.ninputs, -1)

            else:
                raise ValueError("Input vector is the wrong shape")

            # Make sure time dimension of output is the right length
            if self.t.shape[-1] != self.u.shape[-1]:
                raise ValueError("Input vector does not match time vector")

        # Check and store labels, if present
        self.input_labels = _process_labels(
            input_labels, "input", self.ninputs)

        # Figure out if the system is SISO
        if issiso is None:
            # Figure out based on the data
            if self.ninputs == 1:
                issiso = (self.noutputs == 1)
            elif self.ninputs > 1:
                issiso = False
            else:
                # Missing input data => can't resolve
                raise ValueError("Can't determine if system is SISO")
        elif issiso is True and (self.ninputs > 1 or self.noutputs > 1):
            raise ValueError("Keyword `issiso` does not match data")

        # Set the value to be used for future processing
        self.issiso = issiso

        # Keep track of whether to squeeze inputs, outputs, and states
        if not (squeeze is True or squeeze is None or squeeze is False):
            raise ValueError("Unknown squeeze value")
        self.squeeze = squeeze

        # Keep track of whether to transpose for MATLAB/scipy.signal
        self.transpose = transpose

        # Store legacy keyword values (only needed for legacy interface)
        self.return_x = return_x

    def __call__(self, **kwargs):
        response = copy(self)

        # Update any keywords that we were passed
        response.transpose = kwargs.pop('transpose', self.transpose)
        response.squeeze = kwargs.pop('squeeze', self.squeeze)
        response.return_x = kwargs.pop('return_x', self.squeeze)

        # Check for new labels
        input_labels = kwargs.pop('input_labels', None)
        if input_labels is not None:
            response.input_labels = _process_labels(
                input_labels, "input", response.ninputs)

        output_labels = kwargs.pop('output_labels', None)
        if output_labels is not None:
            response.output_labels = _process_labels(
                output_labels, "output", response.noutputs)

        state_labels = kwargs.pop('state_labels', None)
        if state_labels is not None:
            response.state_labels = _process_labels(
                state_labels, "state", response.nstates)

        # Make sure no unknown keywords were passed
        if len(kwargs) != 0:
            raise ValueError("Unknown parameter(s) %s" % kwargs)

        return response

    @property
    def time(self):

        """Time vector.

        Time values of the input/output response(s).

        :type: 1D array"""
        return self.t

    # Getter for output (implements squeeze processing)
    @property
    def outputs(self):
        """Time response output vector.

        Output response of the system, indexed by either the output and time
        (if only a single input is given) or the output, trace, and time
        (for multiple traces).  See :attr:`TimeResponseData.squeeze` for a
        description of how this can be modified using the `squeeze` keyword.

        :type: 1D, 2D, or 3D array

        """
        t, y = _process_time_response(
            self.t, self.y, issiso=self.issiso,
            transpose=self.transpose, squeeze=self.squeeze)
        return y

    # Getter for states (implements squeeze processing)
    @property
    def states(self):
        """Time response state vector.

        Time evolution of the state vector, indexed indexed by either the
        state and time (if only a single trace is given) or the state, trace,
        and time (for multiple traces).  See :attr:`TimeResponseData.squeeze`
        for a description of how this can be modified using the `squeeze`
        keyword.

        :type: 2D or 3D array

        """
        if self.x is None:
            return None

        elif self.squeeze is True:
            x = self.x.squeeze()

        elif self.ninputs == 1 and self.noutputs == 1 and \
             self.ntraces == 1 and self.x.ndim == 3 and \
             self.squeeze is not False:
            # Single-input, single-output system with single trace
            x = self.x[:, 0, :]

        else:
            # Return the full set of data
            x = self.x

        # Transpose processing
        if self.transpose:
            x = np.transpose(x, np.roll(range(x.ndim), 1))

        return x

    # Getter for inputs (implements squeeze processing)
    @property
    def inputs(self):
        """Time response input vector.

        Input(s) to the system, indexed by input (optiona), trace (optional),
        and time.  If a 1D vector is passed, the input corresponds to a
        scalar-valued input.  If a 2D vector is passed, then it can either
        represent multiple single-input traces or a single multi-input trace.
        The optional ``multi_trace`` keyword should be used to disambiguate
        the two.  If a 3D vector is passed, then it represents a multi-trace,
        multi-input signal, indexed by input, trace, and time.

        See :attr:`TimeResponseData.squeeze` for a description of how the
        dimensions of the input vector can be modified using the `squeeze`
        keyword.

        :type: 1D or 2D array

        """
        if self.u is None:
            return None

        t, u = _process_time_response(
            self.t, self.u, issiso=self.issiso,
            transpose=self.transpose, squeeze=self.squeeze)
        return u

    # Getter for legacy state (implements non-standard squeeze processing)
    @property
    def _legacy_states(self):
        """Time response state vector (legacy version).

        Time evolution of the state vector, indexed indexed by either the
        state and time (if only a single trace is given) or the state,
        trace, and time (for multiple traces).

        The `legacy_states` property is not affected by the `squeeze` keyword
        and hence it will always have these dimensions.

        :type: 2D or 3D array

        """

        if self.x is None:
            return None

        elif self.ninputs == 1 and self.noutputs == 1 and \
             self.ntraces == 1 and self.x.ndim == 3:
            # Single-input, single-output system with single trace
            x = self.x[:, 0, :]

        else:
            # Return the full set of data
            x = self.x

        # Transpose processing
        if self.transpose:
            x = np.transpose(x, np.roll(range(x.ndim), 1))

        return x

    # Implement iter to allow assigning to a tuple
    def __iter__(self):
        if not self.return_x:
            return iter((self.time, self.outputs))
        return iter((self.time, self.outputs, self._legacy_states))

    # Implement (thin) getitem to allow access via legacy indexing
    def __getitem__(self, index):
        # See if we were passed a slice
        if isinstance(index, slice):
            if (index.start is None or index.start == 0) and index.stop == 2:
                return (self.time, self.outputs)

        # Otherwise assume we were passed a single index
        if index == 0:
            return self.time
        if index == 1:
            return self.outputs
        if index == 2:
            return self._legacy_states
        raise IndexError

    # Implement (thin) len to emulate legacy testing interface
    def __len__(self):
        return 3 if self.return_x else 2



# Forced response of a linear system
def forced_response(sys, T=None, U=0., X0=0., transpose=False,
                    interpolate=False, return_x=None, squeeze=None):

    # If return_x was not specified, figure out the default
    if return_x is None:
        return_x = False

    A, B, C, D = np.asarray(sys.A), np.asarray(sys.B), np.asarray(sys.C), \
        np.asarray(sys.D)
    # d_type = A.dtype
    n_states = A.shape[0]
    n_inputs = B.shape[1]
    n_outputs = C.shape[0]

    # Convert inputs to numpy arrays for easier shape checking
    if U is not None:
        U = np.asarray(U)
    if T is not None:
        # T must be array-like
        T = np.asarray(T)

    # Set and/or check time vector in discrete time case
    if sys.isdtime():
        if T is None:
            if U is None or (U.ndim == 0 and U == 0.):
                raise ValueError('Parameters ``T`` and ``U`` can\'t both be '
                                 'zero for discrete-time simulation')
            # Set T to equally spaced samples with same length as U
            if U.ndim == 1:
                n_steps = U.shape[0]
            else:
                n_steps = U.shape[1]
            dt = 1. if sys.dt in [True, None] else sys.dt
            T = np.array(range(n_steps)) * dt
        else:
            # Make sure the input vector and time vector have same length
            if (U.ndim == 1 and U.shape[0] != T.shape[0]) or \
                    (U.ndim > 1 and U.shape[1] != T.shape[0]):
                raise ValueError('Parameter ``T`` must have same elements as'
                                 ' the number of columns in input array ``U``')
            if U.ndim == 0:
                U = np.full((n_inputs, T.shape[0]), U)
    else:
        if T is None:
            raise ValueError('Parameter ``T`` is mandatory for continuous '
                             'time systems.')

    # Test if T has shape (n,) or (1, n);
    T = _check_convert_array(T, [('any',), (1, 'any')],
                             'Parameter ``T``: ', squeeze=True,
                             transpose=transpose)

    n_steps = T.shape[0]            # number of simulation steps

    # equally spaced also implies strictly monotonic increase,
    dt = (T[-1] - T[0]) / (n_steps - 1)
    if not np.allclose(np.diff(T), dt):
        raise ValueError("Parameter ``T``: time values must be equally "
                         "spaced.")

    # create X0 if not given, test if X0 has correct shape
    X0 = _check_convert_array(X0, [(n_states,), (n_states, 1)],
                              'Parameter ``X0``: ', squeeze=True)

    # Test if U has correct shape and type
    legal_shapes = [(n_steps,), (1, n_steps)] if n_inputs == 1 else \
        [(n_inputs, n_steps)]
    U = _check_convert_array(U, legal_shapes,
                             'Parameter ``U``: ', squeeze=False,
                             transpose=transpose)

    if sys.bit32:
        if X0.dtype != np.float32:
            X0 = X0.astype(np.float32)
        if U.dtype != np.float32:
            U = U.astype(np.float32)

        xout = np.zeros((n_states, n_steps), dtype=np.float32)
        yout = np.zeros((n_outputs, n_steps), dtype=np.float32)
        dt = dt.astype(np.float32)
    else:
        xout = np.zeros((n_states, n_steps))
        yout = np.zeros((n_outputs, n_steps))

    xout[:, 0] = X0


    # Separate out the discrete and continuous time cases
    if sys.isctime(strict=True):
        # Solve the differential equation, copied from scipy.signal.ltisys.

        # Faster algorithm if U is zero
        # (if not None, it was converted to array above)
        if U is None or np.all(U == 0):
            # Solve using matrix exponential
            expAdt = sp.linalg.expm(A * dt)
            for i in range(1, n_steps):
                xout[:, i] = expAdt @ xout[:, i-1]
            yout = C @ xout

        # General algorithm that interpolates U in between output points
        else:
            # convert input from 1D array to 2D array with only one row
            if U.ndim == 1:
                U = U.reshape(1, -1)  # pylint: disable=E1103

        # Algorithm: to integrate from time 0 to time dt, with linear
            # interpolation between inputs u(0) = u0 and u(dt) = u1, we solve
            #   xdot = A x + B u,        x(0) = x0
            #   udot = (u1 - u0) / dt,   u(0) = u0.
            #
            # Solution is
            #   [ x(dt) ]       [ A*dt  B*dt  0 ] [  x0   ]
            #   [ u(dt) ] = exp [  0     0    I ] [  u0   ]
            #   [u1 - u0]       [  0     0    0 ] [u1 - u0]

            if sys.offline_expm:
                Ad = sys.Ad
                Bd1 = sys.Bd1
                Bd0 = sys.Bd0
            else:
                M = np.block([[A * dt, B * dt, np.zeros((n_states, n_inputs))],
                             [np.zeros((n_inputs, n_states + n_inputs)),
                              np.identity(n_inputs)],
                             [np.zeros((n_inputs, n_states + 2 * n_inputs))]])
                expM = sp.linalg.expm(M)
                Ad = expM[:n_states, :n_states]
                Bd1 = expM[:n_states, n_states+n_inputs:]
                Bd0 = expM[:n_states, n_states:n_states + n_inputs] - Bd1

            if sys.bit32:
                Ad = Ad.astype(np.float32)
                Bd1 = Bd1.astype(np.float32)
                Bd0 = Bd0.astype(np.float32)

            if sys.use_cuda:
                Ad = cp.asarray(Ad)
                Bd1 = cp.asarray(Bd1)
                Bd0 = cp.asarray(Bd0)
                C = cp.asarray(C)
                D = cp.asarray(D)
                U = cp.asarray(U)
                xout = cp.asarray(xout)
                yout = cp.asarray(yout)

            for i in range(1, n_steps):
                if sys.use_cuda:
                    xout[:, i] = cp.dot(Ad, xout[:, i - 1]) + cp.dot(Bd0, U[:, i - 1]) + cp.dot(Bd1, U[:, i])
                else:
                    #xout[:, i] = (Ad @ xout[:, i-1] + Bd0 @ U[:, i-1] + Bd1 @ U[:, i])
                    xout[:, i] = np.dot(Ad, xout[:, i - 1]) + np.dot(Bd0, U[:, i - 1]) + np.dot(Bd1, U[:, i])
            #yout = C @ xout + D @ U
            if sys.use_cuda:
                yout = cp.dot(C, xout) + cp.dot(D, U)
            else:
                yout = np.dot(C, xout) + np.dot(D, U)

            if sys.use_cuda:
                U = cp.asnumpy(U)
                xout = cp.asnumpy(xout)
                yout = cp.asnumpy(yout)
                del D, C, Bd0, Bd1, Ad
                cp.get_default_memory_pool().free_all_blocks()

        tout = T

    else:
        # Discrete type system => use SciPy signal processing toolbox

        # sp.signal.dlsim assumes T[0] == 0
        spT = T - T[0]

        if sys.dt is not True and sys.dt is not None:
            # Make sure that the time increment is a multiple of sampling time

            # First make sure that time increment is bigger than sampling time
            # (with allowance for small precision errors)
            if dt < sys.dt and not np.isclose(dt, sys.dt):
                raise ValueError("Time steps ``T`` must match sampling time")

            # Now check to make sure it is a multiple (with check against
            # sys.dt because floating point mod can have small errors
            if not (np.isclose(dt % sys.dt, 0) or
                    np.isclose(dt % sys.dt, sys.dt)):
                raise ValueError("Time steps ``T`` must be multiples of "
                                 "sampling time")
            sys_dt = sys.dt

            # sp.signal.dlsim returns not enough samples if
            # T[-1] - T[0] < sys_dt * decimation * (n_steps - 1)
            # due to rounding errors.
            # https://github.com/scipyscipy/blob/v1.6.1/scipy/signal/ltisys.py#L3462
            scipy_out_samples = int(np.floor(spT[-1] / sys_dt)) + 1
            if scipy_out_samples < n_steps:
                # parantheses: order of evaluation is important
                spT[-1] = spT[-1] * (n_steps / (spT[-1] / sys_dt + 1))

        else:
            sys_dt = dt         # For unspecified sampling time, use time incr

        xout = np.transpose(xout)
        yout = np.transpose(yout)
        U = np.transpose(U)
        xout[0, :] = np.asarray(X0)
        stoptime = spT[-1]
        out_samples = int(np.floor(stoptime / dt)) + 1
        tout = np.linspace(0.0, stoptime, num=out_samples)
        tout = tout + T[0]
        if sys.bit32:
            tout = tout.astype(np.float32)

        if sys.use_cuda:
            A = cp.asarray(A)
            B = cp.asarray(B)
            C = cp.asarray(C)
            D = cp.asarray(D)
            U = cp.asarray(U)
            xout = cp.asarray(xout)
            yout = cp.asarray(yout)

        # Simulate the system
        for i in range(0, out_samples - 1):
            if sys.use_cuda:
                xout[i + 1, :] = (cp.dot(A, xout[i, :]) + cp.dot(B, U[i, :]))
                yout[i, :] = (cp.dot(C, xout[i, :]) + cp.dot(D, U[i, :]))
            else:
                xout[i + 1, :] = (np.dot(A, xout[i, :]) + np.dot(B, U[i, :]))
                yout[i, :] = (np.dot(C, xout[i, :]) + np.dot(D, U[i, :]))

        # Last point
        if sys.use_cuda:
            yout[out_samples - 1, :] = (cp.dot(C, xout[out_samples - 1, :]) + cp.dot(D, U[out_samples - 1, :]))
        else:
            yout[out_samples - 1, :] = (np.dot(C, xout[out_samples - 1, :]) + np.dot(D, U[out_samples - 1, :]))

        if sys.use_cuda:
            U = cp.asnumpy(U)
            xout = cp.asnumpy(xout)
            yout = cp.asnumpy(yout)
            del A, B, C, D
            cp.get_default_memory_pool().free_all_blocks()

        if not interpolate:
            # If dt is different from sys.dt, resample the output
            inc = int(round(dt / sys_dt))
            tout = T            # Return exact list of time steps
            yout = yout[::inc, :]
            xout = xout[::inc, :]
        else:
            # Interpolate the input to get the right number of points
            U = sp.interpolate.interp1d(T, U)(tout)

        # Transpose the output and state vectors to match local convention
        xout = np.transpose(xout)
        yout = np.transpose(yout)
        U = np.transpose(U)

    if sys.use_cuda:
        cp.cuda.Device().synchronize()

    return TimeResponseData(
        tout, yout, xout, U, issiso=sys.issiso(),
        transpose=transpose, return_x=return_x, squeeze=squeeze)
