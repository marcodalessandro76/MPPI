"""
This module collects some tools to deal with the physics of Two Level Systems (TLS).
This modeling is particularly useful to describe both electromagnetic transition
between two *optically active* bands and magnetic system in which a time-independent
magnetic field induces oscillations between two (non degenerate) states with different
spin.

The feature of the model, as well as the derivation of the equation used in this
implementation of the physics of the TLS are described in this notebook_.

.. _notebook: tutorials/Model_TLS_optical_absorption.ipynb

"""

from scipy.integrate import odeint

def solveBlochEq_Xbasis(x0, time, Omega_abs, delta = 0):
    """
    Solve the Bloch equations in the rotating frame (in the basis of the x_i variables).

    Args:
        x0 (:py:class:`array`) : initial condition in the x_i basis
        time (:py:class:`array`) : array with the time values
        Omega_abs (`function`) : function with the module of the time dependent Rabi
            coupling expressed in the rotating frame, so only the envelope of the pulse
            has to be provided
        delta (:py:class:`float`) : difference between pulse energy and the energy shift
            of the two levels of the system. It has to be provided in the inverse dimensions
            of the time variable

    Returns:
        (:py:class:`array`) : array with the solution of the Bloch equations in the x_i basis.
                The index shape of the array is (3,len(time))

    """
    def Bloch_Eq(x, t, delta, Omega):
        dxdt = [delta*x[1],-delta*x[0]-Omega(t)*x[2],Omega(t)*x[1]]
        return dxdt
    x = odeint(Bloch_Eq, x0, time, args=(delta, Omega))
    return x.transpose()

def convertToXbasis(uprime, Omega0, invert = False):
    """
    Convert the Bloch vector u'_i components to the x_i basis.
    If the ``invert`` argument is True perform the inverse change of basis (from the x_i to the u'_i).

    Args:
        uprime (:py:class:`array`) : array with the Bloch vector. The function assumes that the shape of
            the array is (3,:) where the second index runs over time
        Omega0 (:py:class:`float`) : value of the complex Rabi coupling
        invert (:py:class:`bool`) : if True perform the inverse transformation

    Returns:
        (:py:class:`array`) : array with the transformed Bloch vector

    """
    x = np.zeros([3,len(uprime[0])])
    R_real = Omega0.real/abs(Omega0)
    R_imag = Omega0.imag/abs(Omega0)
    if invert: alpha = -1
    else: alpha = 1
    x[0,:] = R_imag*uprime[0,:] + alpha*R_real*uprime[1,:]
    x[1,:] = - alpha*R_real*uprime[0,:] + R_imag*uprime[1,:]
    x[2,:] = uprime[2,:]
    return x

def solveBlochEq(uprime0, time, Omega_abs, Omega0 = 1j, delta = 0):
    """
    Solve the Bloch equations in the rotating frame.

    Args:
        uprime0 (:py:class:`array`) : initial condition.
        time (:py:class:`array`) : array with the time values
        Omega_abs (`function`) : function with the module of the time dependent Rabi
            coupling expressed in the rotating frame, so only the envelope of the pulse
            has to be provided
        Omega0 (:py:class:`float`) : value of the complex Rabi coupling. It used to perform
            the transformation from the u'_i to the x_i (and its inverse after the equations)
            are solved. If the default value is used the solution is identical in both the
            formulation
        delta (:py:class:`float`) : difference between pulse energy and the energy shift
            of the two levels of the system. It has to be provided in the inverse dimensions
            of the time variable

    """
    bloch0 = np.zeros([3,1]) # add a fake time index to use the convertToXbasis function
    bloch0[:,0] = uprime0
    x0 = convertToXbasis(bloch0,Omega0)[:,0]
    x = solveBlochEq_Xbasis(x0,time,Omega_abs,delta=delta)
    uprime = convertToXbasis(x,Omega0,invert=True)
    return uprime

def convertToRotatingFrame(omega, time, u, invert = False):
    """
    Convert the Bloch vector components to the rotating frame (from the u_i to the u'_i).
    If the ``invert`` argument is True perform the inverse change of basis (from the u'_i to the u_i).

    Args:
        omega (:py:class:`float`) : angular frequency of the pulse, expressed in the
            reciprocal unit of the time variable
        time (:py:class:`array`) : array with the time values
        u (:py:class:`array`) : array with the Bloch vector in the original frame. The
            function assumes that the array has two indices, the first is the component
            and the second one is the time index
        invert (:py:class:`bool`) : if True perform the inverse transformation

    Returns:
        (:py:class:`array`) : array with the transformed Bloch vector

    """
    uprime = np.zeros([3,len(time)])
    if invert: alpha = -1
    else: alpha = 1
    uprime[0,:] = np.cos(omega*time)*u[0,:] + alpha*np.sin(omega*time)*u[1,:]
    uprime[1,:] = -alpha*np.sin(omega*time)*u[0,:] + np.cos(omega*time)*u[1,:]
    uprime[2,:] = u[2,:]
    return uprime
