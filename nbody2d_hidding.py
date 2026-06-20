from __future__ import annotations

import numpy as np
import cft
from dataclasses import dataclass
from scipy.integrate import quad
from cft import Box
from abc import ABC, abstractmethod
from typing import Generic, TypeVar, Callable
from functools import partial

@dataclass
class Cosmology:
    H0 : float
    OmegaM : float
    OmegaL : float

    # ~\~ begin <<lit/index.md|cosmology-methods>>[0]
    @property
    def OmegaK(self):
        return 1 - self.OmegaM - self.OmegaL
    # ~\~ end
    # ~\~ begin <<lit/index.md|cosmology-methods>>[1]
    @property
    def G(self):
        return 3./2 * self.OmegaM * self.H0**2
    # ~\~ end
    # ~\~ begin <<lit/index.md|cosmology-methods>>[2]
    def da(self, a):
        return self.H0 * a * np.sqrt(
                  self.OmegaL \
                + self.OmegaM * a**-3 \
                + self.OmegaK * a**-2)
    # ~\~ end
    # ~\~ begin <<lit/index.md|cosmology-methods>>[3]
    def growing_mode(self, a):
        if isinstance(a, np.ndarray):
            return np.array([self.growing_mode(b) for b in a])
        elif a <= 0.001:
            return a
        else:
            return self.factor * self.adot(a)/a \
                * quad(lambda b: self.adot(b)**(-3), 0.00001, a)[0] + 0.00001

def md_cic(B: Box, X: np.ndarray) -> np.ndarray:
    """Takes a 2*M array of particle positions and returns an array of shape
    `B.shape`. The result is a density field computed by cloud-in-cell method."""
    f  = X - np.floor(X)

    rho = np.zeros(B.shape, dtype='float64')
    rho_, x_, y_ = np.histogram2d(X[:,0]%B.N, X[:,1]%B.N, bins=B.shape,
                        range=[[0, B.N], [0, B.N]],
                        weights=(1 - f[:,0])*(1 - f[:,1]))
    rho += rho_
    rho_, x_, y_ = np.histogram2d((X[:,0]+1)%B.N, X[:,1]%B.N, bins=B.shape,
                        range=[[0, B.N], [0, B.N]],
                        weights=(f[:,0])*(1 - f[:,1]))
    rho += rho_
    rho_, x_, y_ = np.histogram2d(X[:,0]%B.N, (X[:,1]+1)%B.N, bins=B.shape,
                        range=[[0, B.N], [0, B.N]],
                        weights=(1 - f[:,0])*(f[:,1]))
    rho += rho_
    rho_, x_, y_ = np.histogram2d((X[:,0]+1)%B.N, (X[:,1]+1)%B.N, bins=B.shape,
                        range=[[0, B.N], [0, B.N]],
                        weights=(f[:,0])*(f[:,1]))
    rho += rho_

    return rho

class Interp2D:
    "Reasonably fast bilinear interpolation routine"
    def __init__(self, data):
        self.data = data
        self.shape = data.shape

    def __call__(self, x):
        X1 = np.floor(x).astype(int) % self.shape
        X2 = np.ceil(x).astype(int) % self.shape
        xm = x % 1.0
        xn = 1.0 - xm

        f1 = self.data[X1[:,0], X1[:,1]]
        f2 = self.data[X2[:,0], X1[:,1]]
        f3 = self.data[X1[:,0], X2[:,1]]
        f4 = self.data[X2[:,0], X2[:,1]]

        return  f1 * xn[:,0] * xn[:,1] + \
                f2 * xm[:,0] * xn[:,1] + \
                f3 * xn[:,0] * xm[:,1] + \
                f4 * xm[:,0] * xm[:,1]
        
def gradient_2nd_order(F, i):
    return   1./12 * np.roll(F,  2, axis=i) - 2./3  * np.roll(F,  1, axis=i) \
           + 2./3  * np.roll(F, -1, axis=i) - 1./12 * np.roll(F, -2, axis=i)

class VectorABC(ABC):
    @abstractmethod
    def __add__(self, other: Vector) -> Vector:
        raise NotImplementedError

    @abstractmethod
    def __rmul__(self, other: float) -> Vector:
        raise NotImplementedError

VectorABC.register(np.ndarray)
Vector = TypeVar("Vector", bound=VectorABC)
        
@dataclass
class State(Generic[Vector]):
    time : float
    timestep : int
    position : Vector
    momentum : Vector

    # ~\~ begin <<lit/index.md|state-methods>>[0]
    def kick(self, dt: float, h: HamiltonianSystem[Vector]) -> State[Vector]:
        self.momentum += dt * h.momentumEquation(self)
        return self

    def drift(self, dt: float, h: HamiltonianSystem[Vector]) -> State[Vector]:
        self.position += dt * h.positionEquation(self)
        return self

    def wait(self, dt: float) -> State[Vector]:
        self.time += dt
        self.timestep += 1
        return self

class HamiltonianSystem(ABC, Generic[Vector]):
    @abstractmethod
    def positionEquation(self, s: State[Vector]) -> Vector:
        raise NotImplementedError

    @abstractmethod
    def momentumEquation(self, s: State[Vector]) -> Vector:
        raise NotImplementedError

Solver = Callable[[HamiltonianSystem[Vector], State[Vector]], State[Vector]]
Stepper = Callable[[State[Vector]], State[Vector]]
HaltingCondition = Callable[[State[Vector]], bool]

def leap_frog(dt: float, h: HamiltonianSystem[Vector], s: State[Vector]) -> State[Vector]:
    return s.kick(dt, h).wait(dt/2).drift(dt, h).wait(dt/2)

def iterate_step(step: Stepper, halt: HaltingCondition, init: State[Vector]) -> State[Vector]:
    state = init
    while not halt(state):
        state = step(state)
        #fn = 'data/x.{0:05d}.npy'.format(int(round(state.time*1000)))
        #with open(fn, 'wb') as f:
        #    np.save(f, state.position)
        #    np.save(f, state.momentum)
    return state

class PoissonVlasov(HamiltonianSystem[np.ndarray]):
    def __init__(self, box, cosmology, particle_mass):
        self.box = box
        self.cosmology = cosmology
        self.particle_mass = particle_mass
        #self._g = gp.Gnuplot(persist=True)
        #self._g("set cbrange [0.2:50]", "set log cb", "set size square",
        #        "set xrange [0:{0}] ; set yrange [0:{0}]".format(box.N))
        #self._g("set term qt")

    # ~\~ begin <<lit/index.md|position-equation>>[0]
    def positionEquation(self, s: State[np.ndarray]) -> np.ndarray:
        a = s.time
        da = self.cosmology.da(a)
        return s.momentum / (s.time**2 * da)
    # ~\~ end
    # ~\~ begin <<lit/index.md|momentum-equation>>[0]
    def momentumEquation(self, s: State[np.ndarray]) -> np.ndarray:
        a = s.time
        da = self.cosmology.da(a)
        x_grid = s.position / self.box.res
        delta = md_cic(self.box, x_grid) * self.particle_mass - 1.0
        #plt.scatter(s.position[:,0],s.position[:,1],s=1,c='C0',alpha=0.01)
        #self._g(gp.plot_data(gp.array(delta.T+1, "t'' w image")))
        delta_f = np.fft.fftn(delta)
        kernel = cft.Potential()(self.box.K)
        phi = np.fft.ifftn(delta_f * kernel).real * self.cosmology.G / a
        acc_x = Interp2D(gradient_2nd_order(phi, 0))
        acc_y = Interp2D(gradient_2nd_order(phi, 1))
        acc = np.c_[acc_x(x_grid), acc_y(x_grid)] / self.box.res
        return -acc / da
    # ~\~ end
# ~\~ end

# ~\~ begin <<lit/index.md|initialization>>[0]
def a2r(B, X):
    return X.transpose([1,2,0]).reshape([B.N**2, 2])

def r2a(B, x):
    return x.reshape([B.N, B.N, 2]).transpose([2,0,1])

class Zeldovich:
    def __init__(self, B_mass: Box, B_force: Box, cosmology: Cosmology, phi: np.ndarray):
        self.bm = B_mass
        self.bf = B_force
        self.cosmology  = cosmology
        self.u = np.array([-gradient_2nd_order(phi, 0),
                           -gradient_2nd_order(phi, 1)]) / self.bm.res

    def state(self, a_init: float) -> State[np.ndarray]:
        X = a2r(self.bm, np.indices(self.bm.shape) * self.bm.res + a_init * self.u)
        P = a2r(self.bm, a_init * self.u)
        return State(time=a_init, timestep=0, position=X, momentum=P)

    @property
    def particle_mass(self):
        return (self.bf.N / self.bm.N)**self.bm.dim
# ~\~ end
