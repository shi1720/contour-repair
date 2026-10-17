import itertools
import unittest
import numpy as np
from src.contour_recovery import recover,dense_rule,error_norm,spectrum,observed_error_bound


class RecoveryTests(unittest.TestCase):
    def test_all_256_masks_against_independent_dense_solve(self):
        for mask in range(256):
            e=[j for j in range(8) if mask>>j&1]
            r=recover(8,.7,e)
            d=dense_rule(8,.7,e)
            np.testing.assert_allclose(r.weights,d,rtol=2e-12,atol=2e-13)
            self.assertLess(r.normal_residual,1e-12)
            base=2*.7**8/(1+.7**8)
            self.assertAlmostEqual(r.error**2,base+r.loss,places=12)
            for j,g in r.gains.items():
                self.assertAlmostEqual(g,r.error**2-recover(8,.7,[k for k in e if k!=j]).error**2,places=12)

    def test_error_against_direct_fourier_series(self):
        rng=np.random.default_rng(451)
        n=17;q=.8;theta=2*np.pi*np.arange(n)/n
        w=rng.normal(size=n)
        direct=(1-w.sum())**2
        for k in range(1,400):
            direct+=2*q**k*abs(w@np.exp(1j*k*theta))**2
        self.assertAlmostEqual(error_norm(w,q)**2,direct,places=10)

    def test_vector_noisy_bound(self):
        rng=np.random.default_rng(764)
        n=24;q=.75;r=recover(n,q,[0,1,8,11]);theta=2*np.pi*np.arange(n)/n
        modes=np.arange(-40,41)
        a=(rng.normal(size=(81,3))+1j*rng.normal(size=(81,3)))*q**(abs(modes)[:,None]/2)
        a/=np.sqrt(np.sum(abs(a)**2*q**(-abs(modes))[:,None]))
        y=np.exp(1j*theta[:,None]*modes)@a
        noise=(rng.normal(size=y.shape)+1j*rng.normal(size=y.shape))*1e-5
        y+=noise
        bound=observed_error_bound(r.weights,q,1,np.linalg.norm(noise,axis=1))
        self.assertLessEqual(np.linalg.norm(r.weights@y-a[40]),bound)

    def test_invalid(self):
        for e in [[-1],[8],[1,1],[1.2]]:
            with self.assertRaises(ValueError):recover(8,.7,e)
        with self.assertRaises(ValueError):recover(4096,.5,[])

    def test_extremizer_attains_bound(self):
        n=12;q=.6;e=[0,1,5];r=recover(n,q,e)
        theta=2*np.pi*np.arange(n)/n
        k=np.arange(-180,181)
        # Residual representer: coefficients q^|k| (delta_k0 - sum_j w_j exp(-ik theta_j)).
        a=q**abs(k)*((k==0)-(r.weights@np.exp(-1j*theta[:,None]*k)))/r.error
        values=np.exp(1j*theta[:,None]*k)@a
        available=[j for j in range(n) if j not in e]
        self.assertLess(np.max(abs(values[available])),1e-12)
        self.assertAlmostEqual(abs(a[180]-r.weights@values),r.error,places=12)
        self.assertAlmostEqual(float(np.sum(abs(a)**2*q**(-abs(k)))),1.,places=12)


if __name__=='__main__':unittest.main()
