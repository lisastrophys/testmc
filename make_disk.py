"""Write grid.inp and dust_density.inp for a disk with its midplane at x = 0
(the observer at -x sees it face-on), radius measured in the y-z plane:

    Sigma(r) = sigma0 * (r/AU)**-1 * exp(-r/r_c)          [g cm-2], for r > r_inn
    h(r)     = h_asp * r * (r/AU)**psi                     [AU]
    rho      = Sigma / (sqrt(2 pi) h) * exp(-x**2 / (2 h**2))

The vertical (x) axis gets fine cells near the midplane to resolve the scale height.
"""
import numpy as np

AU = 1.495978707e13          # cm
sigma0, r_c = 1e3, 20.0      # g cm-2 at 1 AU (dust), taper radius [AU]
h_asp, psi  = 0.1, 0.1
r_inn       = 1.0            # AU

# grid edges [AU]
def vertical_edges(dx_fine=0.02, x_fine=1.0, x_max=50.5, n_out=30):
    fine = np.arange(0, x_fine + 0.5*dx_fine, dx_fine)            # 0 .. x_fine
    out  = np.geomspace(x_fine, x_max, n_out + 1)[1:]               # x_fine .. x_max
    half = np.concatenate([fine, out])
    return np.concatenate([-half[::-1], half[1:]])

xb = vertical_edges()
yb = np.arange(-50.5, 51.5, 1.0)
zb = np.arange(-50.5, 51.5, 1.0)

x = 0.5*(xb[1:] + xb[:-1]); y = 0.5*(yb[1:] + yb[:-1]); z = 0.5*(zb[1:] + zb[:-1])
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')                       # rho[x, y, z]
r = np.sqrt(Y**2 + Z**2)
with np.errstate(divide='ignore', invalid='ignore'):
    sigma = np.where(r > r_inn, sigma0 / r * np.exp(-r / r_c), 0.0)
    h = h_asp * r * r**psi
    rho = np.where(r > r_inn, sigma / (np.sqrt(2*np.pi) * h * AU) * np.exp(-X**2 / (2*h**2)), 0.0)

with open('grid.inp', 'w') as f:
    f.write('# cell edges in AU: x (vertical), y, z\n')
    for e in (xb, yb, zb):
        f.write(' '.join('%.8g' % v for v in e) + '\n')
# layout of dust_temperature3d.out: blocks of constant z, lines from high to low y, x along a line
np.savetxt('dust_density.inp', rho[:, ::-1, :].transpose(2, 1, 0).reshape(-1, rho.shape[0]),
           fmt='%.5e', header='dust mass density [g cm-3]')
print('cells: nx=%d ny=%d nz=%d' % rho.shape)
