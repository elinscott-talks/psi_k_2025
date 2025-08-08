import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class ECurve:

   def __init__(self, label, integer_energies, curvature = None):
      self.label = label
      self.integer_energies = integer_energies
      self.integer_x = np.arange(len(self.integer_energies))

      # curvature can be a float or a list; make sure self.curvature is a list
      if isinstance(curvature, float) or isinstance(curvature, int):
         self.curvature = [curvature for _ in range(len(integer_energies) - 1)]
      else:
         self.curvature = curvature

      self.curve = None
      self.extrapolations = {}
      self.Delta_E_arrow = None
      self.Delta_E_text = None
      self.error_arrow = None
      self.error_text = None
      self.extrap_arrow = None
      self.extrap_text = None
      self.x_variable = None
      

   _resolution = 100

   @property
   def x(self):
      if self.curvature is None:
         return self.integer_x
      else:
         l = max(self.integer_x)
         return np.linspace(0, l, l*self._resolution + 1)

   @property
   def y(self):
      if self.curvature is None:
         return self.integer_energies
      else:
         # Generate parabola going through the integer values 
         # (x0, y0) and (x1, y1) with the curvature c
         # y = 1/2 c (x - x0)(x - x1) 
         #      + (y1 - y0)(x - x0)/(x1 - x0)
         #       + y0
         y = np.empty(len(self.x))
         for i, (x0, x1, y0, y1, c) in enumerate(zip(self.integer_x[:-1],
            self.integer_x[1:], self.integer_energies[:-1],
            self.integer_energies[1:], self.curvature)):
            x_window = self.x[i*self._resolution: (i + 1)*self._resolution + 1]
            y[i*self._resolution: (i + 1)*self._resolution + 1] = \
                1/2*c*(x_window - x0)*(x_window - x1) \
                + (y1 - y0)*(x_window - x0)/(x1 - x0) \
                + y0
         return y

   @property
   def color(self):
      if self.curve is None:
         return None
      else:
         return self.curve.get_color()

   def plot(self, ax, **kwargs):
      if self.curvature is None:
         marker = 'o'
      else:
         marker = ''
      curve = ax.plot(self.x, self.y, marker=marker, label=self.label, **kwargs)
      self.curve = curve[0]
      ax.set_xlabel(f'${self.x_variable}$')
      ax.set_ylabel(f'$E({self.x_variable})$', rotation=90)
      sns.despine()
      ax.set_yticks([])
      ax.set_xticks(self.integer_x)
      ax.legend()

   def extrapolate(self, i_orb, forward=True, length=1):
      x0 = self.integer_x[i_orb]
      y0 = self.integer_energies[i_orb]
      if forward:
         x1 = self.integer_x[i_orb + 1]
         y1 = self.integer_energies[i_orb + 1]
      else:
         x1 = self.integer_x[i_orb - 1]
         y1 = self.integer_energies[i_orb - 1]
      
      if self.curvature is None:
         c = 0
      else:
         c = self.curvature[i_orb]
      # Gradient of y at x = x0 is given by 1/2 c (x0 - x1) + (y1 - y0) / (x1 - x0)
      m = 1/2*c*(x0 - x1) + (y1 - y0) / (x1 - x0)

      # Scale
      x1 = x0 + (x1 - x0)*length
      y1 = y0 + (y1 - y0)*length

      if forward:
         y1_extrap = y0 + m * length
      else:
         y1_extrap = y0 - m * length

      return x0, y0, x1, y1, y1_extrap

   def plot_gradient(self, ax, i_orb, forward=True, label=None, length=0.2, label_dx=0.0, label_dy=0.0, va='center', ha='center'):
      x0, y0, x1, _, y1_extrap = self.extrapolate(i_orb, forward, length)

      line = ax.annotate('',xy=(x0,y0), xytext=(x1,y1_extrap), xycoords='data', textcoords='data', color=self.color, 
                  arrowprops={'arrowstyle':'<-', 'color': self.color})

      self.extrapolations[(i_orb, forward)] = line
      
      if label:
         text = ax.text((x0 + x1) / 2 + label_dx, (y0 + y1_extrap) / 2 + label_dy, label, color=self.color, va=va, ha=ha, fontsize=12)


   def plot_extrapolation(self, ax, i_orb, forward=True, text_label=None, label_error=None, label_extrap=None, **kwargs):

      x0, y0, x1, y1, y1_extrap = self.extrapolate(i_orb, forward)

      plotkwargs = {'ls': '--'}
      plotkwargs.update(kwargs)
      line = ax.plot([x0, x1], [y0, y1_extrap], color=self.color, alpha=0.5, **plotkwargs)

      self.extrapolations[(i_orb, forward)] = line[0]

      if label_error or label_extrap:
         if forward:
            dx = -0.05
            ha = 'right'
         else:
            dx = 0.05
            ha = 'left'

         if label_extrap:
            arrow = ax.annotate('',xy=(x1,y0), xytext=(x1,y1_extrap), xycoords='data', textcoords='data', color=self.color, 
                        arrowprops={'arrowstyle':'<->', 'color': self.color})
            text = ax.text(x1 + dx, (y0 + y1_extrap) / 2, label_extrap, color=self.color, va='center', ha=ha, fontsize=12)
            self.extrap_arrow = arrow
            self.extrap_text = text

         if label_error:
            arrow = ax.annotate('',xy=(x1,y1), xytext=(x1,y1_extrap), xycoords='data', textcoords='data', color=self.color, 
                        arrowprops={'arrowstyle':'<->', 'color': self.color})
            text = ax.text(x1 + dx, (y1 + y1_extrap) / 2, label_error, color=self.color, va='center', ha=ha)
            self.error_arrow = arrow
            self.error_text = text

      if text_label:
         ax.text(x1, y1, text_label, color=self.color)

   def label_Delta_E(self, ax, i_orb, delta_orb=-1, label=r'$\Delta E$', **kwargs):
      # Add a label showing the change in total energy
      x0 = self.integer_x[i_orb]
      x1 = self.integer_x[i_orb + delta_orb]
      y0 = self.integer_energies[i_orb]
      y1 = self.integer_energies[i_orb + delta_orb]

      if delta_orb > 0:
         dx = 0.05
         ha = 'left'
      else:
         dx = -0.05
         ha = 'right'
      # [arrow] = ax.plot([x1, x1, x0], [y1, y0, y0], color=self.color)
      arrow = ax.annotate('',xy=(x1+dx,y0), xytext=(x1+dx,y1), xycoords='data', textcoords='data', color=self.color, 
                  arrowprops={'arrowstyle':'<->', 'color': self.color})
      text = ax.text(x1 + 2*dx, (y0 + y1) / 2, label, color=self.color, va='center', ha=ha)

      self.Delta_E_arrow = arrow
      self.Delta_E_text = text

   def label_Delta_e(self, i_orb, **kwargs):
      # Add a label showing the eigen-energy
      return

   def remove_annotations(self):
      objs_to_remove = list(self.extrapolations.values()) + [self.Delta_E_arrow, self.Delta_E_text, self.error_arrow, self.error_text, self.extrap_arrow, self.extrap_text]
      for obj in objs_to_remove:
         if obj is not None:
            obj.remove()

def save(ax, name):
   ax.set_xlim([0.5, 3.5])
   ax.set_ylim([-0.5, 3])
   ax.set_xticks(range(1,4))
   ax.set_xticklabels(['$N-1$','$N$','$N+1$'])
   ax.set_xlabel('total number of electrons')
   plt.savefig(name)

class ENCurve(ECurve):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.x_variable = 'N'

class EMCurve(ECurve):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.x_variable = r'\mu'

if __name__ == '__main__':

   sns.set_style('white')
   sns.set_context('talk')

   exact_e = np.array([5, 1.5, 0, 0.5, 1.8])

   i_orb = 2

   exact = ENCurve('exact', exact_e)
   hf = ENCurve('Hartree-Fock', exact_e + 0.3*np.random.rand(len(exact_e)), curvature = -2)
   sl = ENCurve('semi-local', exact_e - 0.3*np.random.rand(len(exact_e)), curvature = 3)

   f, ax = plt.subplots(1,1)
   plt.subplots_adjust(bottom=0.2)

   # Plot of exact piecewise-linear behaviour only
   exact.plot(ax, ls='--')
   save(ax, 'fig_en_curve_exact.pdf')

   # Add annotations
   exact.label_Delta_E(ax, i_orb)
   exact_epsilon_label = ax.text(i_orb - 0.2, exact_e[i_orb] + 0.5, r'$\left.\frac{dE}{dN}\right|_{N = N^-} = \varepsilon_\mathrm{HO}$', color=exact.color, rotation=0)
   save(ax, 'fig_en_curve_exact_annotated.pdf')

   exact_epsilon_label.remove()
   exact.remove_annotations()

   # Plot of Hartree-Fock
   hf.plot(ax)
   save(ax, 'fig_en_curve_hf.pdf')

   # With annotations
   hf.plot_extrapolation(ax, i_orb, forward=False, label_error='error')
   hf_epsilon_label = ax.text(i_orb - 0.2, hf.integer_energies[i_orb] + 0.8, r'$\varepsilon_\mathrm{HO} = \left.\frac{dE}{dN}\right|_{N = N^-}$', color=hf.color)
   hf.label_Delta_E(ax, i_orb)
   save(ax, 'fig_en_curve_hf_annotated.pdf')

   # With semi-local DFT, keeping annotations
   sl.plot(ax)
   save(ax, 'fig_en_curve_with_all_hf_annotated.pdf')

   hf.remove_annotations()
   hf_epsilon_label.remove()

   # With semi-local DFT
   sl.plot(ax)
   save(ax, 'fig_en_curve_with_all.pdf')

   # Stepping through the Koopman's correction
   hf.curve.remove()

   # Step 0: the exact solution DFT
   f, ax = plt.subplots(1,1)
   plt.subplots_adjust(bottom=0.2)
   exact.plot(ax, ls='--')
   ax._get_lines.get_next_color()
   sl.plot(ax)
   save(ax, 'fig_en_curve_koopmans_step0.pdf')

   # Step 1: the starting point
   f_i = 2.58
   point = ax.plot(f_i, sl.y[[abs(x - f_i) < 0.00001 for x in sl.x]], 'o', markersize=15)
   point_text = ax.text(f_i, sl.y[[abs(x - f_i) < 0.00001 for x in sl.x]], '1', color='w', va='center', ha='center', fontsize=8)
   point = point[0]
   save(ax, 'fig_en_curve_koopmans_step1.pdf')

   # Step 2: subtracting the curve
   # point.remove()
   point = ax.plot(np.floor(f_i), sl.y[[abs(x - np.floor(f_i)) < 0.00001 for x in sl.x]], 'o', color=point.get_color(), markersize=15)
   point = point[0]
   point_text = ax.text(np.floor(f_i), sl.y[[abs(x - np.floor(f_i)) < 0.00001 for x in sl.x]], '+2', color='w', va='center', ha='center', fontsize=8)
   mask = [x < f_i and x >= np.floor(f_i) for x in sl.x]
   ax.plot(sl.x[mask], sl.y[mask], color=point.get_color())
   save(ax, 'fig_en_curve_koopmans_step2.pdf')

   # Step 3: adding the piecewise-linear behaviour
   fixed_sl = ENCurve(None, sl.integer_energies, [0 for _ in range(len(sl.integer_energies) - 1)])
   # point.remove()
   point = ax.plot(f_i, fixed_sl.y[[abs(x - f_i) < 0.00001 for x in fixed_sl.x]], 'o', color=point.get_color(), markersize=15)
   point_text = ax.text(f_i, fixed_sl.y[[abs(x - f_i) < 0.00001 for x in fixed_sl.x]], '+3', color='w', va='center', ha='center', fontsize=8)
   point = point[0]
   ax.plot(fixed_sl.x[mask], fixed_sl.y[mask], color=point.get_color())
   save(ax, 'fig_en_curve_koopmans_step3.pdf')
