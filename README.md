# Escaig split

A script to compute the shear stresses and orientation factor in the Friedel-Escaig model of cross-slip in fcc crystals in tension (reverse the values for compression).

```escaig-split.py``` draws:

- the Schmid factor (SF) for the primary slip plane, i.e. tau/sigma:

- the difference in the SF of partial dislocations in the primary plane (tau'_d/sigma)

- the diffrence in the SF of partial dislocations in the cross-slip plane (tau_d/sigma)

- the orientation factor prop to (-2/3 tau_d +tau_d')/tau

![](img.png)

ref: Bonneville, J., Esacaig B., Acta Metallurgica, 27, 1477-1486, 1979.
