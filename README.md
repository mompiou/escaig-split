# Escaig split

A script to compute the shear stresses and orientation factor in the Friedel-Escaig model of cross-slip in fcc crystals in tension (reverse the values for compression).

Run ```escaig-split.py``` to compute the map and ```escaig-plot.py``` to draw. Use ```escaig-split-1point.py``` to get value for 1 orientation and one slip system.

```escaig-plot.py``` draws:

- the Schmid factor (SF) for the primary slip plane, i.e. tau/sigma:

- the difference in the SF of partial dislocations in the primary plane (tau'_d/sigma), for intrinsic fault (reverse the value otherwise)

- the diffrence in the SF of partial dislocations in the cross-slip plane (tau_d/sigma), for intrinsic fault (reverse otherwise)

- the orientation factor prop to (tau_d -tau_d')/tau

![](img.png)

ref: Bonneville, J., Escaig B., Acta Metallurgica, 27, 1477-1486, 1979.
