
# Table of Contents

1.  [Background](#org28a10a7)
2.  [The Toolkit](#org8b6a8f9)
3.  [Current State](#orga3ff997)
4.  [Recommended First Steps](#org4ad2eb9)
5.  [Feedback](#org8f139e4)


<a id="org28a10a7"></a>

# Background

The simulation of sequences or sequence blocks is a recurrent task in magnetic resonance imaging (MRI).
Nonetheless, the chosen method is often situational and strongly depends on the underlying assumptions.
While Bloch(-Torrey) equations are completely sufficient to handle the microscopic scale,
their use on the larger voxel scale is rather cumbersome and provides little insight. 
Here, **extended phase graphs (EPG)** are much better suited, as they allow to predict
the appearance of echoes and their amplitudes. Nonetheless, they have some limitations as well, 
e.g., when it comes to

-   quantify susceptibility effects
-   determine the proper strength of crusher gradients
-   handle non-periodic sequences

Origin of these difficulties is the somewhat fuzzy definition of EPG states as such, which
are usually motivated by an early voxel scale Fourier integral.

As shown in the [documentation](doc/configuration_model.pdf), these technical limitations can be overcome in 
the microscopic **configuration model (CM)**, which is closely related to the EPG formalism,
but postpones the transition to the voxel scale, until it is actually needed.


<a id="org8b6a8f9"></a>

# The Toolkit

`CoMoTk` is a Matlab simulation tool, which implements the configuration model. It can be used for a 
wide range of situations, e.g. simulation of

-   idealized sequences (instantaneous RF pulses, neglecting actual gradients, &#x2026;)
-   complex sequences (as played out on actual scanners, including finite and selective RF pulses)
-   sequence block, e.g. selectice RF pulses (pulse design)
-   susceptibility and/or diffusion effects
-   &#x2026;

For installation, just add the `matlab/` folder with subdirectories to your Matlab path.

**Important**: A Matlab version of **R2016b** or later is required, since implicit expansion is used a lot.


<a id="orga3ff997"></a>

# Current State

`CoMoTk` is currently work in progress and should be functional to the extent as shown in the `examples/`
folder. The following table gives a brief overview on the actual state.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">What</th>
<th scope="col" class="org-left">Implemented</th>
<th scope="col" class="org-left">Tested</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">Bloch equations</td>
<td class="org-left">yes</td>
<td class="org-left">yes</td>
</tr>


<tr>
<td class="org-left">Diffusion</td>
<td class="org-left">yes</td>
<td class="org-left">yes</td>
</tr>


<tr>
<td class="org-left">Derivatives</td>
<td class="org-left">yes</td>
<td class="org-left">yes</td>
</tr>
</tbody>
</table>

It will be updated.


<a id="org4ad2eb9"></a>

# Recommended First Steps

-   The [script](doc/configuration_model.pdf) in the `doc/` folder is crucial to understand the theoretical background on the configuration model and also gives a brief overview on how to use `CoMoTk`.
-   The scripts in the folders `examples/` and `test/` show typical usage scanarios and should serve as good starting points.


<a id="org8f139e4"></a>

# Feedback

Comments? Wishes? Bugs? - Please let me know via the [issue tracker](https://github.com/cganter/CoMoTk/issues) or write an email to
[comotk.rad.med@tum.de](mailto:comotk.rad.med@tum.de).

With respect to bugs: Most helpful are minimum size example scripts, which generate the error.

