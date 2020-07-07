
# Table of Contents

1.  [Background](#orgb4a0678)
2.  [Installation](#org2b9fe90)
3.  [Small Example](#org2525229)
4.  [Current State](#orgea046bb)
5.  [Feedback](#org869e72b)


<a id="orgb4a0678"></a>

# Background

The simulation of sequences or sequence blocks is a recurrent task in magnetic resonance imaging (MRI).
How the calculation is actually performed, is highly variable and requires a decision about the desired degree of realism.

The following non-exhaustive list shows a few possible options and alternatives:

-   RF pulses 
    -   constant flip angles vs. slice profile
    -   instantaneous vs. finite duration
    -   RF pulse design
-   gradients
    -   ideal vs. real crusher gradients (suppression of unwanted echoes)
    -   gradient shapes
    -   eddy currents
-   tissues
    -   pure vs. mixture (e.g. water/fat)
    -   bulk off-resonance
    -   susceptibility effects
    -   diffusion effects (isotropic, directional)
    -   bulk motion
    -   magnetization exchange/transfer
-   sequence
    -   periodic vs. non-periodic
    -   transient phase vs. steady-state

`CoMoTk` is a general purpose simulation tool (written in Matlab), designed to handle essentially all the aforementioned cases. 
It is based upon a generalized version of the so-called configuration model, which will be described elsewhere.


<a id="org2b9fe90"></a>

# Installation

Just add the `matlab/` folder with subdirectories to your Matlab path.

**Important**: A Matlab version of **R2016b** or later is required, since implicit expansion is used a lot.


<a id="org2525229"></a>

# Small Example

The following example should give a first impression of how to work with the toolkit, while the [User Guide](doc/CoMoTk_UserGuide.pdf) gives a complete overview of the available options. It is also recommended to look at the scripts in the `test` and `examples` folders, which cover typical application scenarios and may serve as templates for own projects. To make full use of the toolkit, it is important to understand the theory behind the configuration model, which will be presented in an article (link will be given, as soon as available).

    % Example: sample FID of SSFP during transient phase 
    
    % initialize class
    cm = CoMoTk;      
    
    % set tissue parameters
    cm.R1 = R1;                     % longitudinal relaxation rate
    cm.R2 = R2;                     % transverse relaxation rate
    cm.D = 0;                       % we neglect diffusion in this example
    
    % set RF pulse parameters
    rf.FlipAngle = flip_angle;      % flip angle [rad]
    rf.Phase = phase;               % phase [rad]
    
    % set time interval from RF pulse to echo
    te.mu = 1;                      % unique index of time interval
    te.tau = TE;                    % duration
    
    % set time interval from echo to RF pulse (includes crusher)
    crusher.mu = 2;                 % unique index of time interval
    crusher.tau = TR - TE;          % duration
    
    % desired number of samples
    n_TR = 100;
    
    % allocate space for results
    m_xy = zeros( n_TR, 1 );
    
    % set initial state (here: pure longitudinal magnetization)
    cm.init_configuration ( [ 0; 0; 1 ] );  
    
    % execute sequence loop
    for i = 1 : n_TR                
    
      % RF pulse
      cm.RF( rf );                  
    
      % time to echo
      cm.time( te );               
    
      % FID == all configurations with no (0) dephasing by crusher
      fid.b_n = cm.find( crusher, 0 );
    
      % extract result
      res = cm.sum( fid );
      m_xy( i ) = res.xy;
    
      % time to next RF pulse (includes crusher)
      cm.time( crusher );           
    
    end
    % done


<a id="orgea046bb"></a>

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


<a id="org869e72b"></a>

# Feedback

Please use the [issue tracker](https://github.com/cganter/CoMoTk/issues) or write an email to [comotk.rad.med@tum.de](mailto:comotk.rad.med@tum.de).

Reporting bugs: Most helpful are small example scripts, which generate the error.

