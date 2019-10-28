
# Table of Contents

1.  [Background](#org6d452be)
2.  [The Toolkit](#orge96d733)
    1.  [Scope](#org682792c)
3.  [Usage](#org98e679e)
    1.  [Initialize `CoMoTk`](#org2374442)
    2.  [RF pulse](#org4ffab4e)
    3.  [Time interval](#orgb15d64f)
    4.  [Get results](#org30650e4)
4.  [Current State](#orgd6e8ef9)
5.  [Recommended First Steps](#org41b1e6c)
6.  [Feedback](#org3b21848)


<a id="org6d452be"></a>

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


<a id="orge96d733"></a>

# The Toolkit


<a id="org682792c"></a>

## Scope

`CoMoTk` is a Matlab simulation tool, which implements the configuration model. It can be used for a 
wide range of situations, e.g. simulation of

-   idealized sequences (instantaneous RF pulses, neglecting actual gradients, &#x2026;)
-   complex sequences (as played out on actual scanners, including finite and selective RF pulses)
-   sequence block, e.g. selectice RF pulses (pulse design)
-   susceptibility and/or diffusion effects
-   &#x2026;

For installation, just add the `matlab/` folder with subdirectories to your Matlab path.

**Important**: A Matlab version of **R2016b** or later is required, since implicit expansion is used a lot.


<a id="org98e679e"></a>

# Usage

To familiarize with the toolkit, the scripts in the `test` and `examples` folders are recommended as a good
starting point. Following, a brief overview of the basic functionality.


<a id="org2374442"></a>

## Initialize `CoMoTk`

First, we need to create an instance of `CoMoTk`:

    cm = CoMoTk;

Next, we have to setup a few mandatory tissue parameters. In the simpliest case
of a 1-peak model without diffusion, this can look like this:

    cm.R1 = 0.01;       % longitudinal relaxation rate
    cm.R2 = 0.1;        % transverse relaxation rate
    cm.D = 0;           % apparent diffusion coefficient

It is also possible, to specify more complex n-peak models with variable relative weighting
(\( \propto \) proton density), chemical shift and diffusivity. Here, a 2-peak example:

    cm.R1 = [ 0.01, 0.02 ];
    cm.R2 = [ 0.1, 0.2 ];
    cm.D = [ 3, 1.5 ];
    cm.w = [ 0.6, 0.4 ];    % relative weight (does not need to add to 1)
    cm.dom = [ 0, -0.1 ];   % chemical shift (angular frequency)

Setting `w` and `dom` is optional. If not set, they default to 1 and 0, respectively.

Another optional variable is the relative \( B_1^+ \) field (default = 1):

    cm.B1 = 0.8;

There are further options, for which it is possible to modify the default settings. The most important
one is the desired accuracy. A nonzero value of `epsilon` is set, if the number of stored
configurations needs to be restricted, e.g. due to memory or performance restrictions.

    options = cm.options;    % get default options
    
    options.alloc_d = 3;     % allocated number of dimensions
    options.alloc_n = 10000; % allocated number of configurations
    options.epsilon = 0;     % for maximal accuracy (enough memory needed)
    options.verbose = true;  % for more output
    options.debug = true;    % for debugging purposes (look into CoMoTk.m)
    
    cm.options = options;    % activate new options

Now we have to specify the initial configuration vector, corresponding to \( \bm{n} = 0 \). 
It is supplied as a real vector (row or column) in the usual convention \( \left( m_x, m_y, m_z \right) \):

    cm.init_configuration ( [ 0; 0; 1 ] );  % longitudinal magnetization

In addition to \( m_\nu^{\left(\bm{n}\right)} \), knowledge of certain partial derivatives 
\( \partial m_\nu^{\left(\bm{n}\right)}/\partial \xi \) is sometimes desired as well, 
e.g. for numerical optimization.
This is supported in `CoMoTk` for \( \xi \in \left\{ R_1, R_2, D, B_1, \alpha_\mu, \varphi_\mu, \tau_\mu,
   \bm{p}_\mu, \bm{s}_\mu, \bm{S}_\mu \right\} \). The following command (a full example, involving every
possible variable) is used, to inform `CoMoTk`, which derivatives need to be calculated:

      cm.set_derivatives ( ...
    'R1', [ 1, 3 ], ...           % tissue handles (n-peak model)
    'R2', 2, ...                  % for a 1-peak model the handle (= 1) 
    'D', [ 2, 3 ], ...            % must be supplied as well
    'B1', ...                     % only B1 has no handle
    'FlipAngle', [ 2, 3 ], ...    % non-tissue handles can be chosen
    'Phase', [ 2, 1004, 12 ], ... % arbitrarily
    'tau', [ 3, 4, 6 ], ...       % time interval (= \mu)
    'p', [ 1, 2, 4 ], ...         % gradient moment (see below)
    's', [ 4, 7 ] ...             % gradient shape (see below)
    );

Of course, only the desired subset of parameter/handle combinations needs to be supplied.
If the command is not given, no derivatives are calculated by default.

After having initialized everything, we proceed with the actual sequence. Here, we apply instantaneous 
RF pulses and time intervals, typically in an alternate fashion.


<a id="org4ffab4e"></a>

## RF pulse

Calling an RF pulse is as simple as:

    cm.RF( flip, phase );  % RF pulse with flip angle and phase [rad]

If partial derivatives with respect to flip angle(s) and/or phase(s) have been set before, 
the associated handles can be additionally supplied according to the following format:

    cm.RF( flip, phase, 'FlipAngle', flip_handle, 'Phase', phase_handle );


<a id="orgb15d64f"></a>

## Time interval

Whether the time derivatives are needed or not, executing the time interval always needs specification
of the (otherwise arbitrary) index \( \texttt{mu} := \mu \), which is defined as
in section XXX. At least in the first call of each distinct interval, the duration 
\( \texttt{tau} := \tau_\mu \) must be supplied also:

    cm.time( mu, 'tau', tau );

If we also require the gradient moment \( \texttt{p} := \bm{p}_\mu \) for the simulations 
(diffusion effects) or the results (e.g. slice profile), we have to add this parameter as well:

    cm.time( mu, 'tau', tau, 'p', p );

Called like this, a constant gradient shape according to \eqref{eq:pnulin} and  \eqref{eq:Fnconstant} is assumed.

For arbitrary shapes, the variable \( \texttt{s} := 
   \left( \bm{s}_\nu,\operatorname{tr}\left(\bm{S}_\nu\right) \right) \) must be added too:

    cm.time( mu, 'tau', tau, 'p', p, 's', s );

In later calls, only the handle is required:

    cm.time( mu );

Note, however, that this shorthand form is not safe for nonzero `epsilon`, since the dimension `mu` could
have been eliminated (together with any knowledge about `tau` and `p`) by a previous call of `meltdown()`.
In case of doubt, always specify all parameters.


<a id="org30650e4"></a>

## Get results

After any RF pulse or time interval, we can obtain the sum \eqref{eq:ansatz} like this:

    res = cm.sum( param );

`param` is a structure with optional fields:

-   **`omega` =:** Local angular off-resonance frequency \( \omega\left(\bm{x}\right) \)
-   **`x` =:** Position \( \bm{x} \) according to the definition \eqref{eq:accphasedyn}
-   **`b_n` =:** Subset from the set of stored configuration orders, \( S_\nu^\pm \), to be included in the result.
-   **`w_n` =:** explicit weighting factors (`length( w_n ) = sum( b_n )`)

For unset fields, the following defaults are assumed:

-   `omega` \( = 0 \)
-   `x` \( = 0 \)
-   `b_n` \( = \) `cm.b_n` (= whole sum in Eq. \eqref{eq:ansatz})
-   `w_n` \( = 1 \)

The result is returned separately as transverse (`res.xy`) and longitudinal (`res.z`) component.
Calculated derivatives with respect to `X` are returned as `res.dm_dX.xy` and `res.dm_dX.z`, where
\( \texttt{X} \) is any member of the set \( \left\{ \texttt{R1}, \texttt{R2}, \texttt{D}, \texttt{B1},
   \texttt{FlipAngle}, \texttt{Phase}, \texttt{tau}, \texttt{p}, \texttt{s} \right\} \). We have
`size(res.dm_dX) = [a,b]`, where `a = 1`, except for `X = p` with `a = 3` and `X = s` with `a = 4`.
The second dimension `b` is just the number of derivatives to be calculated for each parameter as defined in
`set_derivatives()`.

For a restriction of the summation, the `find` method can be used to generate `b_n`:

    b_n = cm.find( mu, n );

If `mu` and `n` are scalars, we have 
\( \texttt{b\_n} = \left\{ \bm{n}\in S_\nu^\pm: n_{\texttt{mu}} = \texttt{n} \right\} \).
If no matching configurations are found, \( \texttt{b\_n = []} \) is returned.

`mu` and `n` can also be 1d-arrays of equal length. 
This can be used as a shorthand for a combination with the `|` operator. For example,
we obtain for arrays of length 3 a result equivalent to

      b_n = ...
    cm.find( mu( 1 ), n( 1 ) ) | ...
    cm.find( mu( 2 ), n( 2 ) ) | ...
    cm.find( mu( 3 ), n( 3 ) );

If we want to single out a specific configuration \( \bm{n} \in S_\nu^\pm \), we have to use
the `&` operator explicitly. For example, we get for \( d = 3 \):

      b_n = ...
    cm.find( cm.mu( 1 ), n( 1 ) ) & ...
    cm.find( cm.mu( 2 ), n( 2 ) ) & ...
    cm.find( cm.mu( 3 ), n( 3 ) ) );

It is also possible, to restrict the mask `b_n` directly. For example, to extract all stored configurations,
which satisfy \( \bm{p}_{\bm{n}} = \bm{0} \), cf. Eq. \eqref{eq:zeroselect}, one could use

    b_n = cm.b_n & reshape( ~any( cm.p_n ), size( cm.b_n ) );


<a id="orgd6e8ef9"></a>

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


<a id="org41b1e6c"></a>

# Recommended First Steps

-   The [script](doc/configuration_model.pdf) in the `doc/` folder is crucial to understand the theoretical background on the configuration model and also gives a brief overview on how to use `CoMoTk`.
-   The scripts in the folders `examples/` and `test/` show typical usage scanarios and should serve as good starting points.


<a id="org3b21848"></a>

# Feedback

Comments? Wishes? Bugs? - Please let me know via the [issue tracker](https://github.com/cganter/CoMoTk/issues) or write an email to
[comotk.rad.med@tum.de](mailto:comotk.rad.med@tum.de).

With respect to bugs: Most helpful are minimum size example scripts, which generate the error.

