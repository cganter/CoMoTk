

# Background

The continuous configuration model (CCM) (link follows) constitutes a representation
of arbitrary magnetic resonance imaging (MRI) sequences and tissues, which allows for
a rigorous and transparent treatment of signal localization by
selective excitation and/or spatial encoding.

The configuration model toolkit (CoMoTk), specifically the Matlab class `CoMo`, implements the CCM
for Bloch-Torrey and Bloch-McConnell equations.
The `examples/` folder contains all scripts needed to reproduce the results from the article and
the associated supporting information.


# Installation

Just add the `matlab/` folder with subdirectories to your Matlab path.

**Important**: A Matlab version of **R2016b** or later is required, since implicit expansion is used a lot.


# Small Example

The following example should give a first impression of how to work with the toolkit.
Further information is provided in the documented code and the scripts in the `examples/` folder cover typical application 
scenarios. To fully benefit from the toolkit, it is indispensable to take a closer look at the linked article though.

    % Example: sample FID of SSFP during transient phase 
    % sequence and tissue parameters
    TR = 1;                         % repetition time [ms]
    TE = 0.1;                       % echo time [ms]
    flip_angle = pi / 2;            % flip angle [rad]
    phase = 0;                      % phase [rad]
    
    T1 = 100;                       % longitudinal relaxation time [ms]
    T2 = 10;                        % transverse relaxation time [ms]
    D = 0;                          % we neglect diffusion in this example
    
    % initialize class
    cm = CoMo;      
    
    % set tissue parameters
    cm.R1 = 1 / T1;                 % longitudinal relaxation rate [1/ms]
    cm.R2 = 1 / T2;                 % transverse relaxation rate [1/ms]
    cm.D = D;                       % diffusion constant [um^2/ms]
    
    % discretization of configuration space
    % (ideally the greatest common divisor of all time intervals and/or gradient moments)
    cm.d_tau = TE;
    
    % set RF pulse parameters
    rf.FlipAngle = flip_angle;      % flip angle [rad]
    rf.Phase = phase;               % phase [rad]
    
    % set time interval from RF pulse to echo
    te.tau = TE;                    % duration
    
    % set time interval from echo to RF pulse (includes crusher)
    crusher.tau = TR - TE;          % duration
    
    % desired number of samples
    n_TR = 100;
    
    % allocate space for results
    m_xy = zeros( n_TR, 1 );
    
    % execute sequence loop
    for i = 1 : n_TR                
    
      cm.RF( rf );                  % RF pulse
    
      cm.time( te );                % time to echo
    
      % only magnetization pathways with no dephasing by crushers contribute to the FID signal.
      % this corresponds to cm.tau == 0, as explained in the article, linked above. 
      % to avoid possible rounding problems, we replace the condition cm.tau == 0 by a more robust one.
    
      fid = [];
      fid.b_n = cm.occ & ( abs( cm.tau ) < 0.1 * cm.d_tau );   % cm.occ denotes occupied configurations
    
      res = cm.sum( fid );          % extract FID
      m_xy( i ) = res.xy;           % store transverse component
    
      cm.time( crusher );           % time to next RF pulse (assumed to include a crusher gradient)
    
    end
    % done


# Please Note

CoMoTk constitutes a proof-of-concept implementation of the CCM.
Its numerical performance is limited, as it only supports single-core computation.
For serious projects consider porting it to a parallel architecture, since the time-consuming
part of the CCM iterations will benefit from it. 


# Releases

To supplement the submission of the background article about the CCM, a packaged, citable version of CoMoTk has been created:

![img](https://zenodo.org/badge/DOI/10.5281/zenodo.5347194.svg) (v0.43 - Aug. 2021)


# Feedback

Please use the [issue tracker](https://github.com/cganter/CoMoTk/issues) or write an email to [comotk.rad.med@tum.de](mailto:comotk.rad.med@tum.de).

Reporting bugs: Most helpful are small example scripts, which generate the error.

