%% Compares numerical implementation of magnetization exchange (ME) and transfer (MT)
% with exact analytical solution for the steady-state of balanced SSFP
% published in: Malik SJ et al., Magn Reson Med (2018) 80, 767-779

%% Simulation parameters
% the tissue parameters correspond to myelin water exchange (ME) and white matter (MT)
% as given in Table 1 of Malik et al. 

par = [];
opt = [];
str = [];

par.TR = 5;
par.fa = 10;
par.n_omega = 501;
par.prep = 5;

opt.TR = [];
opt.fa = [];
opt.n_omega = [];
opt.prep = [];

str.TR = '[ms] repetition time';
str.fa = '[deg] flip angle';
str.n_omega = 'number of off-resonance frequencies to be displayed';
str.prep = '[max(T1_a, T1_b)] duration of preparation phase'; 

while ( true )
    
    [ par, sel ] = sfv( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % variable input
    
    in = par;
    
    % constant parameters
    
    in.ME_T1a = 1000;
    in.ME_T2a = 100;
    in.ME_T1b = 500;
    in.ME_T2b = 20;
    in.ME_ka = 0.002;
    in.ME_f = 0.2;
    in.ME_delta_b = 12.8;
    
    in.ME_R1a = 1 / in.ME_T1a;
    in.ME_R2a = 1 / in.ME_T2a;
    in.ME_R1b = 1 / in.ME_T1b;
    in.ME_R2b = 1 / in.ME_T2b;
    in.ME_kb = ( 1 - in.ME_f ) * in.ME_ka / in.ME_f;
    
    in.MT_T1a = 779;
    in.MT_T2a = 45;
    in.MT_T1b = 779;
    in.MT_T2b = 0;
    in.MT_ka = 0.0043;
    in.MT_f = 0.117;
    in.MT_B1 = 13;
    in.MT_G0 = 15.1;
    
    in.MT_R1a = 1 / in.MT_T1a;
    in.MT_R2a = 1 / in.MT_T2a;
    in.MT_R1b = 1 / in.MT_T1b;
    in.MT_R2b = Inf;
    in.MT_kb = ( 1 - in.MT_f ) * in.MT_ka / in.MT_f;
   
    % now we calculate the RF saturation in MT
    % flip angle in radian
    
    in.fa_rad = par.fa * pi / 180;
    
    % gyromagnetic ratio [1 / (us * uT)]
    
    gamma = 42.58 * 1e-6;
    
    % rotation frequency corresponding to constant RF amplitude [rad/us]
    
    omega_rf = 2 * pi * gamma * in.MT_B1;
    
    % duration of RF pulse [us]
    
    tau_rf = in.fa_rad / omega_rf;
    
    % saturation of bound compartment in MT
    
    in.MT_sat = exp( - pi * omega_rf^2 * in.MT_G0 * tau_rf );

    %% set up configuration model
    
    ME_cm0 = CoMoTk;            % uncoupled spins
    ME_cm = CoMoTk;             % with ME or MT

    MT_cm0 = CoMoTk;            % uncoupled spins
    MT_cm = CoMoTk;             % with ME or MT
    
    % mandatory tissue parameters
    
    ME_cm0.R1 = [ in.ME_R1a, in.ME_R1b ];
    ME_cm0.R2 = [ in.ME_R2a, in.ME_R2b ];
    ME_cm0.D = [ 0, 0 ];
    ME_cm.R1 = [ in.ME_R1a, in.ME_R1b ];
    ME_cm.R2 = [ in.ME_R2a, in.ME_R2b ];
    ME_cm.D = [ 0, 0 ];
    
    MT_cm0.R1 = [ in.MT_R1a, in.MT_R1b ];
    MT_cm0.R2 = [ in.MT_R2a, in.MT_R2b ];
    MT_cm0.D = [ 0, 0 ];
    MT_cm.R1 = [ in.MT_R1a, in.MT_R1b ];
    MT_cm.R2 = [ in.MT_R2a, in.MT_R2b ];
    MT_cm.D = [ 0, 0 ];
    
    % proton density
    
    ME_cm0.mu = [ 1 - in.ME_f, in.ME_f ];
    ME_cm.mu = [ 1 - in.ME_f, in.ME_f ];
    
    MT_cm0.mu = [ 1 - in.MT_f, in.MT_f ];
    MT_cm.mu = [ 1 - in.MT_f, in.MT_f ];
    
    % off-resonance frequency [rad/ms]
    
    ME_cm0.dom = [ 0, 0.001 * 2 * pi * in.ME_delta_b ]; 
    ME_cm.dom = [ 0, 0.001 * 2 * pi * in.ME_delta_b ]; 
    
    % exchange rate matrix [1/ms]
    % k_ab corresponds to transition b -> a (== k_b, if defined as above)
    % only off-diagonal elements need to be specified, the diagonal element
    % are determined by the sum rule related to particle conservation
    
    ME_cm.k = [ 0, in.ME_kb; in.ME_ka, 0 ];
    
    MT_cm.k = [ 0, in.MT_kb; in.MT_ka, 0 ];
    
    % default CoMoTk options are ok, no need to change them
    
    % start with longitudinal magnetization
    
    ME_cm0.init_configuration ( [ 0, 0; 0, 0; 1 - in.ME_f, in.ME_f ] );
    ME_cm.init_configuration ( [ 0, 0; 0, 0; 1 - in.ME_f, in.ME_f ] );

    MT_cm0.init_configuration ( [ 0, 0; 0, 0; 1 - in.MT_f, in.MT_f ] );
    MT_cm.init_configuration ( [ 0, 0; 0, 0; 1 - in.MT_f, in.MT_f ] );

    %% number of TR cycles to approach steady state
    
    num_TR = ceil( par.prep * max( [ in.ME_T1a, in.ME_T1b, in.MT_T1a, in.MT_T1b ] ) / par.TR );
        
    %% prepare time between two RF pulses
    
    % unique index
    
    lam_t = 1;
    
    %% approach steady state
    
    for i = 0 : num_TR
        
        if ( i > 0 )
        
            % time interval
            
            param = [];
            param.lambda = lam_t;
            param.tau = par.TR;
            
            ME_cm0.time( param );
            ME_cm.time( param );
                        
            MT_cm0.time( param );
            MT_cm.time( param );
                        
        end
        
        % excitation pulse with alternating phase
        
        param = [];
        param.FlipAngle = in.fa_rad;
        param.Phase = mod( i, 2 ) * pi;
                
        ME_cm0.RF( param );
        ME_cm.RF( param );

        param.MT_sat = [ 0, in.MT_sat ];
        
        MT_cm0.RF( param );
        MT_cm.RF( param );

    end
    
    %% get and compare results
    
    % off-resonance frequencies
    
    if ( par.n_omega == 1 )
        
        in.omega = 0;
        
    else
        
        in.omega = linspace( - pi / par.TR, pi / par.TR, par.n_omega );
    
    end
        
    % allocate space for results
    
    ME_m0_ss = zeros( size( in.omega ) );
    ME_m_ss = zeros( size( in.omega ) );

    MT_m0_ss = zeros( size( in.omega ) );
    MT_m_ss = zeros( size( in.omega ) );

    % extract configurations for n = - n_max : n_max - 1
    
    param = [];
    
    for i = 1 : par.n_omega
        
        param.omega = in.omega( i );
        
        res = ME_cm0.sum( param );
        ME_m0_ss( i ) = res.xy;
        
        res = ME_cm.sum( param );
        ME_m_ss( i ) = res.xy;
        
        res = MT_cm0.sum( param );
        MT_m0_ss( i ) = res.xy;
        
        res = MT_cm.sum( param );
        MT_m_ss( i ) = res.xy;
        
    end
    
    %% calculate analytical result
    
    [ ME_m_ss_theory, MT_m_ss_theory ] = ME_MT_bSSFP_Malik( in );
    
    %% Show results
    
    theta = in.TR .* in.omega;
    
    ax = subplot( 1, 2, 1 );
    plot( theta, abs( MT_m0_ss ), 'g', theta, abs( MT_m_ss ), 'b', theta, abs( MT_m_ss_theory ), 'r--', 'LineWidth', 3 );

    xlim( [ -pi, pi ] );
    ylim( [ 0.03, 0.13 ] );
    xlabel( 'off-resonance [rad/TR]' );
    title( 'Magnetization Transfer' );
    legend( 'no coupling', 'CM', 'Theory', 'Interpreter', 'latex', 'Location', 'north' );
    
    width = 14;
    height = 7.5;
    
    set( gcf, 'Units', 'centimeters' );
    set( gcf, 'Position', [ 0, 0, width, height ] );
    set( gcf, 'Color', 'w' );
    
    delta = 0.01;
    
    ti = ax.TightInset;
    left = ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.5 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    ax = subplot( 1, 2, 2 );
    plot( theta, abs( ME_m0_ss ), 'g', theta, abs( ME_m_ss ), 'b', theta, abs( ME_m_ss_theory ), 'r--', 'LineWidth', 3 );
    
    xlim( [ -pi, pi ] );
    ylim( [ 0.03, 0.17 ] );
    xlabel( 'off-resonance [rad/TR]' );
    title( 'Magnetization Exchange' );
    legend( 'no coupling', 'CM', 'Theory', 'Interpreter', 'latex', 'Location', 'north' );
    
    ti = ax.TightInset;
    left = 0.5 + ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.5 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];

end

function [ ME_m_ss, MT_m_ss ] = ME_MT_bSSFP_Malik( in )
% ME_MT_bSSFP_Malik
%
% Calculates analytical solution for balanced steady-state SSFP
% according to Malik et al. (citation given above)
%
% In:
%
% in = structure with the following fields:
%
% sim = 'What to simulate: magnetization exchange (ME) or transfer (MT)';
% T1_a = '[ms] (MT: free protons)';
% T2_a = '[ms]';
% T1_b = '[ms] (MT: bound protons)';
% T2_b = '[ms]';
% k_a = '[1/ms] exchange rate a -> b';
% f = 'fraction of compartment b';
% delta_b = '[Hz] chemical shift of compartment b';
% TR = '[ms] repetition time';
% fa = '[deg] flip angle';
% omega = [rad/ms] vector of off-resonance frequencies
%
% Out:
%
% ME_m_ss = ME transverse steady-state magnetization, length == length( in.omega )
% MT_m_ss = same for MT

%% common parameters for ME and MT

c = cos( in.fa_rad );
s = sin( in.fa_rad );
c2 = cos( 0.5 * in.fa_rad )^2;
s2 = sin( 0.5 * in.fa_rad )^2;

% off-resonance frequencies

om = in.omega;
n_om = length( om );

% RF rotation matrix for single compartment

T11 = c2;
T12 = s2;
T13 = - 1i * s;
T21 = s2;
T22 = c2;
T23 = 1i * s;
T31 = - 0.5 * 1i * s;
T32 = 0.5 * 1i * s;
T33 = c;
    
%% specific settings

% [rad/ms] off-resonance frequency of compartment b

ME_id = 1i * 2 * pi * 0.001 * in.ME_delta_b; 

% transverse components contributing to signal

ME_Sig = [ 1, 0, 1, 0, 0, 0 ];

MT_Sig = [ 1, 0, 0, 0 ];

% two compartment RF matrix

ME_Th = [ T11, T12, 0, 0, T13, 0; ...
    T21, T22, 0, 0, T23, 0; ...
    0, 0, T11, T12, 0, T13; ...
    0, 0, T21, T22, 0, T23; ...
    T31, T32, 0, 0, T33, 0; ...
    0, 0, T31, T32, 0, T33 ];

MT_Th = [ ...
    T11, T12, T13, 0; ...
    T21, T22, T23, 0; ...
    T31, T32, T33, 0; ...
    0, 0, 0, in.MT_sat ];

% 180° rotation around z-axis

ME_D = [ -1, 0, 0, 0, 0, 0; ...
    0, -1, 0, 0, 0, 0; ...
    0, 0, -1, 0, 0, 0; ...
    0, 0, 0, -1, 0, 0; ...
    0, 0, 0, 0, 1, 0; ...
    0, 0, 0, 0, 0, 1 ];

MT_D = [ ...
    -1, 0, 0, 0; ...
    0, -1, 0, 0; ...
    0, 0, 1, 0; ...
    0, 0, 0, 1 ];

% repolarization term

ME_C = [ 0; 0; 0; 0; in.ME_R1a * ( 1 - in.ME_f ); in.ME_R1b * in.ME_f ];

MT_C = [ 0; 0; in.MT_R1a * ( 1 - in.MT_f ); in.MT_R1b * in.MT_f ];

% unit matrix

ME_I = eye( 6 );

MT_I = eye( 4 );
    
% allocate space for result

ME_m_ss = zeros( size( om ) );
MT_m_ss = zeros( size( om ) );

% calculate steady-state for all frequencies

for i = 1 : n_om
    
    iom = 1i * om( i );
    
    % relaxation, precession and exchange
    
    ME_A = [ ...
        - in.ME_R2a - in.ME_ka - iom, 0, in.ME_kb, 0, 0, 0; ...
        0, - in.ME_R2a - in.ME_ka + iom, 0, in.ME_kb, 0, 0; ...
        in.ME_ka, 0, - in.ME_R2b - in.ME_kb - ME_id - iom, 0, 0, 0; ...
        0, in.ME_ka, 0, - in.ME_R2b - in.ME_kb + ME_id + iom, 0, 0; ...
        0, 0, 0, 0, - in.ME_R1a - in.ME_ka, in.ME_kb;
        0, 0, 0, 0, in.ME_ka, - in.ME_R1b - in.ME_kb ];
    
    ME_exp_TR_A = expm( in.TR * ME_A );
    
    MT_A = [ ...
        - in.MT_R2a - iom, 0, 0, 0; ...
        0, - in.MT_R2a + iom, 0, 0; ...
        0, 0, - in.MT_R1a - in.MT_ka, in.MT_kb;
        0, 0, in.MT_ka, - in.MT_R1b - in.MT_kb ];
    
    MT_exp_TR_A = expm( in.TR * MT_A );
    
    % steady-state
    
    ME_m_ss( i ) = ME_Sig * ( ( ME_I - ME_Th * ME_D * ME_exp_TR_A ) \ ( ME_Th * ( ME_exp_TR_A - ME_I ) * ( ME_A \ ME_C ) ) );

    MT_m_ss( i ) = MT_Sig * ( ( MT_I - MT_Th * MT_D * MT_exp_TR_A ) \ ( MT_Th * ( MT_exp_TR_A - MT_I ) * ( MT_A \ MT_C ) ) );

end

end
