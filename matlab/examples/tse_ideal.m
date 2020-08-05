%% single-shot TSE sequence
% This script compares idealized TSE with and without CPMG condition
% additionally, the CM is compared with EPG according to 
% Scheffler, Concepts Magn. Reson., 1999 (11) 291-304

%% Simulation parameters

opt = [];
par = [];
str = [];

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.resolution = [];
opt.n_echoes = [];
opt.dTE = [];
opt.fa_exc = [];
opt.fa_ref = [];
opt.B1 = [];
opt.show = { 'Magnitude', 'Phase', 'Real Part', 'Imaginary Part' };

par.T1 = 100;
par.T2 = 10;
par.D = 0;
par.resolution = 1;
par.n_echoes = 10;
par.dTE = 1;
par.fa_exc = 90;
par.fa_ref = 180;
par.B1 = 0.8;
par.show = opt.show{ 1 };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.resolution = '[mm] in-plane resolution (determines 2*pi crusher, relevant for D ~= 0 only)';
str.n_echoes = 'number of echoes per excitation';
str.dTE = '[ms] echo spacing';
str.fa_exc = '[deg] excitation flip angle (nominal, without B1+)';
str.fa_ref = '[deg] refocusing flip angle (nominal, without B1+)';
str.B1 = 'relative B1+';
str.show = 'What to display';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % convert to units, as expected by CoMoTk
    
    fa_exc_rad = par.fa_exc * pi / 180;
    fa_ref_rad = par.fa_ref * pi / 180;
    resolution_um = 1000 * par.resolution;
    
    %% initialize configuration model (idealized sequence)
    
    cm_cpmg = CoMoTk;
    cm_no_cpmg = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_cpmg.R1 = 1 / par.T1;
    cm_cpmg.R2 = 1 / par.T2;
    cm_cpmg.D = par.D;
    
    cm_no_cpmg.R1 = 1 / par.T1;
    cm_no_cpmg.R2 = 1 / par.T2;
    cm_no_cpmg.D = par.D;
    
    % further parameters
    
    cm_cpmg.B1 = par.B1;
    cm_no_cpmg.B1 = par.B1;
    
    % get default options
    
    options = cm_cpmg.options;
    
    options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
    options.alloc_d = 2;     % 1 == half echo spacing with crusher
    % 2 == time between last echo readout and next excitation pulse
    %      (includes cpmg spoiler)
    
    options.epsilon = 0;
    
    % set new options
    
    cm_cpmg.options = options;
    cm_no_cpmg.options = options;
    
    % start with longitudinal magnetization
    
    cm_cpmg.init_configuration ( [ 0; 0; 1 ] );
    cm_no_cpmg.init_configuration ( [ 0; 0; 1 ] );
    
    %% prepare time between end of RF pulse and echo
    
    % unique index
    
    lambda_te = 1;
    
    % duration
    
    tau_te = 0.5 * par.dTE;
    
    % gradient moment == crusher
    
    p_te = [ 2 * pi / resolution_um; 0; 0 ];
    
    %% allocate space for results
    
    m_cpmg = zeros( par.n_echoes, 1 );
    m_no_cpmg = zeros( par.n_echoes, 1 );
    
    %% initiate EPG for comparison
    
    m_EPG_cpmg = zeros( par.n_echoes, 1 );
    m_EPG_no_cpmg = zeros( par.n_echoes, 1 );
    
    n_EPG = 2 * par.n_echoes;
    EPG_zero = n_EPG + 1;
    
    % EPG states are stored as follows:
    % EPG( 1, n + EPG_zero ) = F( n )
    % EPG( 2, n + EPG_zero ) = conj( F( -n ) )
    % EPG( 3, n + EPG_zero ) = Z( n )
    
    EPG_cpmg = zeros( 3, 2 * n_EPG + 1 );
    EPG_no_cpmg = zeros( 3, 2 * n_EPG + 1 );
    
    % rotation matrix for excitation ...
    
    f_ = fa_exc_rad * par.B1;
    c = cos( f_ );
    s = sin( f_ );
    c22 = cos( 0.5 * f_ )^2;
    s22 = sin( 0.5 * f_ )^2;
    
    R_exc = [ ...
        c22, s22, - 1i * s; ...
        s22, c22, 1i * s; ...
        - 0.5 * 1i * s, 0.5 * 1i * s, c ];
    
    % and refocusing
    
    f_ = fa_ref_rad * par.B1;
    c = cos( f_ );
    s = sin( f_ );
    c22 = cos( 0.5 * f_ )^2;
    s22 = sin( 0.5 * f_ )^2;
    
    R_ref_no_cpmg = [ ...
        c22, s22, - 1i * s; ...
        s22, c22, 1i * s; ...
        - 0.5 * 1i * s, 0.5 * 1i * s, c ];
    
    z_ = 1i;
    
    R_ref_cpmg = [ ...
        c22, z_^2 * s22, - 1i * z_ * s; ...
        conj( z_^2 ) * s22, c22, 1i * conj( z_ ) * s; ...
        - 0.5 * 1i * conj( z_ ) * s, 0.5 * 1i * z_ * s, c ];
    
    % relaxation
    
    E1 = exp( - tau_te / par.T1 );
    E2 = exp( - tau_te / par.T2 );
    
    % initial state
    
    EPG_cpmg( 3, EPG_zero ) = 1;
    EPG_no_cpmg( 3, EPG_zero ) = 1;
    
    %% TSE with instantaneous RF pulse
    
    %% configuration model
    % excitation pulse
    
    param = [];
    param.FlipAngle = fa_exc_rad;
    param.Phase = 0;
    
    cm_cpmg.RF( param );
    cm_no_cpmg.RF( param );
    
    %% EPG
    % excitation pulse
    
    EPG_cpmg = R_exc * EPG_cpmg;
    EPG_no_cpmg = R_exc * EPG_no_cpmg;
    
    % counter
    n_occ = 0;
    
    for j = 1 : par.n_echoes
        
        %% configuration model
        % time to refocusing pulse
        
        param = [];
        param.lambda = lambda_te;
        param.tau = tau_te;
        param.p = p_te;
        
        cm_cpmg.time( param );
        cm_no_cpmg.time( param );
        
        % refocusing pulse
        
        param = [];
        param.FlipAngle = fa_ref_rad;
        param.Phase = 0.5 * pi;
        
        cm_cpmg.RF( param );
        
        param.Phase = 0;
        
        cm_no_cpmg.RF( param );
        
        % time to echo
        
        param = [];
        param.lambda = lambda_te;
        param.tau = tau_te;
        param.p = p_te;
        
        cm_cpmg.time( param );
        cm_no_cpmg.time( param );
        
        % select the correct configuration:
        % only the zero order configuration along the crusher direction contributes to the voxel signal
        
        param = [];
        param.b_n = cm_cpmg.find( lambda_te, 0 );
        
        % calculate the partial sum
        
        res = cm_cpmg.sum( param );
        
        % save the echo
        
        m_cpmg( j ) = res.xy;
        
        % no CPMG
        
        param = [];
        param.b_n = cm_no_cpmg.find( lambda_te, 0 );
        
        % calculate the partial sum
        
        res = cm_no_cpmg.sum( param );
        
        % save the echo
        
        m_no_cpmg( j ) = res.xy;
        
        %% EPG
        
        % time to refocusing pulse
        % relaxation
        
        rng = EPG_zero - n_occ : EPG_zero + n_occ;
        
        EPG_cpmg( 1, rng + 1 ) = E2 * EPG_cpmg( 1, rng );
        EPG_cpmg( 1, rng( 1 ) ) = 0;
        EPG_cpmg( 2, : ) = conj( EPG_cpmg( 1, end : -1 : 1 ) );
        EPG_cpmg( 3, : ) = E1 * EPG_cpmg( 3, : );

        EPG_no_cpmg( 1, rng + 1 ) = E2 * EPG_no_cpmg( 1, rng );
        EPG_no_cpmg( 1, rng( 1 ) ) = 0;
        EPG_no_cpmg( 2, : ) = conj( EPG_no_cpmg( 1, end : -1 : 1 ) );
        EPG_no_cpmg( 3, : ) = E1 * EPG_no_cpmg( 3, : );
        
        % repolarization
        
        EPG_cpmg( 3, EPG_zero ) = EPG_cpmg( 3, EPG_zero ) + 1 - E1;
        EPG_no_cpmg( 3, EPG_zero ) = EPG_no_cpmg( 3, EPG_zero ) + 1 - E1;
        
        % increase the counter
        
        n_occ = n_occ + 1;
        
        % refocusing pulse
        
        EPG_cpmg = R_ref_cpmg * EPG_cpmg;
        EPG_no_cpmg = R_ref_no_cpmg * EPG_no_cpmg;
        
        % time to echo
        % relaxation
        
        rng = EPG_zero - n_occ : EPG_zero + n_occ;
        
        EPG_cpmg( 1, rng + 1 ) = E2 * EPG_cpmg( 1, rng );
        EPG_cpmg( 1, rng( 1 ) ) = 0;
        EPG_cpmg( 2, : ) = conj( EPG_cpmg( 1, end : -1 : 1 ) );
        EPG_cpmg( 3, : ) = E1 * EPG_cpmg( 3, : );

        EPG_no_cpmg( 1, rng + 1 ) = E2 * EPG_no_cpmg( 1, rng );
        EPG_no_cpmg( 1, rng( 1 ) ) = 0;
        EPG_no_cpmg( 2, : ) = conj( EPG_no_cpmg( 1, end : -1 : 1 ) );
        EPG_no_cpmg( 3, : ) = E1 * EPG_no_cpmg( 3, : );
 
        % repolarization
        
        EPG_cpmg( 3, EPG_zero ) = EPG_cpmg( 3, EPG_zero ) + 1 - E1;
        EPG_no_cpmg( 3, EPG_zero ) = EPG_no_cpmg( 3, EPG_zero ) + 1 - E1;
        
        % increase the counter
        
        n_occ = n_occ + 1;
        
        % select the correct EPG state F( 0 ):
        
        % save the echo
        
        m_EPG_cpmg( j ) = EPG_cpmg( 1, EPG_zero );
        m_EPG_no_cpmg( j ) = EPG_no_cpmg( 1, EPG_zero );
        
    end
    
    %% Compare the decays of the last cycle
    
    te = ( 1 : par.n_echoes )' .* par.dTE;
    
    if ( isequal( par.show, 'Magnitude' ) )
        
        cpmg = abs( m_cpmg );
        no_cpmg = abs( m_no_cpmg );
        cpmg_EPG = abs( m_EPG_cpmg );
        no_cpmg_EPG = abs( m_EPG_no_cpmg );
        
    elseif ( isequal( par.show, 'Phase' ) )
        
        cpmg = angle( m_cpmg );
        no_cpmg = angle( m_no_cpmg );
        cpmg_EPG = angle( m_EPG_cpmg );
        no_cpmg_EPG = angle( m_EPG_no_cpmg );
        
    elseif ( isequal( par.show, 'Real Part' ) )
        
        cpmg = real( m_cpmg );
        no_cpmg = real( m_no_cpmg );
        cpmg_EPG = real( m_EPG_cpmg );
        no_cpmg_EPG = real( m_EPG_no_cpmg );
        
    elseif ( isequal( par.show, 'Imaginary Part' ) )
        
        cpmg = imag( m_cpmg );
        no_cpmg = imag( m_no_cpmg );
        cpmg_EPG = imag( m_EPG_cpmg );
        no_cpmg_EPG = imag( m_EPG_no_cpmg );
        
    end
    
    subplot( 1, 1, 1 );
    plot( te, cpmg, 'b', te, cpmg_EPG, 'b+', te, no_cpmg, 'r', te, no_cpmg_EPG, 'r+' );
    legend( 'CM (CPMG)', 'EPG (CPMG)', 'CM (no CPMG)', 'EPG (no CPMG)' );
    title( 'ideal TSE' );
    
end
