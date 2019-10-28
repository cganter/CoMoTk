%% TSE sequence
% This script compares idealized TSE with and without CPMG condition

%% Simulation parameters

opt = [];
par = [];
str = [];

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.resolution = [];
opt.n_echoes = [];
opt.n_TR = [];
opt.TR = [];
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
par.n_TR = 1;
par.TR = 100;
par.dTE = 1;
par.fa_exc = 90;
par.fa_ref = 180;
par.B1 = 1;
par.show = opt.show{ 1 };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.resolution = '[mm] in-plane resolution (determines 2*pi crusher, relevant for D ~= 0 only)';
str.n_echoes = 'number of echoes per excitation';
str.n_TR = 'number of TR cycles (e.g. to establish steady state)';
str.TR = '[ms] repetition time (relevant, if n_TR > 1 only)';
str.dTE = '[ms] echo spacing';
str.fa_exc = '[deg] excitation flip angle (nominal, without B1+)';
str.fa_ref = '[deg] refocusing flip angle (nominal, without B1+)';
str.B1 = 'relative B1+';
str.show = 'What to display';

while ( true )
    
    [ par, sel ] = sfv( par, opt, str );
    
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
    
    mu_te = 1;
    
    % duration
    
    tau_te = 0.5 * par.dTE;
    
    % gradient moment == crusher
    
    p_te = [ 2 * pi / resolution_um; 0; 0 ];
    
    %% prepare time between last echo readout and next excitation pulse
    
    % unique index
    
    mu_tr = 2;
    
    % duration
    
    tau_tr = par.TR - par.n_echoes * par.dTE;
    
    % gradient moment == Inf (cpmg spoiler, even if D == 0)
    
    p_tr = [ Inf; 0; 0 ];
    
    %% allocate space for results
    
    m_cpmg = zeros( par.n_echoes, par.n_TR );
    m_no_cpmg = zeros( par.n_echoes, par.n_TR );
    
    %% TSE with instantaneous RF pulse
    
    for i = 1 : par.n_TR
        
        % excitation pulse
        
        param = [];
        param.FlipAngle = fa_exc_rad;
        param.Phase = 0;
        
        cm_cpmg.RF( param );
        cm_no_cpmg.RF( param );
        
        for j = 1 : par.n_echoes
            
            % time to refucusing pulse
            
            param = [];
            param.mu = mu_te;
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
            param.mu = mu_te;
            param.tau = tau_te;
            param.p = p_te;
            
            cm_cpmg.time( param );
            cm_no_cpmg.time( param );
            
            % select the correct configuration:
            % only the zero order configuration along the crusher direction contributes to the voxel signal
            
            param = [];
            param.b_n = cm_cpmg.find( mu_te, 0 );
            
            % calculate the partial sum
            
            res = cm_cpmg.sum( param );
            
            % save the echo
            
            m_cpmg( j, i ) = res.xy;
            
            % no CPMG
            
            param = [];
            param.b_n = cm_no_cpmg.find( mu_te, 0 );
            
            % calculate the partial sum
            
            res = cm_no_cpmg.sum( param );
            
            % save the echo
            
            m_no_cpmg( j, i ) = res.xy;
            
        end
        
        if ( i < par.n_TR )
            
            % time to next excitation pulse
            
            param = [];
            param.mu = mu_tr;
            param.tau = tau_tr;
            param.p = p_tr;
            
            cm_cpmg.time( param);
            cm_no_cpmg.time( param );
            
        end
        
    end
    
    %% Compare the decays of the last cycle
    
    te = ( 1 : par.n_echoes )' .* par.dTE;
    
    if ( isequal( par.show, 'Magnitude' ) )
        
        cpmg = abs( m_cpmg( :, end ) );
        no_cpmg = abs( m_no_cpmg( :, end ) );
    
    elseif ( isequal( par.show, 'Phase' ) )
        
        cpmg = angle( m_cpmg( :, end ) );
        no_cpmg = angle( m_no_cpmg( :, end ) );
    
    elseif ( isequal( par.show, 'Real Part' ) )
        
        cpmg = real( m_cpmg( :, end ) );
        no_cpmg = real( m_no_cpmg( :, end ) );
    
    elseif ( isequal( par.show, 'Imaginary Part' ) )
        
        cpmg = imag( m_cpmg( :, end ) );
        no_cpmg = imag( m_no_cpmg( :, end ) );
    
    end
        
    subplot( 1, 2, 1 );
    plot( te, cpmg, te, no_cpmg );
    legend( 'CPMG', 'no CPMG' );
    title( 'ideal TSE' );
    
    subplot( 1, 2, 2 );
    semilogy( te, cpmg, te, no_cpmg );
    legend( 'CPMG', 'no CPMG' );
    title( 'ideal TSE' );
    
end
