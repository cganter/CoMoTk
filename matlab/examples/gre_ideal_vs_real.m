%% GRE sequence
% This script compares an idealized versus realistic implementation of a GRE sequence
% no diffusion or spin coupling

%% Simulation parameters

par = [];
opt = [];
str = [];

par.ssfp = 'balanced';
par.T1 = 1000;
par.T2 = 100;
par.D = 0;
par.TR = 5;
par.TE = 2;
par.num_TR = 5;
par.B1 = 1;
par.fa = 50;
par.ph_inc = 180;
par.phd_inc = 0;
par.sl_th = 1;
par.t_rf = 1;
par.supp_rf = 100;
par.qual_rf = 3;
par.filt_rf = 'Hamming';
par.epsilon = 1e-4;
par.verbose = 'False';

opt.ssfp = { 'balanced', 'unbalanced' };
opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.TR = [];
opt.TE = [];
opt.num_TR = [];
opt.B1 = [];
opt.fa = [];
opt.ph_inc = [];
opt.phd_inc = [];
opt.sl_th= [];
opt.t_rf = [];
opt.supp_rf = [];
opt.qual_rf = [];
opt.filt_rf = { 'Hamming', 'None' };
opt.epsilon = [];
opt.verbose = { 'True', 'False' };

str.ssfp = 'SSFP variant to simulate';
str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.TR = '[ms] repetition time';
str.TE = '[ms] echo time (t_rf / 2 <= TE <= TR - t_rf / 2)';
str.num_TR = 'number of TR intervals to simulate';
str.B1 = 'relative B1+';
str.fa = '[deg] flip angle';
str.ph_inc = '[deg] RF phase increment (same as negative off-resonance)' ;
str.phd_inc = '[deg] RF phase difference increment (RF spoiling)';
str.sl_th= '[mm] slice thickness';
str.t_rf = '[ms] RF pulse duration';
str.supp_rf = 'number of RF support points';
str.qual_rf = 'number of SINC pulse zero crossings (on each side)';
str.filt_rf = 'filter RF pulse to reduce wiggles';
str.epsilon = 'discard configuration vectors with L2 norm smaller than this';
str.verbose = 'provide some informal output';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    %% apply settings
    % calculate phase cycling
    
    ph_deg = zeros( par.num_TR, 1 );
    
    for i = 2 : par.num_TR
        
        ph_deg( i ) = ph_deg( i - 1 ) + par.ph_inc + ( i - 1 ) * par.phd_inc;
        
    end
    
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    ph_rad = ph_deg .* ( pi / 180 );
    sl_th_um = 1000 * par.sl_th;
    
    %% set the RF pulse profile
    % as described in the Handbook of MRI for a SINC pulse +/- Hamming filter
    
    % Number of RF support points
    %
    % as in SLR optimization we approximate the RF pulse as an alternating train of
    %
    % *   par.supp_rf + 1 instantaneous RF pulses, separated by
    % *   par.supp_rf     intervals of constant gradient amplitude
    
    % duration of small intervals
    
    tau_rf = par.t_rf / par.supp_rf;
    
    % support points and flip angles of instantaneous RF pulses (\propto B1 profile)
    
    t = par.qual_rf .* linspace( -1, 1, par.supp_rf + 1 );
    
    if ( isequal( par.filt_rf, 'Hamming' ) )
        
        al = ( 0.54 + 0.46 .* cos( pi .* t ./ par.qual_rf ) ) .* sinc( t );
        
    else
        
        al = sinc( t );
        
    end
    
    % normalization: at the center of the slice profile, we need the desired flip angle
    
    al_rad = al .* ( fa_rad / sum( al ) );
    
    % set verbosity
    
    if ( isequal( par.verbose, 'True' ) )
        
        verbose = true;
        
    else
        
        verbose = false;
        
    end
    
    
    %% initialize configuration model (idealized sequence)

    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        % ======== unbalanced ========
        
        cm_SSFP_ideal = CoMoTk;
        
        % mandatory tissue parameters
        
        cm_SSFP_ideal.R1 = 1 / par.T1;
        cm_SSFP_ideal.R2 = 1 / par.T2;
        cm_SSFP_ideal.D = par.D;
        
        % further parameters
        
        cm_SSFP_ideal.B1 = par.B1;
        
        % set options
        
        options = cm_SSFP_ideal.options;
        options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
        options.alloc_d = 2;        % before and after echo
        options.epsilon = par.epsilon;
        cm_SSFP_ideal.options = options;
        
        % start with longitudinal magnetization
        
        cm_SSFP_ideal.init_configuration ( [ 0; 0; 1 ] );
        
        % prepare time between end of RF pulse and echo
        % unique index
        
        mu_SSFP_ideal_pre = 1;
        mu_SSFP_ideal_post = 2;
        
        % duration
        
        tau_SSFP_ideal_pre = par.TE;
        tau_SSFP_ideal_post = par.TR - par.TE;
        
    else
        
        % ======== balanced ========
        
        cm_bSSFP_ideal = CoMoTk;
        
        % mandatory tissue parameters
        
        cm_bSSFP_ideal.R1 = 1 / par.T1;
        cm_bSSFP_ideal.R2 = 1 / par.T2;
        cm_bSSFP_ideal.D = par.D;
        
        % further parameters
        
        cm_bSSFP_ideal.B1 = par.B1;
        
        % set options
        
        options = cm_bSSFP_ideal.options;
        options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
        options.alloc_d = 1;        % TE = TR / 2
        options.epsilon = par.epsilon;
        cm_bSSFP_ideal.options = options;
        
        % start with longitudinal magnetization
        
        cm_bSSFP_ideal.init_configuration ( [ 0; 0; 1 ] );
        
        % prepare time between end of RF pulse and echo
        % unique index
        
        mu_bSSFP_ideal = 1;
        
        % duration
        
        tau_bSSFP_ideal = par.TR / 2;
        
    end
    
    %% initialize configuration model (real sequence)
    
    % ======== prepare RF pulse plateau ========
    
    % unique index
    
    mu_real_rf = 1;
    
    % duration
    
    tau_real_rf = par.t_rf / par.supp_rf;
    
    % RF pulse bandwidth
    
    f = 2 * par.qual_rf / par.t_rf;
    
    % total moment of slice selection gradient
    
    p_sl = 2 * pi * f * par.t_rf / sl_th_um;
    
    % split into par.supp_rf intervals
    
    p_real_rf = [ 0; 0; p_sl / par.supp_rf ];
    
    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        % ======== unbalanced ========
        
        cm_SSFP_real = CoMoTk;
        
        % mandatory tissue parameters
        
        cm_SSFP_real.R1 = 1 / par.T1;
        cm_SSFP_real.R2 = 1 / par.T2;
        cm_SSFP_real.D = par.D;
        
        % further parameters
        
        cm_SSFP_real.B1 = par.B1;
        
        % set options
        
        options = cm_SSFP_real.options;
        options.alloc_n = 10000;  % CoMoTk will allocate more, if needed
        options.alloc_d = 3;          % 1: RF pulse plateau
        % 2,3: before and after echo
        options.epsilon = par.epsilon;
        options.verbose = verbose;
        cm_SSFP_real.options = options;
        
        % start with longitudinal magnetization
        
        cm_SSFP_real.init_configuration ( [ 0; 0; 1 ] );
        
        % prepare times between RF pulse and echo
        % unique index
        
        mu_SSFP_real_pre = 2;
        mu_SSFP_real_post = 3;
        
        % duration
        
        tau_SSFP_real_pre = par.TE - par.t_rf / 2;
        tau_SSFP_real_post = par.TR - par.TE - par.t_rf / 2;
        
        % gradient moment == crusher and rephasing gradient
        
        p_SSFP_real_pre = [ 0; 0; - p_sl / 2 ];
        p_SSFP_real_post = [ 0; 0; - p_sl / 2 ];
        
    else
        
        % ======== balanced ========
        
        cm_bSSFP_real = CoMoTk;
        
        % mandatory tissue parameters
        
        cm_bSSFP_real.R1 = 1 / par.T1;
        cm_bSSFP_real.R2 = 1 / par.T2;
        cm_bSSFP_real.D = par.D;
        
        % further parameters
        
        cm_bSSFP_real.B1 = par.B1;
        
        % set options
        
        options = cm_bSSFP_real.options;
        options.alloc_n = 10000;  % CoMoTk will allocate more, if needed
        options.alloc_d = 2;          % 1: RF pulse plateau
        % 2: TE = TR / 2
        options.epsilon = par.epsilon;
        cm_bSSFP_real.options = options;
        
        % start with longitudinal magnetization
        
        cm_bSSFP_real.init_configuration ( [ 0; 0; 1 ] );
        
        % prepare times between RF pulse and echo
        % unique index
        
        mu_bSSFP_real = 2;
        
        % duration
        
        tau_bSSFP_real = ( par.TR - par.t_rf ) / 2;
        
        % gradient moment == crusher and rephasing gradient
        
        p_bSSFP_real = [ 0; 0; - p_sl / 2 ];
        
        %% allocate space for results
        
        m_SSFP_ideal = zeros( par.num_TR, 1 );
        m_bSSFP_ideal = zeros( par.num_TR, 1 );
        m_SSFP_real = zeros( par.num_TR, 1 );
        m_bSSFP_real = zeros( par.num_TR, 1 );
        
    end
            
    t_id = zeros( par.num_TR, 1 );
    t_re = zeros( par.num_TR, 1 );

    tic;
    t_tmp = toc;

    for i = 1 : par.num_TR
        
        fprintf( 1, 'i = %d / %d\n', i, par.num_TR );
                
        %% (b)SSFP with instantaneous RF pulse
        
        if ( isequal( par.ssfp, 'unbalanced' ) )
                        
            % ======== unbalanced ========
            
            % time after echo
            % placed here, to initialize the crusher time interval
            % otherwise b_n cannot be set for i == 1 below
            % since the initial state is in equilibrium for i == 1, the magnetization is not changed
            
            if ( i > 1 )
                
                param = [];
                param.mu = mu_SSFP_ideal_post;
                param.tau = tau_SSFP_ideal_post;
                
                cm_SSFP_ideal.time( param );
                
            end
            
            % excitation pulse
            
            param = [];
            param.FlipAngle = fa_rad;
            param.Phase = ph_rad( i );
            
            cm_SSFP_ideal.RF( param );
            
            % time to echo
            
            param = [];
            param.mu = mu_SSFP_ideal_pre;
            param.tau = tau_SSFP_ideal_pre;
            
            cm_SSFP_ideal.time( param );
            
            % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
            %         (assuming a crusher was present after the echo - even, if not simulated)
            
            param = [];
            param.b_n = cm_SSFP_ideal.find( mu_SSFP_ideal_post, 0 );
            
            % calculate the partial sum
            
            res = cm_SSFP_ideal.sum( param );
            
            % save the echo
            
            m_SSFP_ideal( i ) = res.xy;
            
        else
            
            % ======== balanced ========
            
            % time after echo
            
            if ( i > 1 )
                
                param = [];
                param.mu = mu_bSSFP_ideal;
                param.tau = tau_bSSFP_ideal;
                
                cm_bSSFP_ideal.time( param );
                
            end
            
            % excitation pulse
            
            param = [];
            param.FlipAngle = fa_rad;
            param.Phase = ph_rad( i );
            
            cm_bSSFP_ideal.RF( param );
            
            % time to echo
            
            param = [];
            param.mu = mu_bSSFP_ideal;
            param.tau = tau_bSSFP_ideal;
            
            cm_bSSFP_ideal.time( param );
            
            % bSSFP : all configurations contribute to the voxel signal
            
            param = [];
            
            res = cm_bSSFP_ideal.sum( param );
            
            % save the echo
            
            m_bSSFP_ideal( i ) = res.xy;
                        
        end

        t_id( i ) = - t_tmp;
        t_tmp = toc;
        t_id( i ) = t_id( i ) + t_tmp;        
        
        %% (b)SSFP with real RF pulse
        
        if ( isequal( par.ssfp, 'unbalanced' ) )
            
            % ======== unbalanced ========
            
            % time after echo
            
            if ( i > 1 )
                
                param = [];
                param.mu = mu_SSFP_real_post;
                param.tau = tau_SSFP_real_post;
                param.p = p_SSFP_real_post;
                
                cm_SSFP_real.time( param );
                
            end
            
            % excitation pulse
            
            % first small pulse
            
            param = [];
            param.FlipAngle = al_rad( 1 );
            param.Phase = ph_rad( i );
            
            cm_SSFP_real.RF( param );
            
            for j = 1 : par.supp_rf
                
                param = [];
                param.mu = mu_real_rf;
                param.tau = tau_real_rf;
                param.p = p_real_rf;
                
                cm_SSFP_real.time( param );
                
                % and the rest of the small pulses
                
                param = [];
                param.FlipAngle = al_rad( j + 1 );
                param.Phase = ph_rad( i );
                
                cm_SSFP_real.RF( param );
                
            end
            
            % time to echo
            
            param = [];
            param.mu = mu_SSFP_real_pre;
            param.tau = tau_SSFP_real_pre;
            param.p = p_SSFP_real_pre;
            
            cm_SSFP_real.time( param );
            
            % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
            %         (assuming a crusher was present after the echo - even, if not simulated)
            
            b_n = cm_SSFP_real.find( mu_SSFP_real_post, 0 );
            
            % slice encoding direction requires a zero gradient moment
            
            param = [];
            param.b_n = b_n & reshape( abs( cm_SSFP_real.p_n( 3, : ) ) < 0.5 * p_real_rf( 3 ), size( b_n ) );
            
            % calculate the partial sum
            
            res = cm_SSFP_real.sum( param );
            
            % save the echo
            
            m_SSFP_real( i ) = res.xy;
            
        else
            
            % ======== balanced ========
            
            % time after echo
            
            if ( i > 1 )
                
                param = [];
                param.mu = mu_bSSFP_real;
                param.tau = tau_bSSFP_real;
                param.p = p_bSSFP_real;
                
                cm_bSSFP_real.time( param );
                
            end
            
            % excitation pulse
            
            % first small pulse
            
            param = [];
            param.FlipAngle = al_rad( 1 );
            param.Phase = ph_rad( i );
            
            cm_bSSFP_real.RF( param );
            
            for j = 1 : par.supp_rf
                
                param = [];
                param.mu = mu_real_rf;
                param.tau = tau_real_rf;
                param.p = p_real_rf;
                
                cm_bSSFP_real.time( param );
                
                % and the rest of the small pulses
                
                param = [];
                param.FlipAngle = al_rad( j + 1 );
                param.Phase = ph_rad( i );
                
                cm_bSSFP_real.RF( param );
                
            end
            
            % time to echo
            
            param = [];
            param.mu = mu_bSSFP_real;
            param.tau = tau_bSSFP_real;
            param.p = p_bSSFP_real;
            
            cm_bSSFP_real.time( param );
            
            % bSSFP : only the slice encoding direction requires a zero gradient moment
            
            param = [];
            param.b_n = cm_bSSFP_real.b_n & reshape( abs( cm_bSSFP_real.p_n( 3, : ) ) < 0.5 * p_real_rf( 3 ), size( cm_bSSFP_real.b_n ) );
            
            % calculate the partial sum
            
            res = cm_bSSFP_real.sum( param );
            
            % save the echo
            
            m_bSSFP_real( i ) = res.xy;
            
        end
        
        t_re( i ) = - t_tmp;
        t_tmp = toc;
        t_re( i ) = t_re( i ) + t_tmp;        
        
        fprintf( 1, 'ideal RF = %9.3f sec\n', t_id( i ) );
        fprintf( 1, 'real RF  = %9.3f sec\n', t_re( i ) );
        
    end
    
    fprintf( 1, 'ideal RF =\t %9.3f sec total\n', sum( t_id ) );
    fprintf( 1, 'real RF  =\t %9.3f sec total\n', sum( t_re ) );
    
    %% Show results
    
    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        te_SSFP = par.TE + ( 0 : par.num_TR - 1 )' .* par.TR;
        sc_SSFP = abs( m_SSFP_ideal( 1 ) / m_SSFP_real( 1 ) );
        
        subplot( 1, 1, 1 );
        plot( te_SSFP, abs( m_SSFP_ideal ), te_SSFP, sc_SSFP .* abs( m_SSFP_real ) );
        legend( 'ideal', 'real' );
        title( 'SSFP' );
        
    else
        
        te_bSSFP = par.TR / 2 + ( 0 : par.num_TR - 1 )' .* par.TR;
        sc_bSSFP = abs( m_bSSFP_ideal( 1 ) / m_bSSFP_real( 1 ) );
        
        subplot( 1, 1, 1 );
        plot( te_bSSFP, abs( m_bSSFP_ideal ), te_bSSFP, sc_bSSFP .* abs( m_bSSFP_real ) );
        legend( 'ideal', 'real' );
        title( 'bSSFP' );
        
    end
    
end
