%% bSSFP sequence, ideal vs real RF pulse
% This script compares the effect of idealized (instantaneous) versus realistic RF pulses
% for a balanced SSFP sequence with TE = TR/2 and a pure tissue.


%% Simulation parameters

par = [];
opt = [];
str = [];

par.T1 = 500;
par.T2 = 50;
par.TR = 5;
par.num_TR = 50;
par.fa = 30;
par.sl_th = 1;
par.t_rf = 1;
par.supp_rf = 50;
par.qual_rf = 1 : 3;
par.filt_rf = 'Hamming';

opt.T1 = [];
opt.T2 = [];
opt.TR = [];
opt.num_TR = [];
opt.fa = [];
opt.sl_th = [];
opt.t_rf = [];
opt.supp_rf = [];
opt.qual_rf = [];
opt.filt_rf = { 'Hamming', 'None' };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.TR = '[ms] repetition time (echo at TR/2)';
str.num_TR = 'number of TR intervals to simulate';
str.fa = '[deg] flip angle';
str.sl_th = '[mm] slice thickness';
str.t_rf = '[ms] RF pulse duration';
str.supp_rf = 'number of RF support points';
str.qual_rf = 'number of SINC pulse zero crossings (on each side)';
str.filt_rf = 'filter RF pulse to reduce wiggles';

while ( true )
    
    [ par, select_conf ] = set_field_values( par, opt, str );
    
    if ( select_conf == -1 )
        
        break;
        
    end
    
    % TE = TR / 2
    
    TE = 0.5 * par.TR;
    
    % phase cycling
    
    ph_rad = zeros( par.num_TR, 1 );
    ph_rad( 2 : 2 : end ) = pi;
    
    %% apply settings
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    sl_th_um = 1000 * par.sl_th;
    
    % slice profile locations
    
    n_sl = 501;

    x = zeros( 3, n_sl );
    x( 3, : ) = linspace( - 1.5 * sl_th_um, 1.5 * sl_th_um, n_sl );
    
    % ideal slice profile at TE = TR / 2
    
    sp_ideal = zeros( 1, n_sl );
    sp_ideal( abs( x( 3, : ) ) < 0.5 * sl_th_um ) = sin( fa_rad );
    
    sp_ideal = sp_ideal .* exp( - 0.5 * par.TR / par.T2 );
                        
    %% ideal sequence with instantaneous RF pulses
    % initialize configuration model

    cm_ideal = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_ideal.R1 = 1 / par.T1;
    cm_ideal.R2 = 1 / par.T2;
    cm_ideal.D = 0;
    
    % allocated support in configuration space
    
    cm_ideal.d_tau = TE;
    cm_ideal.n_tau = 2 * num_TR;
    
    % prepare RF pulses
    
    RF_ideal = cell( 1, par.num_TR );
    
    for idx_TR = 1 : par.num_TR

        RF_ideal{ idx_TR }.FlipAngle = fa_rad;
        RF_ideal{ idx_TR }.Phase = ph_rad( idx_TR );
        
    end
    
    % prepare time periods

    RF_to_Echo_ideal = [];
    RF_to_Echo_ideal.tau = TE; % for balanced SSFP, we consider a centered echo
    
    Echo_to_RF_ideal = RF_to_Echo_ideal; % the two intervals are equivalent (duration & (zero) gradient moment)
                
    % allocate space for results
        
    m_ideal = zeros( par.num_TR, 1 );

    %% initialize timing
    
    t_id = zeros( par.num_TR, 1 );

    tic;
    t_tmp = toc;

    %%  loop over TR
    
    for idx_TR = 1 : par.num_TR
        
        fprintf( 1, 'i = %d / %d\n', idx_TR, par.num_TR );
                
        %% execute (b)SSFP sequence with ideal, instantaneous pulse ...
        
        if ( idx_TR > 1 ) % not needed before first RF pulse
            
            cm_ideal.time( Echo_to_RF_ideal );
            
        end
        
        cm_ideal.RF( RF_ideal{ idx_TR } );
        
        cm_ideal.time( RF_to_Echo_ideal );
          
        %% ... and extract the result
        
        % balanced SSFP has no crusher gradients
        % in absence of gradient dephasing, all occupied configurations contribute
        
        select_conf = [];
            
        % calculate the partial sum
        
        res = cm_ideal.sum( select_conf );
        
        % save the echo
        
        m_ideal( idx_TR ) = res.xy;
            
        % time spent for the ideal part
        
        t_id( idx_TR ) = - t_tmp;
        t_tmp = toc;
        t_id( idx_TR ) = t_id( idx_TR ) + t_tmp;        
        
        fprintf( 1, 'ideal RF = %9.3f sec\n', t_id( idx_TR ) );
        
    end
    
    fprintf( 1, 'ideal RF =\t %9.3f sec total\n', sum( t_id ) );
    
    %% now for the sequence(s) with realistic slice profiles
    % more than one profile with different number of zero-crossings are
    % allowed
    
    zc = sort( par.qual_rf );
    n_prof = length( zc );
    
    % allocate space for results
        
    m_real = zeros( par.num_TR, n_prof );
    
    % as described in the Handbook of MRI for a SINC pulse +/- Hamming filter
    
    % Number of RF support points
    %
    % as in SLR optimization we approximate the RF pulse as an alternating train of
    %
    % *   par.supp_rf + 1 instantaneous RF pulses, separated by
    % *   par.supp_rf     intervals of constant gradient amplitude
    
    % duration of small intervals
    
    tau_rf = par.t_rf / par.supp_rf;
    
    % allocate space for slice profile
    
    sp_real = zeros( n_prof, n_sl );

    % now things depend on the profile
    
    for ip = 1 : n_prof
         
        % support points and flip angles of instantaneous RF pulses (\propto B1 profile)
        
        t = zc( ip ) .* linspace( -1, 1, par.supp_rf + 1 );
        
        if ( isequal( par.filt_rf, 'Hamming' ) )
            
            al = ( 0.54 + 0.46 .* cos( pi .* t ./ zc( ip ) ) ) .* sinc( t );
            
        elseif ( isequal( par.filt_rf, 'None' ) )
            
            al = sinc( t );
            
        end
        
        % normalization: at the center of the slice profile, we need the desired flip angle
        
        al_rad = al .* ( fa_rad / sum( al ) );
        
        % approximate real RF pulses as sequence of small instantaneous pulses ...
        
        RF_real = cell( par.supp_rf + 1, par.num_TR );
        
        for idx_TR = 1 : par.num_TR
            
            for idx_RF = 1 : par.supp_rf + 1
                
                RF_real{ idx_RF, idx_TR }.FlipAngle = al_rad( idx_RF );
                RF_real{ idx_RF, idx_TR }.Phase  = ph_rad( idx_TR );
                
            end
            
        end
        
        % ... separated by small time intervals of constant duration and gradient moment
        
        DeltaTime = [];
        DeltaTime.tau = par.t_rf / par.supp_rf;
        
        % the associated gradient moment is defined by RF bandwidth, pulse shape and slice thickness:
        
        % RF pulse bandwidth
        
        f = 2 * zc( ip ) / par.t_rf;
        
        % total moment of slice selection gradient
        
        p_sl = 2 * pi * f * par.t_rf / sl_th_um;
        
        % split into par.supp_rf intervals
        
        DeltaTime.p = [ 0; 0; p_sl / par.supp_rf ];
        
        % prepare remaining time periods
        
        RF_to_Echo_real = [];
        
        RF_to_Echo_real.tau = ( par.TR - par.t_rf ) / 2; % for balanced SSFP, we consider a centered echo
        RF_to_Echo_real.p = [ 0; 0; - p_sl / 2 ];
        
        Echo_to_RF_real = RF_to_Echo_real; % the two intervals are equivalent (duration & (zero) gradient moment)
            
        %% initialize configuration model (real sequence)
        
        cm_real = CoMoTk;
        
        % mandatory tissue parameters
        
        cm_real.R1 = 1 / par.T1;
        cm_real.R2 = 1 / par.T2;
        cm_real.D = 0;
        
        % allocated support in configuration space
        
        cm_real.d_tau = DeltaTime.tau;
        cm_real.n_tau = 2 * num_TR * round( par.TR / DeltaTime.tau );
        cm_real.d_p = DeltaTime.p;
        cm_real.n_p = [ 0; 0; par.supp_rf * 2 * num_TR ];
        
        %% initialize timing
        
        t_re = zeros( par.num_TR, 1 );
        
        tic;
        t_tmp = toc;
        
        %%  loop over TR
        
        for idx_TR = 1 : par.num_TR
            
            fprintf( 1, 'Profile %d / %d, i = %d / %d\n', ip, n_prof, idx_TR, par.num_TR );
            
            %% execute (b)SSFP with real RF pulse ...
            
            if ( idx_TR > 1 ) % not needed before first RF pulse
                
                cm_real.time( Echo_to_RF_real );
                
            end
            
            for idx_RF = 1 : par.supp_rf
                
                cm_real.RF( RF_real{ idx_RF, idx_TR } );
                
                cm_real.time( DeltaTime );
                
            end
            
            cm_real.RF( RF_real{ par.supp_rf + 1, idx_TR } );
        
            cm_real.time( RF_to_Echo_real );

            %% ... and extract the result
            % in the first interval we extract the slice profile

            if ( idx_TR == 1 )
                                          
                for i = 1 : n_sl
                    
                    % we sum over all configurations (third argument == [])
                    
                    param = [];
                    param.omega = 0;
                    param.x = x( :, i );
                    
                    res = cm_real.sum( param );
                    
                    sp_real( ip, i ) = abs( res.xy );
                    
                end
                
            end
            
            % balanced SSFP has no crusher gradients
            % in absence of gradient dephasing, all occupied configurations contribute
            
            select_conf = [];
                
            % in both cases and unlike the real case, we *further* need to integrate over the slice profile
            
            select_conf.b_n = select_conf.b_n & ...
                reshape( abs( cm_real.p_n( 3, : ) ) < 0.5 * abs( DeltaTime.p( 3 ) ), size( select_conf.b_n ) );
            
            % calculate the partial sum
            
            res = cm_real.sum( select_conf );
            
            % save the echo
            
            m_real( idx_TR, ip ) = res.xy;
            
            % time spent for the real part
            
            t_re( idx_TR ) = - t_tmp;
            t_tmp = toc;
            t_re( idx_TR ) = t_re( idx_TR ) + t_tmp;
            
            fprintf( 1, 'real RF  = %9.3f sec\n', t_re( idx_TR ) );
            
        end
        
        fprintf( 1, 'real RF  =\t %9.3f sec total\n', sum( t_re ) );
        
    end
    
    %% Show results

    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        te = par.TE + ( 0 : par.num_TR - 1 )' .* par.TR;
        
    else
        
        te = par.TR / 2 + ( 0 : par.num_TR - 1 )' .* par.TR;
        
    end
   
    %%
    
    if ( isequal( par.ssfp, 'unbalanced' ) )

        m_ideal = abs( m_ideal );
        m_real = abs( m_real );
        
    else
        
        m_ideal = imag( m_ideal );
        m_real = imag( m_real );

    end
        
    loc = x( 3, : ) .* 0.001;
    
    %%
    
    hold off;
    
    ax = subplot( 1, 2, 1 );

    plot( te, m_ideal, 'DisplayName', 'ideal' );
    
    legend( '-DynamicLegend', 'Location', 'southeast' );
    
    hold all;
    
    for ip = n_prof : -1 : 1

        sc = abs( m_ideal( 1 ) / m_real( 1, ip ) );
        plot( te, sc .* m_real( :, ip ), 'DisplayName', [ 'zc = ', num2str( zc( ip ) ) ] );
        
        legend( 'Location', 'southeast' );

    end
     
    ylabel( '$m_{xy}''''$', 'Interpreter', 'latex' );
    xlabel( '$t$ [ms]', 'Interpreter', 'latex' );
    xlim( [ 0, te( end ) ] );
    ylim( [ -1.15, 0.7 ] );
    
    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        title( 'SSFP', 'Interpreter', 'latex' );
        
    elseif ( isequal( par.ssfp, 'balanced' ) )
        
        title( 'Transient bSSFP Oscillations', 'Interpreter', 'latex' );
        
    end
        
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
    
    hold off;

    ax = subplot( 1, 2, 2 );

    plot( loc, sp_ideal );
    
    % legend( '-DynamicLegend' );
    
    hold all;
    
    for ip = n_prof : -1 : 1

        % plot( loc, sp_real( ip, : ), 'DisplayName', [ 'zc = ', num2str( zc( ip ) ) ] );
        plot( loc, sp_real( ip, : ) );

    end
    
    title( 'Slice Profiles', 'Interpreter', 'latex' );
    xlabel( 'position [mm]', 'Interpreter', 'latex' );
    ylabel( '$\left|m_{xy}\right|$', 'Interpreter', 'latex' );    
    xlim( [ -1.5, 1.5 ] );
    
    ti = ax.TightInset;
    left = 0.5 + ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.5 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    %%
    
end
