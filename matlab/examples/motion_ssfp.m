%% Equivalence between RF spoiling and constant motion along crusher gradient
% We calculate the steady state FID of an unbalanced SSFP sequence
% (instantaneous RF pulses with constant phase, i.e. no RF spoiling) for
% spins moving along the crusher gradient with some constant velocity. We
% demonstrate that the signal amplitude is equivalent to that of a static
% spin in presence of RF spoiling.

%% Simulation parameters
% (by default, all set by the user)

par = [];
opt = [];
str = [];

par.T1 = [ 100, 500 ];
par.T2 = [ 50, 250 ];
par.D = 0;
par.TR = 5;
par.fa = 90;
par.resolution = 1000;
par.v_max = 120;
par.n_v = 20;
par.t_prep = 10;

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.TR = [];
opt.fa = [];
opt.resolution = [];
opt.v_max = [];
opt.n_v = [];
opt.t_prep = [];

str.T1 = '[ms] set of T1 values';
str.T2 = '[ms] set of T2 values';
str.D = '[um^2/ms] isotropic ADC';
str.TR = '[ms] repetition time';
str.fa = '[deg] flip angle';
str.resolution = '[um] voxel size (relevant for the moment of the 2 * pi crusher)';
str.v_max = '[mm/s] maximum (absolute) velocity component in direction of crusher';
str.n_v = 'number of velocities in the range [ 0, v_max ]';
str.t_prep = 'duration of preparation phase in units of T1';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % check that T1 and T2 have the same length
    
    n_T1 = length( par.T1 );
    n_T2 = length( par.T2 );

    if ( n_T1 ~= n_T2 )
        
        fprintf( 1, 'T1 and T2 must have the same length.\n' );
        continue;
        
    end
        
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    
    % to test the implementation, we individually set the directions of the crusher
    % gradient and the velocity at random
    
    ec = randn( 3, 1 );
    ec = ec ./ norm( ec );
    
    ev = randn( 3, 1 );
    ev = ev ./ norm( ev );
    
    % moment of (constant) crusher gradient, effecting 2 * pi dephasing per TR and resolution
    
    pc = ec ./ par.resolution;
    
    % we need the scalar product of 'ec' and 'ev' for the velocity component along the
    % crusher gradient (hoping that it is not exactly zero by chance)
    
    vs = linspace( 0, par.v_max, par.n_v );
    v = ev .* vs ./ ( ec' * ev );
        
    % allocate space for results
    
    xy_CM = zeros( par.n_v, n_T2 );
    xy_TH = zeros( par.n_v, n_T2 );
    
    for i_Tx = 1 : n_T2
        
        %% prepare configuration model
        % To speed up the computation, we assume ideal RF pulses and a zero
        % echo time. The CM has then only one dimension (the time interval,
        % which contains the constant crusher gradient)
                
        % RF parameters
        
        RF_par = [];
        RF_par.FlipAngle = fa_rad;
        RF_par.Phase = 0;
        
        % parameters for time interval
 
        lambda_crusher = 1;     % unique index

        Time_par = [];
        Time_par.lambda = lambda_crusher;
        Time_par.tau = par.TR;
        Time_par.p = pc;
        
        % number of TR cycles to approach steady state
        
        n_TR = ceil( par.t_prep * par.T1( i_Tx ) / par.TR );
       
        for i_v = 1 : par.n_v
        
            fprintf( 1, 'T2 = %d / %d, v = %d / %d\n', i_Tx, n_T2, i_v, par.n_v );
                        
            % set actual velocity and position at beginning of time
            % interval
            
            Time_par.v = v( :, i_v );
            Time_par.x = zeros( 3, 1 );
            
            %% initialize configuration model (idealized sequence)
            
            cm = CoMoTk;        % create instance
            
            % mandatory tissue parameters
            
            cm.R1 = 1 / par.T1( i_Tx );
            cm.R2 = 1 / par.T2( i_Tx );
            cm.D = par.D;
            
            %% approach steady state
            
            for i_TR = 1 : n_TR
                
                % excitation pulse
                
                cm.RF( RF_par );
                
                if ( i_TR == n_TR )
                    
                    % in the last interval we extract the steady-state
                    % we assume a conventional FID sequence with crusher AFTER the echo
                    
                    % the associated configuration must be zero
                    % since the signal from nonzero orders is (hopefully) dephased
                    
                    sel_conf = [];
                    sel_conf.b_n = cm.find( lambda_crusher, 0 );
                    
                    % calculate the partial sum
                    
                    res = cm.sum( sel_conf );
                    xy_CM( i_v, i_Tx ) = abs( res.xy );
                    
                else
                    
                    % otherwise, we execute the time interval
                    
                    cm.time( Time_par );
                    
                    % we move ahead..
                    
                    Time_par.x = Time_par.x + par.TR .* Time_par.v;
                    
                end
                
            end
            
            % now we calculate the corresponding RF spoiled steady-state
            % for static spins
            
            in = [];
            in.FA = par.fa;
            in.TR = par.TR;
            in.R1 = 1 / par.T1( i_Tx );
            in.R2 = 1 / par.T2( i_Tx );
            in.PDI = ( 180 / pi ) * par.TR * ( Time_par.p' * Time_par.v ); % phase difference increment
            in.acc = 3;  % digits of phase difference increment;
            disp( in.PDI );
            
            out = ssfp_fid( in );
            
            xy_TH( i_v, i_Tx ) = abs( out.xy );
                  
        end
    
    end
        
    %% Show results
    
    % colors and line style
    
    col = [ 'b', 'r', 'g', 'k' ];
    
    hold off;
    
    ax = subplot( 1, 1, 1 );
    
    for i_Tx = 1 : n_T2
        
        plot( vs, xy_CM( :, i_Tx ), col( i_Tx ), 'DisplayName', [ 'CM: T1 / T2 = ', num2str( par.T1( i_Tx ) ), ' / ', num2str( par.T2( i_Tx ) ) ] );
        
        legend( 'Interpreter', 'latex' );
        
        if ( i_Tx == 1 )
            
            legend( '-DynamicLegend' );
            hold all;
    
        end
                
        plot( vs, xy_TH( :, i_Tx ), [ col( i_Tx ), '+' ], 'DisplayName', [ 'TH: T1 / T2 = ', num2str( par.T1( i_Tx ) ), ' / ', num2str( par.T2( i_Tx ) ) ] );
        
        legend( 'Interpreter', 'latex' );
        
    end
    
    title( 'Bulk Motion in Unbalanced SSFP', 'Interpreter', 'latex' );
    xlabel( '$v$ [mm/s]', 'Interpreter', 'latex' );
    ylabel( '$\left|m_{xy}\right|$', 'Interpreter', 'latex' );
    xlim( [ 0 par.v_max ] );
    
    width = 9;
    height = 9;
    
    set( gcf, 'Units', 'centimeters' );
    set( gcf, 'Position', [ 0, 0, width, height ] );
    set( gcf, 'Color', 'w' );
    
    delta = 0.01;
    
    ti = ax.TightInset;
    left = ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 1 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];

end

function out = ssfp_fid( in )
% SSFP_FID
% Calculates analytical FID steady state signal (immediately after an RF pulse)
% of unbalanced SSFP with optional RF spoiling.
% (based upon: Ganter C, MRM 55:98-107 (2006))
%
% Ideal spoiling can be generated by setting the phase difference increment PDI == Inf.
%
% Used convention the continued fraction:
% Lambda = b_0 + a_0 / ( b_1 + a_1 / ( b_2 + ...
%
% Input:
%
% Required parameters:
%
% in.FA = flip angle [deg]
% in.TR = repetition time
% in.R1 = longitudinal relaxation rate
% in.R2 = transverse relaxation rate (not required for Ernst formula, i.e. if PDI == Inf)
%
% All these parameters can be either scalars or vectors.
% It is possible to mix (i.e. some parameters as scalars, some as vectors),
% but all vectors must have the same length, since (only) vector elements with
% equal index are combined to a parameter set.
%
% Benefit:
% Avoids multiple function calls in case of independent calculations.
%
% Only the phase difference increment must be a scalar:
%
% in.PDI = RF phase difference increment [deg]
%          (set to Inf for Ernst formula, as mentioned above)
%
% Reason:
% Its value defines the periodicity N, which in turn determines the allocated memory.
%
% Optional parameters:
%
% in.rel_err = relative accuracy of Lambda, default = eps
% in.acc = PDI rounded to 'acc' decimals (to limit periodicity N), default = 2
% in.period = periodicity (sets N by hand)
% in.max_iter = upper limit of iterations, default = 1000
%
% Output:
%
% out.xy = transverse magnetization
% out.Lambda = according to theory
%
% out.iter = performed iterations

def.rel_err = eps;
def.acc = 2;
def.max_iter = 1000;

E1 = exp( - in.R1( : ) .* in.TR( : ) );

if ( in.PDI ~= Inf )
    
    E2 = exp( - in.R2( : ) .* in.TR( : ) );
    
end

ca = cosd( in.FA( : ) );
sa = sind( in.FA( : ) );
deg_2_rad = pi / 180;
y = exp( 1i .* deg_2_rad .* in.PDI );

if ( in.PDI == Inf ) % ideal spoiling == Ernst formula
    
    out.Lambda = 0;
    
    for i = 1 : d.n
        
        X = in.derivatives{ i };
        dX = [ 'd', X ];
        
        out.dLambda.( dX ) = 0;
        
    end
    
else
    
    rel_err = def.rel_err;
    max_iter = def.max_iter;
    
    if ( isfield( in, 'rel_err' ) )
        
        rel_err = in.rel_err;
        
    end
    
    if ( isfield( in, 'max_iter' ) )
        
        max_iter = in.max_iter;
        
    end
    
    if ( isfield( in, 'period' ) )
        
        N = in.period;
        
    else
        
        if ( isfield( in, 'acc' ) )
            
            acc = in.acc;
            
        else
            
            acc = def.acc;
            
        end
        
        if ( acc > 0 )
            
            fak = 10^acc;
            N = fak * 360 / gcd( fak * 360 , round( fak * in.PDI ) );
            
        else
            
            N = 360 / gcd( 360 , round( in.PDI ) );
            
        end
        
    end
    
    if ( N == 1 )
        
        N = 2;
        
    end
    
    tiny_num = min( 1e-30, realmin( 'single' ) );
    
    % exploit periodicity
    
    yn = y.^( 1 : N );
    ynsq = yn .* yn;
    
    aln = 1 - ca .* E1 .* yn;
    ben = ca - E1 .* yn;
    gan = 0.5 .* ( aln - ben );
    
    ga = gan ./ aln;
    bg = ben ./ gan;
    
    E2sqy = E2.^2 .* y;
    
    a0 = ( ga( :, 1 ) + bg( :, 1 ) ) .* E2sqy;
    b0 = - bg( :, 1 ) .* E2sqy;
    
    n = 1 : N;
    n1 = mod( n, N ) + 1;
    
    an = - ga .* ( ga( :, n1 ) + bg( :, n1 ) ) .* ynsq .* E2sqy;
    bn = 1 + ga .* bg( :, n1 ) .* ynsq .* E2sqy;
    
    out.Lambda = b0;
    
    C0 = out.Lambda;
    D1 = bn( :, 1 );
    D1_msk = D1 == 0;
    D1( D1_msk ) = tiny_num;
    
    C1 = bn( :, 1 ) + a0 ./ C0;
    C1_msk = C1 == 0;
    C1( C1_msk ) = tiny_num;
    
    D1 = 1 ./ D1;
    
    dj = C1 .* D1;
    
    out.Lambda = out.Lambda .* dj;
    
    C0 = C1;
    D0 = D1;
    
    j = 1;
    j1 = 2;
    
    out.iter = 0;
    
    while ( max( abs( dj - 1 ) ) > rel_err && out.iter < max_iter )
        
        D1 = bn( :, j1 ) + an( :, j ) .* D0;
        D1_msk = D1 == 0;
        D1( D1_msk ) = tiny_num;
        
        C1 = bn( :, j1 ) + an( :, j ) ./ C0;
        C1_msk = C1 == 0;
        C1( C1_msk ) = tiny_num;
        
        D1 = 1 ./ D1;
        
        dj = C1 .* D1;
        
        out.Lambda = out.Lambda .* dj;
        
        C0 = C1;
        D0 = D1;
        
        j = mod( j, N ) + 1;
        j1 = mod( j, N ) + 1;
        
        out.iter = out.iter + 1;
        
    end
    
end

nom_ = - 1i .* sa .* ( 1 - conj( out.Lambda ) ) .* ( 1 - E1 );
den_ = 1 - E1 .* ca - ( 1 - ca ) .* ( 1 + E1 ) .* real( out.Lambda ) + ( E1 - ca ) .* abs( out.Lambda ).^2;

out.xy = nom_ ./ den_;

end
