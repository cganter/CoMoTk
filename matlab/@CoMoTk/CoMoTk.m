classdef CoMoTk < matlab.mixin.Copyable
    % CoMoTk  Configuration Model Toolkit
    %
    % Implementation of the configuration model (CM), which has
    % been extensively described in a separate manuscript.
    % (link available on: https://github.com/cganter/CoMoTk)
    % 
    % Also check out the specific documentation and example scripts.
    %
    % Requires Matlab R2018b (or later).
    %
    % Carl Ganter 2020
    
    properties ( Constant )
        
        % frequently used parameters
        
        sqrt_2 = sqrt( 2 );
        sqrt_0p5 = sqrt( 0.5 );
        
        % diffusion tensor indices

        Tensor_idx = struct( ...
			'xx',      1, ...
			'yy',      2, ...
			'zz',      3, ...
			'xy',      4, ...
			'xz',      5, ...
			'yz',      6 ...
			);

        % diffusion variants
        
        No_Diff = 1;
        Iso_Diff = 2;
        Tensor_Diff = 3;
        
    end
    
    properties
        
        % function handle for inhomogeneous broadening
        %
        % the format is 'inhomogeneous_decay( t, param )'
        % where 't' is the time and 'param' may contain arbitrary fields with
        % optional parameters (supplied in the equally named argument of
        % the method 'sum')
        % Any user supplied function needs to accept these two parameters
        % (if 'param' is not needed, its value can simply be ignored)
        
        inhomogeneous_decay = [];
        
    end
    
    properties ( Dependent )
        
        % ALL dependent quantities (except for the optional s)
        % have to be set explicitly prior to simulation. Otherwise an error is returned.
        
        %% tissue properties
        % relaxation rates
        
        R1;                                                                % [ 1 / ms ]
        R2;                                                                % [ 1 / ms ]
        
        % diffusion sensitivity (scalar or tensor)
        
        D;                                                                 % [ um^2 / ms ]  i.e. D( Water ) is approximately equal to 3
        
        %% multiple tissues
        % (relative) proton density
        
        pd;
        
        % delta omega (due to chemical shift, e.g. to specify water/fat models)
        
        dom;                                                               % [ rad / ms ]
        
        % transition rates (magnetization transfer and exchange)
        
        k;                                                                 % [ 1 / ms ]
        
        % as unit for bulk motion we assume [ mm / s ]
        % (used in method "time")
        
        %% measurement conditions
        % relative B1
        
        B1;
        
        %% allocated (p,tau) space
        
        n_p;
        d_p;
        
        n_tau;
        d_tau;
                
    end
    
    properties ( SetAccess = private )
        
        b_prep = false;                                                    % is everything prepared?
        
        %% information about the actual state
        
        m = [];                                                            % configuration vector
                                                                           % m = 6d-complex array
                                                                           % dimensions:
                                                                           % 1: magnetization vector components (== 3)
                                                                           % 2: number of tissues (== cm.n_tissues)
                                                                           % 3-5: accumulated moments p(j) (== 2 * cm.n_p(j) + 1)
                                                                           % 6: accumulated transverse time tau (== 2 * cm.n_tau + 1)
                                                                           
        n_tissues = 0;		                                               % number of tissues
        
        d_p_priv = zeros( 1, 3 );                                          % 3-d real, associated small increments (absolute values) 
        n_p_priv = zeros( 1, 3 );                                          % 3-d integer, number of positive gradient moments in each direction
        
        p = { [], [], [] };                                                % associated gradient moments
        
        d_tau_priv = 0;                                                    % real, associated small increment (absolute value)
        n_tau_priv = 0;                                                    % integer, number of positive accumulated times
        
        tau = [];                                                          % associated duration vectors
        
        b_occ = [];                                                        % actually occupied configurations
                                                                           % size( cm.b_occ ) == size( cm.m )
    
        p0 = ones( 1, 3 );                                                 % 3-d index array, corresponding to p == 0
        tau0 = 1;                                                          % index, corresponding to tau == 0
        
        %% mandatory tissue properties
        % general properties
        
        diffusion = CoMoTk.No_Diff;
        
        % cf. dependent counterparts
        
        R1_priv = [];
        R2_priv = [];
        D_priv = [];
           
        %% optional tissue properties and measurement conditions
        % cf. dependent counterparts
        
        pd_priv = 1;
        dom_priv = 0;
        k_priv = [];
        B1_priv = 1;
                
        % auxiliary variable for tensor diffusion
        
        D_vec = [];
        
    end
    
    methods
        
        %% constructor
        
        function cm = CoMoTk ( )
            % constructor (no arguments, no action)
            
        end
        
        %% set / get methods (mandatory variables)
        
        function set.R1 ( cm, R1 )
            % set relaxation rate R1 = 1 / T1
            %
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( min( R1 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R1 );
                cm.pd = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( length( R1 ) ~= cm.n_tissues )
                
                error( 'length( R1 ) must be equal to number of tissues.' );
                
            end
            
            cm.R1_priv = reshape( R1, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.R1 ( cm )
            % return relaxation rate
            
            res = cm.R1_priv;
            
        end
        
        function set.R2 ( cm, R2 )
            % set relaxation rate R2 = 1 / T2
            %
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( min( R2 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R2 );
                cm.pd = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( length( R2 ) ~= cm.n_tissues )
                
                error( 'length( R2 ) must be equal to number of tissues.' );
                
            end
            
            cm.R2_priv = reshape( R2, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.R2 ( cm )
            % return relaxation rate
            
            res = cm.R2_priv;
            
        end
        
        function set.D ( cm, D )
            % set diffusion constant or tensor
            %
            % determined by size of argument:
            % size( D ) == [ 1, n_tissues ]    (isotropic diffusion constant)
            % size( D ) == [ 3, 3, n_tissues ] (anisotropic diffusion tensor)
            %
            % absence of diffusion: just set all elements D( : ) to zero 
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            sz_D = size( D );
            
            if ( sz_D( 1 ) ~= 3 || sz_D( 2 ) ~= 3 )
                
                cm.diffusion = CoMoTk.Iso_Diff;
                
                n_tis = sz_D( 2 );
                
            else
                
                cm.diffusion = CoMoTk.Tensor_Diff;
                
                if ( length( sz_D ) == 3 )
                    
                    n_tis = sz_D( 3 );
                    
                else
                    
                    n_tis = 1;
                    
                end
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = n_tis;
                cm.pd = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( n_tis ~= cm.n_tissues )
                
                error( 'Supplied D does not match number of tissues.' );
                
            end
            
            if ( cm.diffusion == CoMoTk.Tensor_Diff )
                
                cm.D_vec = zeros( 6, cm.n_tissues );
                
                cm.D_vec( CoMoTk.Tensor_idx.xx, : ) = D( 1, 1, : );
                cm.D_vec( CoMoTk.Tensor_idx.yy, : ) = D( 2, 2, : );
                cm.D_vec( CoMoTk.Tensor_idx.zz, : ) = D( 3, 3, : );
                cm.D_vec( CoMoTk.Tensor_idx.xy, : ) = 2 .* D( 1, 2, : );
                cm.D_vec( CoMoTk.Tensor_idx.xz, : ) = 2 .* D( 1, 3, : );
                cm.D_vec( CoMoTk.Tensor_idx.yz, : ) = 2 .* D( 2, 3, : );
                
            end
            
            cm.D_priv = D;
            
            if ( max( abs( cm.D( : ) ) ) == 0 )
                
                cm.diffusion = CoMoTk.No_Diff;
                
            end
            
        end
        
        function res = get.D ( cm )
            % return diffusion constant or tensor
            
            res = cm.D_priv;
            
        end
        
        function set.d_p ( cm, d_p )
            % set positive gradient moment increment in each direction
            %
            % IN
            %
            % non-negative vector of length 3 
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( length( d_p ) ~= 3 || min( d_p ) < 0 )
                
                error( 'Invalid setting of d_p.' );
                
            end
                        
            cm.d_p_priv = d_p( : );
            
        end
        
        function res = get.d_p ( cm )
            
            res = cm.d_p_priv;
            
        end        
        
        function set.n_p ( cm, n_p )
            % set number of positive gradient moments in each direction
            %
            % IN
            %
            % non-negative vector of length 3 
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( length( n_p ) ~= 3 || min( n_p ) < 0 )
                
                error( 'Invalid setting of n_p.' );
                
            end
                        
            cm.n_p_priv = round( n_p( : )' );
            
            cm.p0 = cm.n_p_priv + 1;
            
        end
        
        function res = get.n_p ( cm )
            
            res = cm.n_p_priv;
            
        end        
        
        function set.d_tau ( cm, d_tau )
            % set positive gradient moment increment in each direction
            %
            % IN
            %
            % positive number 
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( d_tau <= 0 )
                
                error( 'Invalid setting of d_tau.' );
                
            end
                        
            cm.d_tau_priv = d_tau;
            
        end
        
        function res = get.d_tau ( cm )
            
            res = cm.d_tau_priv;
            
        end        
           
        function set.n_tau ( cm, n_tau )
            % set number of positive accumulated times in each direction
            %
            % IN
            %
            % positive number
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( n_tau < 1 )
                
                error( 'Invalid setting of n_tau.' );
                
            end
                        
            cm.n_tau_priv = round( n_tau );
            
            cm.tau0 = cm.n_tau_priv + 1;
            
        end
        
        function res = get.n_tau ( cm )
            
            res = cm.n_tau_priv;
            
        end        
        
        function [ occ, fra ] = n_conf( cm )
        % returns information about stored configurations in (p,tau) space
        %
        % OUT
        %
        % occ : total number of stored configurations
        % fra : fractional occupation in each of the four CS directions
        %
        % remark, regarding the information returned in 'fra':
        %
        % any j = 1, ..., 4, the value occ( j ) takes a value from the
        % interval [ 0, 1 ].
        % Its definition is explained for the tau dimension (j == 4):
        % fra( j ) = max( abs( index of stored tau ) ) / n_tau
        %
        % A small value at the end of the calculation indicates a
        % potentially too cautious setting of n_tau.
        %
        % A value equal to 1 shows that the allocated space has been fully
        % used in this direction. This may indicate that more space should
        % be allocated (depending on the magnitude of the discarded
        % configurations).
        %
        % The value NaN is assigned to unset directions with n_p( j ) == 0
            
            b = any( cm.b_occ( :, 1, :, :, :, : ), 1 );
            
            occ = sum( b, 'all' );
            
            fra = NaN( 4, 1 );

            % p direction
            
            for j = 1 : 3
                
                if ( cm.n_p( j ) ~= 0 )
                
                    b_ = any( b, 2 + [ 1 : j-1, j+1 : 4 ] );
                    rng = - cm.n_p( j ) : cm.n_p( j );
                    fra( j ) = max( abs( rng( b_( : ) ) ) ) ./ cm.n_p( j );
            
                end
                
            end
            
            % tau direction
            
            b_ = any( b, 3 : 5 );
            rng = - cm.n_tau : cm.n_tau;
            fra( 4 ) = max( abs( rng( b_( : ) ) ) ) ./ cm.n_tau;
            
        end
        
        %% set / get methods (optional variables)
        
        % relative B1+
        
        function set.B1 ( cm, B1 )
            % set relative B1+
            %
            % set to 1, if not explicitly specified            
            
            if ( cm.b_prep )
                
                error( 'B1 can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( B1 < 0 )
                
                error( 'B1 must not be negative.' );
                
            end
            
            cm.B1_priv = B1;
            
        end
        
        function res = get.B1 ( cm )
            % get relative B1+
            
            res = cm.B1_priv;
            
        end
        
        % (relative) proton density of complex tissues
        % (default == 1)
        
        function set.pd ( cm, pd )
            % set the (relative) proton density
            %
            % default = 1 (only, if cm.n_tissues == 1)
            % mandatory parameter, if cm.n_tissues > 1
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( pd );
                
            elseif ( length( pd ) ~= cm.n_tissues )
                
                error( 'length( pd ) must be equal to number of tissues.' );
                
            end
            
            cm.pd_priv = reshape( pd, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.pd ( cm )
            % return the (relative) proton density
            
            res = cm.pd_priv;
            
        end
        
        % frequency offset due to chemical shift
        
        function set.dom ( cm, dom )
            % relative resonance frequencies of subtissues 
            % (e.g. due to chemical shift)
            %
            % length must match the number of subtissues (e.g. for fat models)            
            % set to 0, if not explicitly specified
                        
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( dom );
                
            elseif ( length( dom ) ~= cm.n_tissues )
                
                error( 'length( dom ) must be equal to number of tissues.' );
                
            end
            
            cm.dom_priv = reshape( dom, [ cm.n_tissues, 1 ] );
            
        end
        
        function res = get.dom ( cm )
            % return relative resonance frequencies
            
            res = cm.dom_priv;
            
        end
                
        % transition rate matrix (ME/MT)
        
        function set.k ( cm, k )
            % transition rates for magnetization transfer/exchange
            %
            % the following size is required:
            % size( k ) == [ n_tissues, n_tissues ]
            % set to [], if not set (non-exchanging tissues)
            % 
            % Note: 
            % Diagonal elements are calculated according to sum rule
            % (particle conservation). Supplied diagonal values are ignored. 
            
            if ( cm.n_tissues == 0 )
                
                error( 'initialize tissues before setting the transition rates');
                
            end 
            
            sz_k = size( k );
            
            if ( length( sz_k ) ~= 2 || sz_k( 1 ) ~= cm.n_tissues || sz_k( 2 ) ~= cm.n_tissues )
                
                error( 'size( k ) ~= [ n_tissues, n_tissues ].' );
                
            end
                        
            % establish sum rule (particle conservation)

            k = k - diag( diag( k ) );  % first, set the diagonal to zero

            cm.k_priv = k - diag( sum( k ) );   % now we assign the negative row sum to the diagonal
                                                % and store the result
                            
        end
        
        function res = get.k( cm )
            % get transition rates (== [], if not set)
            
            res = cm.k_priv;
            
        end
        
        %% spin dynamics
        
        function RF ( cm, param )
            % Applies an instantaneous RF pulse
            %
            % IN
            %
            % param : structure with the info about the RF pulse,
            % specifically:
            %
            % param.FlipAngle : flip angle [rad] (mandatory parameter)
            %
            % param.Phase     : phase [rad] (mandatory parameter)
            % 
            % param.MT_sat           : parameters for MT saturation pulse
            % (instead of instantaneous rotation) 
            % For proper use see example script and associated passage in the manuscript.
            
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            % check for minimal information
            
            if ( ~isfield( param, 'FlipAngle' ) )
                
                error( 'No flip angle supplied.' );
                
            end
            
            if ( ~isfield( param, 'Phase' ) )
                
                error( 'No phase supplied.' );
                
            end
            
            % this is our actual flip angle
            
            ActualFlipAngle = cm.B1 * param.FlipAngle;
            
            % no action required, if no rotation is performed
            
            if ( mod( ActualFlipAngle, 2 * pi ) == 0 )
                
                return;
                
            end
            
            % calculate the rotation matrix
            
            c_al = cos( ActualFlipAngle );
            s_al = sin( ActualFlipAngle );
            
            ei_ph = cos( param.Phase ) + 1i * sin( param.Phase );
            
            RotMat = zeros( 3 );
            
            RotMat( 1, 1 ) = 0.5 * ( 1 + c_al );
            RotMat( 1, 2 ) = - 1i * s_al * ei_ph * CoMoTk.sqrt_0p5;
            RotMat( 1, 3 ) = 0.5 * ( 1 - c_al ) * ei_ph^2;
            RotMat( 2, 1 ) = - conj( RotMat( 1, 2 ) );
            RotMat( 2, 2 ) = c_al;
            RotMat( 2, 3 ) = - RotMat( 1, 2 );
            RotMat( 3, 1 ) = conj( RotMat( 1, 3 ) );
            RotMat( 3, 2 ) = - RotMat( 2, 1 );
            RotMat( 3, 3 ) = RotMat( 1, 1 );
            
            % check for magnetization transfer (MT)
            
            b_MT = false( 1, cm.n_tissues );
            
            if ( isfield( param, 'MT_sat' ) )
            
                if ( length( param.MT_sat ) ~= cm.n_tissues )
                    
                    error( 'wrong size of saturation due to magnetization transfer' );
                    
                end
                
                % a value > 0 is treated as a bound pool
                
                b_MT( : ) = param.MT_sat > 0;
                                
            end
            
            n_MT = sum( b_MT );
            
            if ( n_MT > 0 )
                
                if ( cm.n_tissues > n_MT ) % should be always the case..
                
                    cm.m( cm.b_occ & ~b_MT ) = ...
                        RotMat * ...
                        reshape( cm.m( cm.b_occ & ~b_MT ), 3, [] );
                    
                end
 
                % MT saturation affects only the longitudinal component

                long_ = [ false; true; false ];
                
                cm.m( long_ & b_MT & cm.b_occ ) = ...
                    reshape( param.MT_sat( b_MT ), [ 1, n_MT ] ) .* ...
                    reshape( cm.m( long_ & b_MT & cm.b_occ ), 1, n_MT, [] );
                
                % since transverse magnetization is never generated, we
                % don't need to set it to zero explicitly

            else
            
                %% update the configuration cm.m
                
                cm.m( cm.b_occ ) = ...
                    RotMat * reshape( cm.m( cm.b_occ ), 3, [] );
                
            end
                               
        end
        
        function time ( cm, param )
            % Applies time interval
            %
            % IN
            %
            % param : structure with the info about the RF pulse,
            % specifically:
            %
            % param.tau    : duration of interval [ ms ]
            % (mandatory parameter)
            % 
            % param.p      : zero-order gradient moment [ rad / um ]
            % (set to zero, if absent)
            % 
            % param.s and param.S : information about gradient shape
            % (relevant for diffusion only, otherwise ignored)
            % 
            % param.v      : (constant) bulk motion velocity vector during
            % interval [ mm / s ] (optional)
            %
            % param.x      : location (vector) at start of interval [ um ]
            % (required in case of bulk motion)
            %
            % param.p1     : first order gradient moment [ rad * ms / um ]
            % (optional, only relevant in case of bulk motion)za

            %% prepare everything, if not done so yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            %% check for duration and zero-order gradient moment
            
            % duration (mandatory)
            
            if ( isfield( param, 'tau' ) )
                
                tau_ = param.tau;
                
            else
                
                error( 'Time interval requires duration.' );
                
            end
            
            % zero order gradient moment
            
            if ( isfield( param, 'p' ) )
                
                p_ = param.p( : );
                
            else
                
                p_ = zeros( 3, 1 );
                
            end
 
            %% prepare relaxation
            
            E = zeros( 3, cm.n_tissues );
            
            E( 1, : ) = exp( - cm.R2 .* tau_ );
            E( 2, : ) = exp( - cm.R1 .* tau_ );
            E( 3, : ) = E( 1, : );
            
            %% prepare diffusion damping
            
            % first, set gradient shape
            
            if ( cm.diffusion ~= CoMoTk.No_Diff )
                
                if ( isfield( param, 's' ) )
                                        
                    s_ = param.s( : );
                    
                else
                    
                    s_ = p_;
                    
                end
                
                if ( isfield( param, 'S' ) )
                    
                    if ( cm.diffusion == CoMoTk.Iso_Diff )
                        
                        S_vec = param.S;
                        
                    elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                        
                        S_vec = zeros( 6, 1 );
                        
                        S_vec( CoMoTk.Tensor_idx.xx ) = param.S( 1, 1 );
                        S_vec( CoMoTk.Tensor_idx.yy ) = param.S( 2, 2 );
                        S_vec( CoMoTk.Tensor_idx.zz ) = param.S( 3, 3 );
                        S_vec( CoMoTk.Tensor_idx.xy ) = param.S( 1, 2 );
                        S_vec( CoMoTk.Tensor_idx.xz ) = param.S( 1, 3 );
                        S_vec( CoMoTk.Tensor_idx.yz ) = param.S( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                                        
                else
                    
                    if ( cm.diffusion == CoMoTk.Iso_Diff )
                        
                        S_vec = ( p_' * p_ ) / 3;
                        
                    elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                        
                        S_vec = zeros( 6, 1 );
                        
                        S_tmp = ( p_ * p_' ) ./ 3;
                        
                        S_vec( CoMoTk.Tensor_idx.xx ) = S_tmp( 1, 1 );
                        S_vec( CoMoTk.Tensor_idx.yy ) = S_tmp( 2, 2 );
                        S_vec( CoMoTk.Tensor_idx.zz ) = S_tmp( 3, 3 );
                        S_vec( CoMoTk.Tensor_idx.xy ) = S_tmp( 1, 2 );
                        S_vec( CoMoTk.Tensor_idx.xz ) = S_tmp( 1, 3 );
                        S_vec( CoMoTk.Tensor_idx.yz ) = S_tmp( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                    
                end
                
                S_ = S_vec( : );

            end
            
            % now calculate the damping factors

            Diff_Damp = [];
            
            if ( cm.diffusion ~= CoMoTk.No_Diff )
                
                if ( cm.diffusion == CoMoTk.Iso_Diff )
                    
                    tau_D = tau_ .* cm.D;
                    
                    pn2 = 0;
                    pns = 0;
                    
                    for i = 1 : 3
                    
                        pn2 = pn2 + cm.p{ i } .* cm.p{ i };
                        pns = pns + cm.p{ i } .* s_( i );
                        
                    end
                                        
                    Diff_Damp_exp = zeros( [ 3, 1, 2 * cm.n_p + 1 ] );
                    
                    Diff_Damp_exp( 1, 1, :, :, : ) = - ( pn2 + pns + S_ );
                    Diff_Damp_exp( 2, 1, :, :, : ) = - pn2;
                    Diff_Damp_exp( 3, 1, :, :, : ) = - ( pn2 - pns + S_ );
                    
                    % finalize exponent
                    
                    Diff_Damp_exp = tau_D .* Diff_Damp_exp + zeros( [ ones( 1, 5 ), 2 * cm.n_tau + 1 ] );
                    
                elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                    
                    pnDpn = 0;
                    pnDs = 0;
                    
                    for i = 1 : 3

                        for j = 1 : 3
                        
                            pnDpn = pnDpn + cm.p{ i } .* cm.D( i, j ) .* cm.p{ j };
                            pnDs = pnDs + cm.p{ i } .* cm.D( i, j ) .* s_( j );
                       
                        end
                        
                    end
                                        
                    DS = reshape( sum( cm.D_vec .* S_, 1 ), [ 1, cm.n_tissues ] );
                    
                    Diff_Damp_exp = zeros( [ 3, 1, 2 * cm.n_p + 1 ] );
                    
                    Diff_Damp_exp( 1, 1, :, :, : ) = - ( pnDpn + pnDs + DS );
                    Diff_Damp_exp( 2, 1, :, :, : ) = - pnDpn;
                    Diff_Damp_exp( 3, 1, :, :, : ) = - ( pnDpn - pnDs + DS );
                                        
                    % finalize exponent
                    
                    Diff_Damp_exp = tau_ .* Diff_Damp_exp + zeros( [ ones( 1, 5 ), 2 * cm.n_tau + 1 ] );
                    
                else
                    
                    error( 'Unexpected value for ''cm.diffusion''.' );
                    
                end
                    
                % update diffusion damping factors
                    
                Diff_Damp = exp( Diff_Damp_exp( cm.b_occ ) );
                    
            end
            
            %% initialize chemical shift
            
            cs_tau = ones( 3, cm.n_tissues );
            
            cs_tau( 1, : ) = exp( - 1i .* cm.dom .* tau_ );
            cs_tau( 3, : ) = conj( cs_tau( 1, : ) );    
           
            %% execute: relaxation, chemical shift, diffusion damping, chemical exchange
            
            if ( isempty( cm.k ) )              % no magnetization exchange
                
                %% update the configuration cm.m
                
                % chemical shift and relaxation                
                
                cm.m( cm.b_occ ) = ( cs_tau .* E ) .* reshape( cm.m( cm.b_occ ), 3, cm.n_tissues, [] );
                            
                % diffusion damping, if relevant
                
                if ( ~isempty( Diff_Damp ) )

                    cm.m( cm.b_occ ) = Diff_Damp .* cm.m( cm.b_occ );
                    
                end
                
            else                              % with magnetization exchange

                b_inf_R2 = isinf( cm.R2 );
                
                if ( any( b_inf_R2 ) ) 
                    
                    % In case of MT, there is no coupling of transverse
                    % magnetization to bound pool(s).
                    %
                    % For the (rare?) case than there is more than one
                    % mobile compartment, we allow for coupling between
                    % them, though. To this end, we define a transition matrix
                    % for the free compartments only and modify the sum
                    % rule accordingly (otherwise the diagonal elements
                    % would leak into the bound pool, which is usually
                    % assumed as uncoupled with respect to the transverse
                    % magnetization).
                    %
                    % NOTE: 
                    % Use with care, since I am not sure, whether this
                    % approach is sound.
                    
                    k_ = cm.k;
                    
                    % eliminate coupling to bound pool(s)
                    
                    k_( :, b_inf_R2 ) = 0;
                    k_( b_inf_R2, : ) = 0;
                    
                    % modify sum rule for transverse components
                    
                    for j = 1 : cm.n_tissues
                        
                        k_( j, j ) = 0;
                        
                        for i = 1 : cm.n_tissues
                            
                            if ( i ~= j )
                                
                                k_( j, j ) = k_( j, j ) - k_( i, j );
                                
                            end
                            
                        end
                        
                    end

                    Y_1 = k_ - diag( 1i .* cm.dom( : ) + cm.R2( : ) );
                    
                    Y_1( :, b_inf_R2 ) = 0;
                    Y_1( b_inf_R2, : ) = 0;
                    
                else
                    
                    Y_1 = cm.k - diag( 1i .* cm.dom( : ) + cm.R2( : ) );

                end
                
                Y_0 = cm.k - diag( cm.R1( : ) );

                exp_Y_tau_1 = expm( Y_1 .* tau_ );
                exp_Y_tau_0 = expm( Y_0 .* tau_ );

                y = cm.R1( : ) .* cm.pd( : );
                
                % relaxation
                  
                co_ = eye( 3 ) ~= 0;    % to access components via logical indexing
                
                cm.m( co_( :, 1 ) & cm.b_occ ) = ...
                    exp_Y_tau_1 * ...
                    reshape( cm.m( co_( :, 1 ) & cm.b_occ ), cm.n_tissues, [] );
                
                cm.m( co_( :, 2 ) & cm.b_occ ) = ...
                    exp_Y_tau_0 * ...
                    reshape( cm.m( co_( :, 2 ) & cm.b_occ ), cm.n_tissues, [] );
                
                cm.m( co_( :, 3 ) & cm.b_occ ) = ...
                    conj( exp_Y_tau_1 ) * ...
                    reshape( cm.m( co_( :, 3 ) & cm.b_occ ), cm.n_tissues, [] );
                
            end
                       
            % bulk motion due to constant velocity during the interval
            
            if ( isfield( param, 'v' ) )
                
                % "v" is the actual velocity vector, units are [ mm / s ]
                % (it does not matter, whether it is a row or column vector)

                if ( ~isfield( param, 'x' ) )
                    
                    error( 'The location should be supplied in case of bulk motion.' );
                    
                end                
                
                if ( isfield( param, 'p1' ) )
                    
                    % If the first order gradient is supplied, we use it.
                    % Unit of p1: [ rad * ms / um ]

                    dth = - param.p1( : )' * param.v( : );

                else
                    
                    % 0therwise, we assume a constant gradient during the
                    % interval. In that case, the angle due to dephasing is just
                    % - (tau/2) * (p * v)
                    
                    % NOTE: nonzero velocities are currently not taken into account, when calculating derivatives!
                    
                    dth = - ( 0.5 * tau_ ) * ( p_' * param.v( : ) );
                    
                end
                
                dth = dth - p_' * param.x( : );
                
                ei_dth = exp( 1i * dth );
                
                co_ = eye( 3 ) ~= 0;    % to access components via logical indexing
     
                cm.m( co_( :, 1 ) & cm.b_occ ) = ...
                    ei_dth .* cm.m( co_( :, 1 ) & cm.b_occ );
                
                cm.m( co_( :, 3 ) & cm.b_occ ) = ...
                    conj( ei_dth ) .* cm.m( co_( :, 3 ) & cm.b_occ );
                
            end
            
            %% execute shift in configuration space
            
            % we do this for each direction
            
            p_shift = zeros( 3, 1 );
            
            for i = 1 : 3
                
                if ( p_( i ) ~= 0 )
                    
                    if ( cm.d_p( i ) == 0 )
                        
                        error( 'Nonzero moment in unexpected direction.' );
                        
                    end
                    
                    p_shift( i ) = round( p_( i ) / cm.d_p( i ) );
                    
                end

                if ( p_shift( i ) ~= 0 )
                    
                    cm.m( 1, :, :, :, :, : ) = circshift( ...
                        cm.m( 1, :, :, :, :, : ), p_shift( i ), i + 2 );
            
                    cm.m( 3, :, :, :, :, : ) = circshift( ...
                        cm.m( 3, :, :, :, :, : ), - p_shift( i ), i + 2 );
            
                    cm.b_occ( 1, :, :, :, :, : ) = circshift( ...
                        cm.b_occ( 1, :, :, :, :, : ), p_shift( i ), i + 2 );
            
                    cm.b_occ( 3, :, :, :, :, : ) = circshift( ...
                        cm.b_occ( 3, :, :, :, :, : ), - p_shift( i ), i + 2 );
            
                end
                
            end
            
            tau_shift = round( tau_ / cm.d_tau );
            
            if ( tau_shift ~= 0 )
                
                cm.m( 1, :, :, :, :, : ) = circshift( ...
                    cm.m( 1, :, :, :, :, : ), tau_shift, 6 );
                
                cm.m( 3, :, :, :, :, : ) = circshift( ...
                    cm.m( 3, :, :, :, :, : ), - tau_shift, 6 );
                
                cm.b_occ( 1, :, :, :, :, : ) = circshift( ...
                    cm.b_occ( 1, :, :, :, :, : ), tau_shift, 6 );
                
                cm.b_occ( 3, :, :, :, :, : ) = circshift( ...
                    cm.b_occ( 3, :, :, :, :, : ), - tau_shift, 6 );
                
            end
            
            % now we must invalidate unoccupied locations for each
            % direction
            
            if ( p_shift( 1 ) > 0 )
               
                cm.m( 1, :, 1 : p_shift( 1 ), :, :, : ) = 0;
                cm.m( 3, :, end - p_shift( 1 ) + 1 : end, :, :, : ) = 0;

                cm.b_occ( 1, :, 1 : p_shift( 1 ), :, :, : ) = false;
                cm.b_occ( 3, :, end - p_shift( 1 ) + 1 : end, :, :, : ) = false;
                
            elseif ( p_shift( 1 ) < 0 )
                
                cm.m( 3, :, 1 : - p_shift( 1 ), :, :, : ) = 0;
                cm.m( 1, :, end + p_shift( 1 ) + 1 : end, :, :, : ) = 0;

                cm.b_occ( 3, :, 1 : - p_shift( 1 ), :, :, : ) = false;
                cm.b_occ( 1, :, end + p_shift( 1 ) + 1 : end, :, :, : ) = false;
                
            end
            
            if ( p_shift( 2 ) > 0 )
               
                cm.m( 1, :, :, 1 : p_shift( 2 ), :, : ) = 0;
                cm.m( 3, :, :, end - p_shift( 2 ) + 1 : end, :, : ) = 0;

                cm.b_occ( 1, :, :, 1 : p_shift( 2 ), :, : ) = false;
                cm.b_occ( 3, :, :, end - p_shift( 2 ) + 1 : end, :, : ) = false;
                
            elseif ( p_shift( 2 ) < 0 )
                
                cm.m( 3, :, :, 1 : - p_shift( 2 ), :, : ) = 0;
                cm.m( 1, :, :, end + p_shift( 2 ) + 1 : end, :, : ) = 0;

                cm.b_occ( 3, :, :, 1 : - p_shift( 2 ), :, : ) = false;
                cm.b_occ( 1, :, :, end + p_shift( 2 ) + 1 : end, :, : ) = false;
                
            end
            
            if ( p_shift( 3 ) > 0 )
               
                cm.m( 1, :, :, :, 1 : p_shift( 3 ), : ) = 0;
                cm.m( 3, :, :, :, end - p_shift( 3 ) + 1 : end, : ) = 0;

                cm.b_occ( 1, :, :, :, 1 : p_shift( 3 ), : ) = false;
                cm.b_occ( 3, :, :, :, end - p_shift( 3 ) + 1 : end, : ) = false;
                
            elseif ( p_shift( 3 ) < 0 )
                
                cm.m( 3, :, :, :, 1 : - p_shift( 3 ), : ) = 0;
                cm.m( 1, :, :, :, end + p_shift( 3 ) + 1 : end, : ) = 0;

                cm.b_occ( 3, :, :, :, 1 : - p_shift( 3 ), : ) = false;
                cm.b_occ( 1, :, :, :, end + p_shift( 3 ) + 1 : end, : ) = false;
                
            end
            
            if ( tau_shift > 0 )
               
                cm.m( 1, :, :, :, :, 1 : tau_shift ) = 0;
                cm.m( 3, :, :, :, :, end - tau_shift + 1 : end ) = 0;

                cm.b_occ( 1, :, :, :, :, 1 : tau_shift ) = false;
                cm.b_occ( 3, :, :, :, :, end - tau_shift + 1 : end ) = false;
                
            elseif ( tau_shift < 0 )
                
                cm.m( 3, :, :, :, :, 1 : - tau_shift ) = 0;
                cm.m( 1, :, :, :, :, end + tau_shift + 1 : end ) = 0;

                cm.b_occ( 3, :, :, :, :, 1 : - tau_shift ) = false;
                cm.b_occ( 1, :, :, :, :, end + tau_shift + 1 : end ) = false;
                
            end
            
            % to make life easier, we do not distinguish the occupation
            % of different components
            
            b_ = any( cm.b_occ, 1 );
            
            sz_ = size( cm.b_occ );
            
            cm.b_occ = reshape( repmat( reshape( b_, 1, [] ), [ sz_( 1 ), 1 ] ), sz_ );
            
            %% repolarization
            
            if ( isempty( cm.k ) )              % no magnetization exchange
                                
                cm.m( 2, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) = ...
                    cm.m( 2, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) + ...
                    cm.pd .* ( 1 - E( 2, : ) );
                
            else                              % with magnetization exchange

                % (here, the relative proton densities come into play)
                
                cm.m( 2, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) = ...
                    cm.m( 2, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) - ...
                    reshape( Y_0 \ ( ( eye( cm.n_tissues ) - exp_Y_tau_0 ) * y ), [ 1, cm.n_tissues ] );
                
            end            
            
        end
        
        function spoiler ( cm )
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            % set all transverse magnetization to zero (instantaneously)
           
            cm.m( [ 1, 3 ], :, :, :, :, : ) = 0;
                        
        end
        
        %% get results
        
        function b_n = find ( cm, lambda, n )
            % search for specific configurations
            %
            % IN
            %
            % lambda : array of CM dimension indices (== mandatory
            % parameter in method 'time')
            % n      : array of configuration orders
            %
            % OUT
            %
            % b_n : location of stored configurations, which match ANY of
            % the supplied condition pairs lambda(i), n(i)
            
            if ( length( lambda ) ~= length( n ) )
                
                error( 'lambda and n must have same length.' );
                
            end
            
            % check, which of the supplied lambda are actually alive (i.e. stored)
            
            b_lambda_2d = cm.b_lambda & ( lambda( : ) == cm.lambda );  % size( b_lambda_2d ) == [ length( lambda ), cm.alloc_d ]
            
            % non-existing lambda must not have non-zero configuration orders
            
            b_lambda_found = any( b_lambda_2d, 2 );            % size( b_lambda_found ) == [ length( lambda ), 1 ]
            
            if ( any( n( ~b_lambda_found ) ) )
                
                b_n = [];
                return;
                
            end

            % n == 0 is always true for not (yet) occupied time index lambda
            
            if ( sum( ~b_lambda_found ) > 0 )

                b_n = cm.b_n;
                return;
                
            end
            
            n_lambda_found = sum( b_lambda_found );
            
            % for the rest, determine the indices in storage
            
            [ idx_lambda, ~ ] = find( b_lambda_2d' );
            
            % we can now determine b_n
            
            b_n = cm.b_n & any( cm.n( :, idx_lambda ) == reshape( n( b_lambda_found ), [ 1, n_lambda_found ] ), 2 );
            
        end
        
        function res = sum ( cm, param )            
            % return the weighted CM sum and calculated derivatives
            %
            % IN
            %
            % structure param (if empty, complete sum over all stored
            % configurations is returned)
            %
            % param.b_occ : over which configurations shall the sum be taken?
            % (optional)
            % (if field is missing, the sum is over all stored
            % configurations)
            %
            % param.omega : bulk off-resonance [ rad / ms ] (optional)
            % (default == 0, if field is missing)
            % may be an array
            %
            % param.x : position vector [ um ] (optional)
            % (== zero vector, if field is missing)
            % may have a second dimension
            %
            % decay due to inhomogeneous broadening is taken into account,
            % if cm.inhomogeneous_decay is nonempty
            % (more detailed information about the format of this function
            % can be found above at the associated property
            %
            % OUT
            %
            % res.xy : transverse component (conventional transverse
            % magnetization)
            %
            % res.z  : longitudinal component
            %
            % res.dxy.dX : derivative of transverse component with respect
            % to X
            %
            % res.dx.dX : same for longitudinal component
            
            if ( ~isempty( param ) && isfield( param, 'b_occ' ) )
                
                b_occ_ = param.b_occ;
                
            else
                
                b_occ_ = cm.b_occ;
                
            end
            
            % restriction to configuration space
            
            b_occ_pure = b_occ_( 1, 1, :, :, :, : );
               
            % initialize weigthing factor
            % could also be used for external weighting (currently not
            % implemented)
            % final expected size (if not equal to one) 
            % == [ 1, n_occ, n_om, n_x ]
            
            w_n_ = 1;

            % off-resonance(s)
            
            if ( ~isempty( param ) && isfield( param, 'omega' ) )
                
                n_om = length( param.omega ); 
                omega_ = reshape( param.omega, [ 1, 1, n_om ] );
                
                % update weights
                % size( w_n_ ) == [ 1, n_occ, n_om ]
                
                tau_ = cm.tau + zeros( size( b_occ_pure ) );
                
                w_n_ = w_n_ .* exp( - 1i .* omega_ .* ...
                    reshape( tau_( b_occ_pure ), 1, [] ) );
            
            end

            % location(s)
            
            if ( ~isempty( param ) && isfield( param, 'x' ) )
                
                n_x = size( param.x, 2 );
                x_ = reshape( param.x, [ 3, ones( 1, 5 ), n_x ] );
                                
                % update weights
                % size( w_n_ ) == [ 1, n_occ, n_om, n_x ]
                
                px_ = 0;
                
                for i = 1 : 3
                    
                    px_ = px_ + cm.p{ i } .* x_( i, 1, 1, 1, 1, 1, : );
                    
                end
                
                px_ = px_ + zeros( [ ones( 1, 5 ), size( b_occ_pure, 6 ) ] );
                
                n_occ = sum( b_occ_pure( : ) );
                
                px_ = reshape( px_( b_occ_pure & true( [ ones( 1, 6 ), n_x ] ) ), [ 1, n_occ, 1, n_x ] );
                
                w_n_ = w_n_ .* exp ( - 1i .* px_ );
                
            end
            
            % effects due to inhomogeneous broadening (if applicable, most commonly due to R2p)
            
            if ( ~isempty( cm.inhomogeneous_decay ) )
                
                w_n_ = w_n_ .* cm.inhomogeneous_decay( cm.tau( b_occ_pure ), param );
                
            end
            
            % calculate the isochromat
            
            b_10 = [ true; true; false ];
            
            tmp_ = reshape( ...
                sum( ...
                reshape( cm.m( b_10 & b_occ_ ), 2, cm.n_tissues, [] ) .* ...
                reshape( w_n_, [ 1, size( w_n_ ) ] ), 2 : 3 ), ...
                [ 2, size( w_n_, 3 : 4 ) ] );
            
            res.xy = CoMoTk.sqrt_2 .* tmp_( 1, :, : );
            res.z = tmp_( 2, :, : );
            
        end
        
        %% Lorentz function
        
        function res = lorentz_decay ( cm, tau, param )
            % calculates Lorentz decay exp( - tau * R2p )
            %
            % param.R2p : relaxation rate R2^\prime [ 1 / ms ] (mandatory parameter)
           
            if ( isfield( param, 'R2p' ) )
                
                R2p = param.R2p;
                
                if ( length( R2p ) ~= cm.n_tissues && length( R2p ) ~= 1 )
                    
                    error( 'length( R2p ) must be equal to the number of tissues or one.' );
                    
                end
                
                if ( length( R2p ) == cm.n_tissues )
                    
                    R2p = reshape( R2p, [ cm.n_tissues, 1 ] );
                    
                else
                    
                    R2p = repmat( R2p, [ cm.n_tissues, 1 ] );
                    
                end
                                
                res = exp( - R2p .* reshape( abs( tau ), [ 1, length( tau ) ] ) );
                               
            else
                
                error( 'R2p not specified in lorentz_decay' );
                
            end
            
        end
        
        function res = spherical_decay ( cm, tau, param )
            % decay due to spherical frequency distribution
            %
            % param.bw_sphere : max. off-resonance frequency magnitude 
            % [ rad / ms ] (mandatory parameter)
            
            if ( isfield( param, 'bw_sphere' ) )
                
                sigma = 0.5 .* param.bw_sphere;
                
                if ( length( sigma ) ~= cm.n_tissues && length( sigma ) ~= 1 )
                    
                    error( 'length( bw_sphere ) must be equal to the number of tissues or one.' );
                    
                end
                
                if ( length( sigma ) == cm.n_tissues )
                    
                    sigma = reshape( sigma, [ cm.n_tissues, 1 ] );
                    
                else
                    
                    sigma = repmat( sigma, [ cm.n_tissues, 1 ] );
                    
                end
                                    
                st_ = sigma .* reshape( tau, [ 1, length( tau ) ] );
                    
                res = besselj( 1, 2 .* st_ ) ./ ( st_ );

                res( st_ == 0 ) = 1;
                    
            else
                
                error( 'bw_sphere not specified in spherical_decay' );
                
            end
            
        end
                
    end
    
    methods ( Access = private )
        % several helper routines

        function check ( cm )
            % Check, whether all mandatory tissue parameters have been set.
            
            if ( isempty( cm.R1 ) || ...
                    isempty( cm.R2 ) || ...
                    isempty( cm.D ) )
                
                error( 'Incomplete specification of mandatory tissue parameters.' );
                
            end
            
            if ( cm.d_tau == 0 )
                
                error( 'd_tau must be positive.' );
                
            end                        
            
            % check proton densities
            
            if ( size( cm.pd, 2 ) ~= cm.n_tissues )
                
                error( 'Proton density not set properly.' );
                
            end
            
            % check correct long time limit of magnetization transfer/exchange
            % (particle conservation is properly incorporated by set.k) 
            
            if ( ~isempty( cm.k ) )
               
                % cm.pd must lie in the null space of cm.k:
                
                k_pd = cm.k * cm.pd( : );
                
                if ( max( abs( k_pd ) > 10 * eps ) )
                    
                    % Maybe there is an error, maybe it is just numerical
                    % inaccuracy. Provide an informal message to the user and let
                    % the program continue.
                    
                    fprintf( 1, 'Check, whether cm.pd lies in null space of cm.k\n' );
                    fprintf( 1, 'Found result for cm.k * cm.pd =\n' );
                    disp( k_pd );
                    
                end
                
            end
            
        end
        
        function prep ( cm )
            % Initializes magnetization and other stuff
            
            % check, whether required parameters have been set properly
            
            cm.check;
            
            % no diffusion AND magnetization exchange
            
            if ( cm.diffusion ~= CoMoTk.No_Diff && ~isempty( cm.k ) )
                
                error( 'Diffusion AND magnetization exchange not supported.')
                
            end
            
            % allocate space set initial values
            
            cm.m = zeros( [ 3, cm.n_tissues, 2 .* cm.n_p + 1, 2 .* cm.n_tau + 1 ] );
            
            % thermal equilibrium:
            cm.m( 2, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) = cm.pd;  
            
            cm.b_occ = false( size( cm.m ) );
            
            % only center of configuration space is occupied, initially:
            cm.b_occ( :, :, cm.p0( 1 ), cm.p0( 2 ), cm.p0( 3 ), cm.tau0 ) = true;
            
            % accumulated gradient moments
            
            p_max = cm.n_p( : ) .* cm.d_p( : );
            
            for i = 1 : 3

                cm.p{ i } = reshape( ...
                    linspace( - p_max( i ), p_max( i ), 2 * cm.n_p( i ) + 1 ), ...
                    [ ones( 1, 1 + i ), 2 * cm.n_p( i ) + 1 ] );
                
            end
            
            % accumulated durations
            
            tau_max = cm.n_tau * cm.d_tau;
           
            cm.tau = reshape( ...
                linspace( - tau_max, tau_max, 2 * cm.n_tau + 1 ), ...
                [ ones( 1, 5 ), 2 * cm.n_tau + 1 ] );
                                               
            % everything is prepared now
            
            cm.b_prep = true;
            
        end        
              
    end
    
end