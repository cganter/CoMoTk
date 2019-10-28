function [ out, sel ] = sfv ( in, varargin )
%SFV set multiple field values interactively

out = in;

opt = [];
str = [];

nVarargs = length( varargin );

if ( nVarargs > 0 && ~isempty( varargin{ 1 } ) )
    
    opt = varargin{ 1 };
    fn = fieldnames( opt );

else
    
    fn = fieldnames( out );
        
end

n_fields = length( fn );

len_fn = numel( int2str( n_fields ) );

max_fn = 0;

for i_ = 1 : n_fields
    
    if ( length( fn{ i_ } ) > max_fn )
        
        max_fn = length( fn{ i_ } );
        
    end
    
end

fs_1 = [ '%', int2str( len_fn ), 'd: %', int2str( max_fn ), 's  =  %' ];
fs_2 = 's\t %s\n';

if ( nVarargs == 2 )
    
    str = varargin{ 2 };
        
end

sel = 1;

while( sel > 0 )
    
    select_field( );
    
    if ( sel > 0 )
        
        set_field( );
        
    end
    
end

% -------------------
% auxiliary functions
% -------------------

    function select_field ( )
        
        sel = -2;
        
        while( ( sel < 1 || sel > n_fields ) && sel ~= 0 && sel ~= -1 )
            
            fprintf( 1, 'Select field to change (ok == 0, quit == -1)\n\n' );
            
            max_val = 0;
            
            for i = 1 : n_fields
                
                if ( isnumeric( out.( fn{ i } ) ) )
                    
                    s_ = num2str( out.( fn{ i } ) );
                    
                else
                                        
                    s_ = out.( fn{ i } );
                    
                end
                
                if ( length( s_ ) > max_val )
                    
                    max_val = length( s_ );
                    
                end
                
            end
            
            fs = [ fs_1, int2str( max_val ), fs_2 ];
            
            for i = 1 : n_fields
                
                s_ = '';
                
                if ( ~isempty( str ) && isfield( str, fn{ i } ) )
                    
                    s_ = str.( fn{ i } );
                    
                end

                if ( isnumeric( out.( fn{ i } ) ) )
                    
                    fprintf( 1, fs, i, fn{ i }, num2str( out.( fn{ i } ) ), s_ );
                    
                else
                                        
                    fprintf( 1, fs, i, fn{ i }, out.( fn{ i } ), s_ );
                    
                end
                
            end
            
            fprintf( 1, '\n' );
            
            sel = input( 'field = ' );
            
        end
        
    end

    function set_field ( )
        
        if ( isempty( opt ) || isempty( opt.( fn{ sel } ) ) )
            
            out.( fn{ sel } ) = input( [ fn{ sel }, ' = ' ] );
            
        else
          
            opt_ = opt.( fn{ sel } );
            n_opt = length( opt_ );

            c_opt = 0;
        
            while( c_opt < 1 || c_opt > n_opt )
                
                fprintf( 1, 'Select option\n\n' );
                
                for i = 1 : n_opt
                    
                    fprintf( 1, '%d: %s\n', i, opt_{ i } );
                    
                end
                
                fprintf( 1, '\n' );
                
                c_opt = input( 'option = ' );
                
            end
            
            out.( fn{ sel } ) = opt_{ c_opt };
            
        end
        
    end

end
