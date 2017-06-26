function [ node_num, node_dim, element_num, element_order ] = ...
  gmsh_size_read ( gmsh_filename )

%*****************************************************************************80
%
%% GMSH_SIZE_READ reads sizes from a GMSH file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 October 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string GMSH_FILENAME, the GMSH filename.
%
%    Output, integer NODE_NUM, the number of nodes.
%
%    Output, integer NODE_DIM, the spatial dimension.
%
%    Output, integer ELEMENT_NUM, the number of elements.
%
%    Output, integer ELEMENT_ORDER, the order of the elements.
%
  r8_big = 1.0E+30;

  node_num = 0;
  node_dim = 0;
  element_num = 0;
  element_order = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;
%
%  Open the file.
%
  input = fopen ( gmsh_filename, 'rt' );

  if ( input < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_SIZE_READ - Error!\n' );
    fprintf ( 1, '  Could not open the input file "%s".\n', gmsh_filename );
    error ( 'GMSH_SIZE_READ - Error!' );
    return
  end

  level = 0;

  while ( 1 )

    text = fgetl ( input );

    if ( text == -1 )
      break
    end

    if ( level == 0 )
      if ( s_begin ( text(1:6), '$Nodes' ) )
        level = 1;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      node_num = value(1);
      level = 2;
    elseif ( level == 2 )
      if ( s_begin ( text(1:9), '$EndNodes' ) )
        break
      else
        [ value, count ] = sscanf ( text, '%d  %f  %f  %f' );
        indx = value(1);
        x = value(2);
        x_min = min ( x_min, x );
        x_max = max ( x_max, x );
        if ( 2 < count )
          y = value(3);
          y_min = min ( y_min, y );
          y_max = max ( y_max, y );
          if ( 3 < count )
            z = value(4);
            z_min = min ( z_min, z );
            z_max = max ( z_max, z );
          end
        end
      end
    end

  end
%
%  Make a very simple guess as to the dimensionality of the data.
%
  node_dim = 3;
  if ( z_max == z_min )
    node_dim = 2;
    if ( y_max == y_min )
      node_dim = 1;
    end
  end
%
%  Now read element information.
%
  level = 0;

  while ( 1 )

    text = fgetl ( input );

    if ( text == -1 )
      fprintf ( 'ran out\n' );
      break
    end

    if ( level == 0 )
      if ( s_begin ( text(1:9), '$Elements' ) )
        level = 1;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      element_num = value(1);
      level = 2;
    elseif ( level == 2 )
      if ( s_begin ( text(1:12), '$EndElements' ) )
        break
      else
        [ value, count ] = sscanf ( text, '%d' );
        element_order = count - 5;
        break
      end
    end

  end

  fclose ( input );

  return
end
