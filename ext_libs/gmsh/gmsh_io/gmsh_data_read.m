function [ node_x, element_node ] = gmsh_data_read ( gmsh_filename, ...
  node_dim, node_num, element_order, element_num )

%*****************************************************************************80
%
%% GMSH_DATA_READ reads data from a GMSH file.
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
%    Input, integer NODE_DIM, the spatial dimension.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Output, real NODE_X(NODE_DIM,NODE_NUM), the node coordinates.
%
%    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
%    the nodes that make up each element.
%
  node_x = zeros ( node_dim, node_num );
  element_node = zeros ( element_order, element_num );
%
%  Open the file.
%
  input = fopen ( gmsh_filename, 'rt' );

  if ( input < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_DATA_READ - Error!\n' );
    fprintf ( 1, '  Could not open the input file "%s".\n', gmsh_filename );
    error ( 'GMSH_DATA_READ - Error!' );
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
        j = 0;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      level = 2;
    elseif ( level == 2 )
      if ( s_begin ( text(1:9), '$EndNodes' ) )
        break
      else
        j = j + 1;
        [ value, count ] = sscanf ( text, '%d  %f  %f  %f' );
        indx = value(1);
        node_x(1,j) = value(2);
        if ( 2 < count )
          node_x(2,j) = value(3);
          if ( 3 < count )
            node_x(3,j) = value(4);
          end
        end
      end
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
        j = 0;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      level = 2;
    elseif ( level == 2 )
      if ( s_begin ( text(1:12), '$EndElements' ) )
        break
      else
        j = j + 1;
        [ value, count ] = sscanf ( text, '%d' );
        for i = 1 : element_order
          element_node(i,j) = value(5+i);
        end
      end
    end

  end

  fclose ( input );

  return
end
