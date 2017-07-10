function gmsh_mesh3d_write ( gmsh_filename, m, node_num, node_x, ...
  element_order, element_num, element_node )

%*****************************************************************************80
%
%% GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.
%
%  Discussion:
%
%    The node ordering for the 20 node element is not standard.
%
%    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:
%
%    1:    a
%    2:        b
%    3:            c
%    4:                d
%    5: (2*a  +b        )/3
%    6: (  a+2*b        )/3
%    7: (    2*b+  c    )/3
%    8: (      b+2*c    )/3
%    9: (  a    +2*c    )/3
%   10: (2*a    +  c    )/3
%   11: (2*a        +  d)/3
%   12: (  a        +2*d)/3
%   13: (     b     +2*d)/3
%   14: (   2*b     +  d)/3
%   15: (       +  c+2*d)/3
%   16: (       +2*c+  d)/3
%   17: (  a+  b+  c    )/3
%   18: (  a+  b    +  d)/3
%   19: (      b+  c+  d)/3
%   20: (  a+      c+  d)/3
%
%    Leo Rebholz used the following ordering:
%
%    1:    a
%    2:        b
%    3:            c
%    4:                d
%    5: (2*a  +b        )/3
%    6: (2*a    +  c    )/3
%    7: (  a+2*b        )/3
%    8: (  a    +2*c    )/3
%    9: (  a+  b+  c    )/3
%   10: (    2*b+  c    )/3
%   11: (      b+2*c    )/3
%   12: (2*a        +  d)/3
%   13: (   2*b     +  d)/3
%   14: (       +2*c+  d)/3
%   15: (  a+  b    +  d)/3
%   16: (      b+  c+  d)/3
%   17: (  a+      c+  d)/3
%   18: (  a        +2*d)/3
%   19: (     b     +2*d)/3
%   20: (       +  c+2*d)/3
%
%    Since the only 20 node data we have is from Leo, we will assume that
%    all 20 node input data is in Leo's format, and needs to be converted
%    to the Gmsh convention.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 October 2014
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Christophe Geuzaine, Jean-Francois Remacle,
%    Gmsh: a three-dimensional finite element mesh generator with
%    built-in pre- and post-processing facilities,
%    International Journal for Numerical Methods in Engineering,
%    Volume 79, Number 11, pages 1309-1331, 2009.
%
%  Parameters:
%
%    Input, string GMSH_FILENAME, the name of the Gmsh file.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, real NODE_X(M,NODE_NUM), the node coordinates.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the nodes
%    that make up each element.
%
  leo_to_gmsh = [...
     1,  2,  3,  4,  5, ...
     7, 10, 11,  8,  6, ...
    12, 18, 19, 13, 20, ...
    14,  9, 15, 16, 17 ];
%
%  Enforce 1-based indexing.
%
  element_node = mesh_base_one ( node_num, element_order, element_num, ...
    element_node );
%
%  Open the file.
%
  gmsh = fopen ( gmsh_filename, 'wt' );

  if ( gmsh < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_MESH3D_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'GMSH_MESH3D_WRITE - Error!' );
  end
%
%  Write the data.
%
  fprintf ( gmsh, '$MeshFormat\n' );
  fprintf ( gmsh, '2.2 0 8\n' );
  fprintf ( gmsh, '$EndMeshFormat\n' );

  fprintf ( gmsh, '$Nodes\n' );
  fprintf ( gmsh, '%d\n', node_num );
  for node = 1 : node_num
    fprintf ( gmsh, '%d', node );
    for dim = 1 : 3
      if ( dim <= m )
        fprintf ( gmsh, '  %g', node_x(dim,node) );
      else
        fprintf ( gmsh, '  %g', 0.0 );
      end
    end
    fprintf ( gmsh, '\n' );
  end
  fprintf ( gmsh, '$EndNodes\n' );
%
%  These are the Gmsh codes for 4, 10 and 20 node tetrahedral elements.
%
  if ( element_order == 4 )
    element_type = 4;
  elseif ( element_order == 10 )
    element_type = 11;
  elseif ( element_order == 20 )
    element_type = 29;
  end

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, '$Elements\n' );
  fprintf ( gmsh, '%d\n', element_num );
  for element = 1 : element_num
    fprintf ( gmsh, '%d  %d  %d  %d  %d', ...
      element, element_type, tag_num, tag1, element );
    for vertex = 1 : element_order
      if ( element_order == 20 )
        v = leo_to_gmsh(vertex);
      else
        v = vertex;
      end
      fprintf ( gmsh, '  %d', element_node(v,element) );
    end
    fprintf ( gmsh, '\n' );
  end
  fprintf ( gmsh, '$EndElements\n' );

  fclose ( gmsh );

  return
end
