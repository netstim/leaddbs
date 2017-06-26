function [ element_num, element_order ] = gmsh_mesh2d_element_size_example ( )

%*****************************************************************************80
%
%% GMSH_MESH2D_ELEMENT_SIZE_EXAMPLE: element size information for the example.
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
%    Output, integer ELEMENT_NUM, the number of elements.
%
%    Output, integer ELEMENT_ORDER, the order of the elements.
%
  element_num = 24;
  element_order = 3;

  return
end
