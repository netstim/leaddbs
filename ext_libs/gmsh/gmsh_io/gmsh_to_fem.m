function gmsh_to_fem ( prefix )

%*****************************************************************************80
%
%% MAIN is the main program for GMSH_TO_FEM.
%
%  Discussion:
%
%    GMSH_TO_FEM converts mesh data from GMSH to FEM format.
%
%  Usage:
%
%    gmsh_to_fem prefix
%
%    where 'prefix' is the common filename prefix:
%
%    * 'prefix'.msh contains the GMSH mesh data.
%    * 'prefix'_nodes.txt will contain the node coordinates.
%    * 'prefix'_elements.txt will contain the element node connectivity.
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
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'GMSH_TO_FEM\n' );
  fprintf ( 1, '  MATLAB version:\n' );
  fprintf ( 1, '  Read a mesh description created by GMSH:\n' );
  fprintf ( 1, '  * "prefix".msh, contains the GMSH mesh data.\n' );
  fprintf ( 1, '  Write two simple FEM files listing nodes and elements.\n' );
  fprintf ( 1, '  * "prefix"_nodes.txt, node coordinates.\n' );
  fprintf ( 1, '  * "prefix"_elements.txt, element connectivity.\n' );
%
%  Get the filename prefix.
%
  if ( nargin < 1 )

    prefix = input ( 'Enter the filename prefix:  ' );

  end
%
%  Create the filenames.
%
  gmsh_filename = strcat ( prefix, '.msh' );
  fem_node_filename = strcat ( prefix, '_nodes.txt' );
  fem_element_filename = strcat ( prefix, '_elements.txt' );
%
%  Read the GMSH sizes.
%
  [ node_num, m, element_num, element_order ] = gmsh_size_read ( ...
    gmsh_filename );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Size information from GMSH files:\n' );
  fprintf ( 1, '  Spatial dimension M = %d\n', m );
  fprintf ( 1, '  Number of nodes NODE_NUM = %d\n', node_num );
  fprintf ( 1, '  Number of elements ELEMENT_NUM = %d\n', element_num );
  fprintf ( 1, '  Element order ELEMENT_ORDER = %d\n', element_order );
%
%  Read GMSH data.
%
  [ node_x, element_node ] = gmsh_data_read ( gmsh_filename, m, node_num, ...
    element_order, element_num );
%
%  Write FEM data.
%
  r8mat_write ( fem_node_filename, m, node_num, node_x );

  i4mat_write ( fem_element_filename, element_order, element_num, ...
    element_node );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'GMSH_TO_FEM:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
function c = ch_cap ( c )

%*****************************************************************************80
%
%% CH_CAP capitalizes a single character.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 April 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character C, the character to capitalize.
%
%    Output, character C, the capitalized character.
%
  if ( 'a' <= c && c <= 'z' )
    c = c + 'A' - 'a';
  end

  c = char ( c );

  return
end
function truefalse = ch_eqi ( c1, c2 )

%*****************************************************************************80
%
%% CH_EQI is a case insensitive comparison of two characters for equality.
%
%  Example:
%
%    CH_EQI ( 'A', 'a' ) is TRUE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 July 2000
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character C1, C2, the characters to compare.
%
%    Output, logical TRUEFALSE, is TRUE (1) if the characters are equal.
%
  if ( ch_cap ( c1 ) == ch_cap ( c2 ) )
    truefalse = 1;
  else
    truefalse = 0;
  end

  return
end
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
function i4mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% I4MAT_WRITE writes an I4MAT file.
%
%  Discussion:
%
%    An I4MAT is an array of I4's.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, integer TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'I4MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %d', round ( table(i,j) ) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% R8MAT_WRITE writes an R8MAT file.
%
%  Discussion:
%
%    An R8MAT is an array of R8's.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  Alternative print statements include:
%
%     fprintf ( output_unit, '  %24.16e', table(i,j) );
%     fprintf ( output_unit, '  %14.6e', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %g', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function value = s_begin ( s1, s2 )

%*****************************************************************************80
%
%% S_BEGIN is TRUE if one string matches the beginning of the other.
%
%  Discussion:
%
%    The strings are compared, ignoring blanks and capitalization.
%
%  Example:
%
%     S1              S2      S_BEGIN
%
%    'Bob'          'BOB'     TRUE
%    '  B  o b '    ' bo b'   TRUE
%    'Bob'          'Bobby'   TRUE
%    'Bobo'         'Bobb'    FALSE
%    ' '            'Bob'     FALSE    (Do not allow a blank to match
%                                       anything but another blank string.)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 January 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character S1(*), S2(*), the strings to be compared.
%
%    Output, logical S_BEGIN, is TRUE if the strings match up to
%    the end of the shorter string, ignoring case.
%
  len1 = s_len_trim ( s1 );
  len2 = s_len_trim ( s2 );
%
%  If either string is blank, then both must be blank to match.
%  Otherwise, a blank string matches anything, which is not
%  what most people want.
%
  if ( len1 == 0 || len2 == 0 )

    if ( len1 == 0 && len2 == 0 )
      value = 1;
    else
      value = 0;
    end

    return

  end

  i1 = 0;
  i2 = 0;
%
%  Find the next nonblank in S1.
%
  while ( 1 )

    while ( 1 )

      i1 = i1 + 1;

      if ( len1 < i1 )
        value = 1;
        return
      end

      if ( s1(i1) ~= ' ' )
        break
      end

    end
%
%  Find the next nonblank in S2.
%
    while ( 1 )

      i2 = i2 + 1;

      if ( len2 < i2 )
        value = 1;
        return
      end

      if ( s2(i2) ~= ' ' )
        break
      end

    end
%
%  If the characters match, get the next pair.
%
    if ( ~ch_eqi ( s1(i1), s2(i2) ) )
      break
    end

  end

  value = 0;

  return
end
function len = s_len_trim ( s )

%*****************************************************************************80
%
%% S_LEN_TRIM returns the length of a character string to the last nonblank.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 June 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string S, the string to be measured.
%
%    Output, integer LEN, the length of the string up to the last nonblank.
%
  len = length ( s );

  while ( 0 < len )
    if ( s(len) ~= ' ' )
      return
    end
    len = len - 1;
  end

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
