function [ component_num, c ] = i4block_components ( l, m, n, a )

%*****************************************************************************80
%
%% I4BLOCK_COMPONENTS assigns contiguous nonzero pixels to a common component.
%
%  Discussion:
%
%    On input, the A array contains values of 0 or 1.
%
%    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
%    into connected components.
%
%    The pixel A(I,J,K) is "connected" to the pixels:
%
%      A(I-1,J,  K  ),  A(I+1,J,  K  ),
%      A(I,  J-1,K  ),  A(I,  J+1,K  ),
%      A(I,  J,  K-1),  A(I,  J,  K+1),
%
%    so most pixels have 6 neighbors.
%
%    On output, COMPONENT_NUM reports the number of components of nonzero
%    data, and the array C contains the component assignment for
%    each nonzero pixel, and is 0 for zero pixels.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 February 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer L, M, N, the order of the array.
%
%    Input, integer A(L,M,N), the pixel array.
%
%    Output, integer COMPONENT_NUM, the number of components
%    of nonzero data.
%
%    Output, integer C(L,M,N), the component array.
%

%
%  Initialization.
%
  c = zeros ( l, m, n );
  component_num = 0;
%
%  P is simply used to store the component labels.  The dimension used
%  here is, of course, usually an absurd overestimate.
%
  p = 1 : l * m * n;
%
%  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
%  west neighbor already has a label, the current pixel inherits it.
%  In case the labels disagree, we need to adjust the P array so we can
%  later deal with the fact that the two labels need to be merged.
%
  for i = 1 : l

    for j = 1 : m

      for k = 1 : n

        if ( i == 1 )
          north = 0;
        else
          north = c(i-1,j,k);
        end

        if ( j == 1 )
          west = 0;
        else
          west = c(i,j-1,k);
        end

        if ( k == 1 )
          up = 0;
        else
          up = c(i,j,k-1);
        end

        if ( a(i,j,k) ~= 0 )
%
%  New component?
%
          if ( north == 0 && west == 0 && up == 0 )
            component_num = component_num + 1;
            c(i,j,k) = component_num;
%
%  One predecessor is labeled.
%
          elseif ( north ~= 0 && west == 0 && up == 0 )
            c(i,j,k) = north;
          elseif ( north == 0 && west ~= 0 && up == 0 )
            c(i,j,k) = west;
          elseif ( north == 0 && west == 0 && up ~= 0 )
            c(i,j,k) = up;
%
%  Two predecessors are labeled.
%
          elseif ( north == 0 && west ~= 0 && up ~= 0 )
            c(i,j,k) = min ( west, up );
            c1 = min ( p(west), p(up) );
            p(west) = c1;
            p(up) = c1;
          elseif ( north ~= 0 && west == 0 && up ~= 0 )
            c(i,j,k) = min ( north, up );
            c1 = min ( p(north), p(up) );
            p(north) = c1;
            p(up) = c1;
          elseif ( north ~= 0 && west ~= 0 && up == 0 )
            c(i,j,k) = min ( north, west );
            c1 = min ( p(north), p(west) );
            p(north) = c1;
            p(west) = c1;
%
%  Three predecessors are labeled.
%
          elseif ( north ~= 0 && west ~= 0 && up ~= 0 )
            c(i,j,k) = min ( north, min ( west, up ) );
            c1 = min ( p(north), min ( p(west), p(up) ) );
            p(north) = c1;
            p(west) = c1;
            p(up) = c1;
          end

        end

      end

    end

  end
%
%  When a component has multiple labels, have the higher labels
%  point to the lowest one.
%
  for component = component_num : -1 : 1
    b = component;
    while ( p(b) ~= b )
      b = p(b);
    end
    p(component) = b;
  end
%
%  Locate the minimum label for each component.
%  Assign these mininum labels new consecutive indices.
%
  q = zeros ( 1, component_num );
  i = 0;
  for component = 1 : component_num
    if ( p(component) == component )
      i = i + 1;
      q(component) = i;
    end
  end

  component_num = i;
%
%  Replace the labels by consecutive labels.
%
  i = find ( c ~= 0 );
  c(i) = q ( p ( c(i) ) );

  return
end
