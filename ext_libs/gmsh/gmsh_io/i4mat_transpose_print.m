function i4mat_transpose_print ( m, n, a, title )

%*****************************************************************************80
%
%% I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 September 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns.
%
%    Input, integer A(M,N), an M by N matrix to be printed.
%
%    Input, string TITLE, a title.
%
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return
end
