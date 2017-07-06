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
