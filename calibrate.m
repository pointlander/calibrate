# Copyright 2020 The Calibrate Authors. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

# compute the calibration values for an accelerometer

clear all;

global a = dlmread( "accelerometer.csv" );

function parameters = find_parameters( current )
   global a;
   parameters = zeros( 6, 5 );
   [ rows, columns ] = size( a );
   for c = 1 : 6
      start = current( c, 4 ) + 1 - floor( 3 * rand() );
      stop = current( c, 5 ) + 1 - floor( 3 * rand() );
      if ( start < 1 )
         start = 1;
      endif
      if ( stop > rows )
         stop = rows
      endif
      if ( stop < start )
         stop = start;
      endif
      count = 0;
      for d = start : stop
         parameters( c, 1:3 ) += a( d, : );
         count++;
      endfor
      parameters( c, 1:3 ) /= count;
      parameters( c, 4 ) = start;
      parameters( c, 5 ) = stop;
   endfor
endfunction

global G = 9.80665^2;
function sum = squaredDifferenceSum( calibration )
   global a;
   global G;
   [ rows, columns ] = size( a );
   sum = 0;
   for c = 1 : rows
      sum = sum + ( ( ( a( c, 1 ) - calibration( 1 ) ) / calibration( 2 ) )^2 +
                    ( ( a( c, 2 ) - calibration( 3 ) ) / calibration( 4 ) )^2 +
                    ( ( a( c, 3 ) - calibration( 5 ) ) / calibration( 6 ) )^2 - G );
   endfor
   sum = abs( sum / rows );
endfunction

global P;
function y = f( x )
   global G;
   global P;
   y = zeros( 1, 6 );
   for i = 1 : 6
      y( i ) = ( ( P( i, 1 ) - x( 1 ) ) / x( 2 ) )^2 + ...
               ( ( P( i, 2 ) - x( 3 ) ) / x( 4 ) )^2 + ...
               ( ( P( i, 3 ) - x( 5 ) ) / x( 6 ) )^2 - G;
   endfor
endfunction

current_calibration = [ 0, 1, 0, 1, 0, 1 ];
current_sum = squaredDifferenceSum( current_calibration );
current = zeros( 6, 5 );
d = 1;
for c = 1 : 3
   [ m, i ] = max( a( :, c ) );
   current( d, 4 ) = current( d, 5 ) = i;
   d++;

   [ m, i ] = min( a( :, c ) );
   current( d, 4 ) = current( d, 5 ) = i;
   d++;
endfor

best_sum = 1000;
T = 400;
while ( T > 0 )
   P = find_parameters( current );
   [ calibration, information, message ] = fsolve( "f", current_calibration );
   sum = squaredDifferenceSum( calibration );
   if ( sum < best_sum )
      best_sum = sum;
      best_calibration = calibration;
      best = P;
   endif
   acceptance = 1;
   if ( sum > current_sum )
      acceptance = exp( ( current_sum - sum ) / T );
   endif
   if ( acceptance > rand() )
      current = P;
      current_calibration = calibration;
      current_sum = sum;
   endif
   T--;
endwhile

best
best_calibration
best_sum
